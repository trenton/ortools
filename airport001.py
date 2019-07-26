import sys
import functools
from math import radians, cos, sin, asin, sqrt
import csv
from collections import namedtuple, defaultdict

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

from vincenty import vincenty

from haversine import haversine, Unit

from geopy.distance import geodesic, great_circle


Airport = namedtuple('Airport', 'code lat_rad lon_rad county state lat_deg lon_deg')

all_states = (
    "AL AZ AR CA CO CT DE FL GA ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT "
    "NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY "
)
ALLOWED_STATES = set(all_states.split())
ALLOWED_STATES = set("ID MT NV UT WY ".split())

DEFAULT_AIRPORTS = [
    Airport("PAE", 172465.2000, 440213.7000, "KING", "WA", None, None),
    Airport("GEG", 171428.5000, 423126.8000, "SPOKANE", "WA", None, None),
    Airport("W56", 164474.2080, 441078.6710, "CLARK", "WA", None, None),
    Airport("VUO", 164233.6290, 441563.3580, "CLARK", "WA", None, None),
    Airport("2S0", 174062.1070, 432338.2230, "OKANOGAN", "WA", None, None),
]

AIRPORTS = []
DISJUNCTIONS = {}

ALLOW_ANY_START = True


def print_solution(manager, routing, assignment):
    """Prints assignment on console."""
    print('distance: {} miles'.format(assignment.ObjectiveValue()))
    index = routing.Start(0)
    plan_output = 'Route:\n'
    route_distance = 0
    while not routing.IsEnd(index):
        city_name = AIRPORTS[manager.IndexToNode(index)].code
        #plan_output += ' {} ->'.format(manager.IndexToNode(index))
        plan_output += f' {city_name} '
        previous_index = index
        index = assignment.Value(routing.NextVar(index))
        leg_distance = routing.GetArcCostForVehicle(previous_index, index, 0)
        #plan_output += f" - {leg_distance} - "
        route_distance += leg_distance
    plan_output += ' {}\n'.format(AIRPORTS[manager.IndexToNode(index)].code)
    print(plan_output)
    plan_output += 'Route distance: {}miles\n'.format(route_distance)

@functools.lru_cache(maxsize=None)
def do_vincenty(lon1, lat1, lon2, lat2):
    return vincenty((lat1, lon1), (lat2, lon2), miles=True) / 1.151

@functools.lru_cache(maxsize=None)
def do_great_circle(lon1, lat1, lon2, lat2):
    return great_circle((lat1, lon1), (lat2, lon2)).miles / 1.151

@functools.lru_cache(maxsize=None)
def do_haversine2(lon1, lat1, lon2, lat2):
    return haversine((lat1, lon1), (lat2, lon2), unit=Unit.NAUTICAL_MILES)

@functools.lru_cache(maxsize=None)
def do_haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)

    https://stackoverflow.com/questions/4913349

    Requires args in radians. The library one assumes degrees
    """
    if lon1 == lon2 and lat1 == lat2:
        return 0

    # convert decimal degrees to radians 
    #lon1, lat1, lon2, lat2 = map(radians, [lon1/3600, lat1/3600, lon2/3600, lat2/3600])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 3440 # radius of earth in chosen units
    return c * r


def main(depot=0):
    """Entry point of the program."""

    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(AIRPORTS), 1, depot)

    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)


    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)

        from_airport = AIRPORTS[from_node]
        to_airport = AIRPORTS[to_node]

        if all([from_airport.lon_rad, from_airport.lat_rad, to_airport.lon_rad, to_airport.lat_rad]):
            # :TODO: more accurate, but waaaay slower
            #return do_great_circle(from_airport.lon_deg, from_airport.lat_deg, to_airport.lon_deg, to_airport.lat_deg)
            #return do_vincenty(from_airport.lon_deg, from_airport.lat_deg, to_airport.lon_deg, to_airport.lat_deg)
            #return do_haversine2(from_airport.lon_deg, from_airport.lat_deg, to_airport.lon_deg, to_airport.lat_deg)
            return do_haversine(from_airport.lon_rad, from_airport.lat_rad, to_airport.lon_rad, to_airport.lat_rad)
        else:
            return 0

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Put a huge penalty in for all disjunction members
    for k, v in DISJUNCTIONS.items():
        routing.AddDisjunction(v, 1000000000)

    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

    # Solve the problem.
    assignment = routing.SolveWithParameters(search_parameters)

    # Print solution on console.
    if assignment:
        print_solution(manager, routing, assignment)


def load_airports_tsv(airport_data_file):
    with open(airport_data_file, 'r') as f:
        airport_data = csv.DictReader(f, delimiter='\t')

        airports = []
        for row in airport_data:
            if row['Type'] != 'AIRPORT' or row['Use'] != 'PU':
                continue
            
            if row['State'] not in ALLOWED_STATES:
                continue

            #if not row['IcaoIdentifier']:
            #    continue

            #if "100" not in row['FuelTypes']:
            #    continue

            # prefer ICAO airport code
            if row['IcaoIdentifier']:
                code = row['IcaoIdentifier']
            else:
                code = row['LocationID'].strip("'")
                
            ns_hemi = 1 if row['ARPLatitudeS'].endswith("S") else -1
            ew_hemi = 1 if row['ARPLongitudeS'].endswith("E") else -1

            lat = float(row['ARPLatitudeS'].strip("NS")) / 3600
            lat *= ns_hemi

            lon = float(row['ARPLongitudeS'].strip("EW")) / 3600
            lon *= ew_hemi

            airport = Airport(
		code,
		#radians(float(row['ARPLatitudeS'].strip("NS"))/3600) * ns_hemi,
		radians(lat),
		#radians(float(row['ARPLongitudeS'].strip("EW"))/3600) * ew_hemi,
		radians(lon),
		row['County'],
		row['State'],
                lat,
                lon,
            )

            airports.append(airport)

        return airports


if __name__ == '__main__':
    if len(sys.argv) >= 2:
        print(f"Loading airports from {sys.argv[1]}... ", end="")
        airports = load_airports_tsv(sys.argv[1])

        if len(sys.argv) >= 3:
            base_index = [x.code for x in airports].index(sys.argv[2])
        else:
            # insert a free dummy airport and always start from there
            dummy = Airport('DUMMY', None, None, None, None, None, None)
            airports.insert(0, dummy)
            base_index = 0

        disjunctions = defaultdict(list)
        for i, airport in enumerate(airports):
            #disjunctions[airport.county].append(i)
            disjunctions[airport.state].append(i)

        # setup globals
        AIRPORTS = airports
        DISJUNCTIONS = disjunctions
        #print(DISJUNCTIONS)
        print(f"loaded {len(AIRPORTS)} and {len(DISJUNCTIONS)} disjunctions")
        print(f"from: {','.join(ALLOWED_STATES)}")

        main(base_index)

    else:
        AIRPORTS = DEFAULT_AIRPORTS
        main()
