import sys
import math
import random
from AntColony_GTSP import AntColony
from AntGraph_GTSP import AntGraph

fuel_to_air = (1, 1, 1)
gravitational_constant = (9.81, 9.81, 9.81)
air_density = (1.2041, 1.2041, 1.2041)
coef_rolling = (0.01, 0.01, 0.01)
diesel_engine_efficiency = (0.45, 0.45, 0.45)
heating_value_diesel = (44, 44, 44)
vehicle_speed = (22.2, 22.2, 22.2)
conversion_factor = (737, 737, 737)
road_angle = (0, 0, 0)
kerb_weight = (3500, 5500, 13154)
max_payload = (4000, 12500, 17236)
engine_friction_factor = (0.25, 0.20, 0.15)
engine_speed = (38.3, 36.7, 30.2)
engine_displacement = (4.50, 6.90, 6.66)
aero_drag = (0.6, 0.7, 0.7)
frontal_SA = (7.0, 8.0, 9.8)
vehicle_drive_eff = (0.45, 0.45, 0.50)
u = 2.63
# supposedly nonconstant stuff
speed_ij = (60, 60, 60)
# weight carried by the car
F_ijkpt = (200, 300, 400)

def create_carbon_emission(distance): #0, 1, 2 for vehicle type
    # Generate a random integer from the set {0, 1, 2}
    vehicle_type = random.choice([0, 1, 2])
    
    # extra constant stuff
    lambda_s = fuel_to_air[vehicle_type]/(heating_value_diesel[vehicle_type]*conversion_factor[vehicle_type])
    s = gravitational_constant[vehicle_type] * (math.sin(road_angle[vehicle_type]) + coef_rolling[vehicle_type]*math.cos(road_angle[vehicle_type]))
    gamma_k = 1/(1000*diesel_engine_efficiency[vehicle_type]*vehicle_drive_eff[vehicle_type])
    beta_k = 0.5*aero_drag[vehicle_type]*air_density[vehicle_type]*frontal_SA[vehicle_type]
    y_k = engine_friction_factor[vehicle_type] * engine_speed[vehicle_type] * engine_displacement[vehicle_type]

    engine_module = (lambda_s * y_k * distance)/speed_ij[vehicle_type]
    speed_module = lambda_s*gamma_k*beta_k*distance*(speed_ij[vehicle_type]**2)
    weight_module = lambda_s * gamma_k * s * distance * (kerb_weight[vehicle_type] + F_ijkpt[vehicle_type])

    carbon_emission = (engine_module + speed_module + weight_module) * u

    return (carbon_emission) # remove int 

if __name__ == "__main__":   
    # help eventually establish individual subfolders for each instance
    file_path = '../GTSP_InstancesText/'
    file_name = '20kroA100.txt'

    f = open(file_path+file_name, 'r', encoding='iso-8859-1')

    nodes_temp = f.readline().strip().split()
    num_nodes = int(nodes_temp[1])

    clusters_temp = f.readline().strip().split()
    num_clusters = int(clusters_temp[1])

    # skip the symmetric amd triangle for now
    f.readline()
    f.readline()

    clusters_mat = []
    for i in range(num_clusters):
        clusters_list = f.readline().strip().split()
        clusters_num_list = []
        for j, val in enumerate(clusters_list):
            if j != 0:
                clusters_num_list.append(int(val)-1) # reindex 1-n to 0-(n-1)
        clusters_mat.append(clusters_num_list)
    print(clusters_mat)

    if num_nodes <= 10:
        num_ants = 20
        num_iterations = 12
        num_repetitions = 1
    else:
        num_ants = 40
        # 4 cluster -> 20
        # 8 cluster -> 20-30 iteration
        # for 16 -> 50-60 iteration
        # we should set a termination critera if the optimal path doesn't change after 15 iteration
        num_iterations = 100 
        num_repetitions = 1
    
    cities = [str(i) for i in range(num_nodes)]
    print(cities)

    cost_mat = []
    carbon_mat = []
    while True:
        line = f.readline().strip()
        if not line:
            break
        row = line.strip().split()
        number_row = []
        carbon_row = []
        for num_str in row:
            number = int(num_str)
            number_row.append(number)
            carbon_val = create_carbon_emission(number)
            carbon_row.append(carbon_val)
        cost_mat.append(number_row)
        carbon_mat.append(carbon_row)
    f.close()

    print(cost_mat)
    print()

    if num_nodes < len(cost_mat):
        cost_mat = cost_mat[0:num_nodes]
        for i in range(0, num_nodes):
            cost_mat[i] = cost_mat[i][0:num_nodes]

    print (cost_mat)

    if num_nodes < len(carbon_mat):
        carbon_mat = carbon_mat[0:num_nodes]
        for i in range(0, num_nodes):
            carbon_mat[i] = carbon_mat[i][0:num_nodes]
    
    print(carbon_mat)
    print()

    try:
        graph = AntGraph(file_name[:-4], num_ants, num_nodes, cost_mat, carbon_mat, clusters_mat)
        best_path_cost = sys.maxsize
        best_path_vec = None
        for i in range(0, num_repetitions):
            graph.reset_tau()
            ant_colony = AntColony(graph, num_ants, num_iterations)
            ant_colony.start()
            if ant_colony.best_path_cost < best_path_cost:
                best_path_cost = ant_colony.best_path_cost
                best_path_vec = ant_colony.best_path_vec
                c_cost = graph.carbon_cost(best_path_vec, carbon_mat)

        print ("\n************************************************************")
        print ("                   Final Results                              ")
        print ("**************************************************************")
        print ("\nBest path found = %s" % (best_path_vec,))
        for node in best_path_vec:
            print (cities[node] + " ",)
        print ("\nBest path cost = %s\n" % (best_path_cost,))
        print ("\nBest path carbon cost = %s\n" % (c_cost,))
        graph.create_image('final', best_path_vec, best_path_cost)
    
    except Exception as e:
        print ("exception: " + str(e))
