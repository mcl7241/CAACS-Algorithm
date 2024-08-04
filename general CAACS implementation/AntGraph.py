from threading import Lock
import numpy as np
from sklearn.manifold import MDS
import matplotlib.pyplot as plt
import os

class AntGraph:
    def __init__(self, instance_name, num_ants, num_nodes, delta_mat, carbon_mat, scaled_mat, clusters_mat, tau_mat=None):
        #print (len(delta_mat))
        if len(delta_mat) != num_nodes:
            raise Exception("len(delta) != num_nodes")

        self.instance_name = instance_name
        self.num_ants = num_ants
        self.num_nodes = num_nodes
        self.delta_mat = delta_mat # matrix of node distance deltas
        self.carbon_mat = carbon_mat
        self.scaled_mat = scaled_mat
        self.clusters_mat = clusters_mat
        self.lock = Lock()

        # symmetrizing the matrix by taking the average of the matrix and its transpose
        # can alternatively be done with checking if a matrix is symmetric
        temp_delta_mat = np.array(delta_mat)
        symmetric_distances = (temp_delta_mat + temp_delta_mat.T)/2

        # applying MDS
        mds = MDS(n_components=2, dissimilarity="precomputed", random_state=6)
        positions = mds.fit_transform(symmetric_distances)

        self.positions = []
        for i in range(len(positions)):
            self.positions.append((positions[i, 0], positions[i, 1]))

        # tau mat contains the amount of phermone at node x,y
        if tau_mat is None:
            self.tau_mat = []
            for i in range(0, num_nodes):
                self.tau_mat.append([0]*num_nodes)

    def delta(self, r, s):
        return self.delta_mat[r][s]

    def carbon(self, r, s):
        return self.carbon_mat[r][s]
    
    def scaled_emission(self, r, s):
        return self.scaled_mat[r][s]

    def tau(self, r, s):
        return self.tau_mat[r][s]

    # 1 / delta = eta or etha 
    def etha(self, r, s):
        if self.delta(r, s) == 0.0:
            return float('inf')
        return 1.0 / self.delta(r, s)

    # inner locks most likely not necessary
    def update_tau(self, r, s, val):
        lock = Lock()
        lock.acquire()
        self.tau_mat[r][s] = val
        lock.release()

    def reset_tau(self):
        lock = Lock()
        lock.acquire()
        avg = self.average_delta()

        # initial tau 
        self.tau0 = 1.0 / (self.num_nodes * 0.5 * avg)

        #print ("Average = %s" % (avg,))
        #print ("Tau0 = %s" % (self.tau0))

        for r in range(0, self.num_nodes):
            for s in range(0, self.num_nodes):
                self.tau_mat[r][s] = self.tau0
        lock.release()

    # average delta in delta matrix
    def average_delta(self):
        return self.average(self.delta_mat)

    # average tau in tau matrix
    def average_tau(self):
        return self.average(self.tau_mat)

    # average val of a matrix
    def average(self, matrix):
        sum = 0
        for r in range(0, self.num_nodes):
            for s in range(0, self.num_nodes):
                sum += matrix[r][s]

        avg = sum / (self.num_nodes * self.num_nodes)
        return avg
    
    # get the cluster id of a given node
    def get_cluster(self, s):
        for i in range(len(self.clusters_mat)):
            cluster = self.clusters_mat[i]
            for j in range(len(cluster)):
                if cluster[j] == s:
                    return i
    
    def carbon_cost(self, best_path_vec, carbon_mat):
        carbon_total = 0
        for i in range(len(best_path_vec)):
            if i+1 != len(best_path_vec):
                city_one = best_path_vec[i]
                city_two = best_path_vec[i+1]
                carbon_total += carbon_mat[city_one][city_two]
            else: 
                city_one = best_path_vec[i]
                city_two = best_path_vec[0]
                carbon_total += carbon_mat[city_one][city_two]
        return carbon_total
    
    # create image
    def create_image(self, iteration, best_path_vec, best_path_cost):
        # hard-coded colors, can change later
        # (R, G, B, Opacity)
        colors = ['blue', 'orange', 'purple', 'gray', 'yellow', 
                  'cyan', 'magenta', 'teal', 'pink', 'lime', 'black',
                  (1, 0, 0, 0.5), (0, 1, 0, 0.5), (0, 0, 1, 0.5),
                  (0.5, 0.5, 0, 0.5), (0, 0.5, 0.5, 0.5), (0.5, 0, 0.5, 0.5), 
                  (0.3, 0.3, 0.3, 0.5), (0.7, 0.2, 0.8, 0.6), (0.4, 0.7, 0.1, 0.8), 
                  (0.5, 0, 0, 0.5), (0, 0.5, 0, 0.5), (0, 0, 0.5, 0.5), 
                  (0.1, 0.2, 0.3), (0.1, 0.3, 0.2)]

        # create 2 by 2 size plot
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.set_facecolor('white')
        
        # plot points in each cluster with different colors
        for i in range(len(self.clusters_mat)):
            color = 'black'
            if i <= len(colors)-1:
                color = colors[i]
            #(j, val) represent the index and val of a node within one cluster 
            for j, val in enumerate(self.clusters_mat[i]):
                x, y = self.positions[val]
                ax.plot(x, y, marker=(4, 0, 90), markersize=3, linestyle='None', 
                        color=color, markeredgecolor='black', markeredgewidth=0.3)
                # maybe add ax.circle around the clusters
        # set aspect of plot to be equal
        ax.set_aspect('equal', adjustable='datalim')

        # turn off the axes
        ax.axis('off')

        if (iteration != 'None'):
            path = []
            for i, val in enumerate(best_path_vec):
                x, y = self.positions[val]
                path.append((x, y))
                plt.annotate(str(i), (x, y), textcoords="offset points", xytext=(2, 2), fontsize=4)
            # append the source node
            x, y = self.positions[best_path_vec[0]]
            path.append((x, y))

            path_x, path_y = zip(*path) # seperates the x and y coordinates
            if iteration == 'F':
                ax.plot(path_x, path_y, 'g-')
            else:
                ax.plot(path_x, path_y, 'r-') # red line
            
            '''
            ax.set_xlabel('')
            ax.set_ylabel('')
            '''
            total_carbon = 0
            total_cost = 0

            for i in range(len(best_path_vec)):
                r = best_path_vec[i]
                s = best_path_vec[0]
                if (i + 1) != len(best_path_vec):
                    s = best_path_vec[i + 1]
                total_carbon += self.carbon(r, s)
                total_cost += self.delta(r, s)
            #print("Carbon = ", total_carbon, " Cost = ", total_cost)

            ax.set_title('I:' + iteration + ',' + 'Cost:' + str(int(best_path_cost)) + ',' 
                        + 'Carbon:' + str(int(total_carbon)), fontsize=6)
        # save the plot to a file
        directory_path = './image/' + self.instance_name + '/'
        # Check if the directory exists
        if not os.path.exists(directory_path):
            # If the directory does not exist, create it
            os.makedirs(directory_path)
        plt.savefig(directory_path + 'pic_iter_' + iteration + "_num_ants_" + str(self.num_ants) +'.png', format='png', dpi=200, bbox_inches='tight', pad_inches=0)
        plt.close('all')


