import math
import random
import sys
from threading import *

class Ant(Thread):
    def __init__(self, ID, start_node, colony):
        Thread.__init__(self)
        self.ID = ID
        self.start_node = start_node
        self.colony = colony

        self.curr_node = self.start_node
        self.graph = self.colony.graph
        self.path_vec = []
        self.path_vec.append(self.start_node)
        self.path_cost = 0

        # same meaning as in standard equations
        self.Beta = 1
        #self.Q0 = 1  # Q0 = 1 works just fine for 10 city case (no explore)
        self.Q0 = 0.5
        self.Rho = 0.99

        # store the nodes remaining to be explored here
        self.nodes_to_visit = {}

        for i in range(0, self.graph.num_nodes):
            if i != self.start_node:
                self.nodes_to_visit[i] = i

        # store the clusters remaining to be explored here
        self.clusters_to_visit = {}
        for i in range(len(self.graph.clusters_mat)):
            if i != self.graph.get_cluster(self.start_node):
                self.clusters_to_visit[i] = i

        # create n X n matrix 0'd out to start
        self.path_mat = []

        for i in range(0, self.graph.num_nodes):
            self.path_mat.append([0]*self.graph.num_nodes)

    # overide Thread's run()
    def run(self):
        graph = self.colony.graph
        while not self.end():
            # we need exclusive access to the graph
            graph.lock.acquire()
            new_node = self.state_transition_rule(self.curr_node)
            self.path_cost += graph.delta(self.curr_node, new_node)

            self.path_vec.append(new_node)
            self.path_mat[self.curr_node][new_node] = 1  #adjacency matrix representing path

            print ("Ant %s : %s, %s" % (self.ID, self.path_vec, self.path_cost,))
            
            self.local_updating_rule(self.curr_node, new_node)
            graph.lock.release()

            self.curr_node = new_node

        # don't forget to close the tour
        self.path_cost += graph.delta(self.path_vec[-1], self.path_vec[0])

        # send our results to the colony
        self.colony.update(self)
        print ("Ant thread %s terminating." % (self.ID,))

        # allows thread to be restarted (calls Thread.__init__)
        self.__init__(self.ID, self.start_node, self.colony)

    def end(self):
        return not self.clusters_to_visit 

    # described in report -- determines next node to visit after curr_node
    def state_transition_rule(self, curr_node):
        graph = self.colony.graph
        q = random.random()
        max_node = -1

        # sanity check
        if graph.get_cluster(curr_node) in self.clusters_to_visit.values():
            del self.clusters_to_visit[graph.get_cluster(curr_node)]

        if q < self.Q0:
            # missing alpha parameter compared to standard ACO
            print ("Exploitation")
            max_val = -1
            val = None

            for node in self.nodes_to_visit.values():
                if graph.tau(curr_node, node) == 0:
                    raise Exception("tau = 0")
                
                if graph.get_cluster(node) not in self.clusters_to_visit.values():
                    continue

                val = graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta)
                if val > max_val:
                    max_val = val
                    max_node = node
        else:
            # change here to stocastic probability rule
            print ("Exploration")
            probabilities = []
            nodes = []
            sum = 0.0

            for node in self.nodes_to_visit.values():
                if graph.get_cluster(node) not in self.clusters_to_visit:
                    continue
                tau = graph.tau(curr_node, node)
                etha = graph.etha(curr_node, node)
                if tau == 0:
                    raise Exception("tau = 0")
                
                probability = tau * math.pow(etha, self.Beta)
                probabilities.append(probability)
                nodes.append(node)
                sum += probability

            if sum == 0:
                raise Exception("sum = 0")

            # normalize prob
            probabilities = [p / sum for p in probabilities]

            # Select the next node based on the probabilities
            max_node = random.choices(nodes, weights=probabilities, k=1)[0]
            '''
            sum = 0
            node = -1

            for node in self.nodes_to_visit.values():
                if graph.tau(curr_node, node) == 0:
                    raise Exception("tau = 0")
                sum += graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta)
            if sum == 0:
                raise Exception("sum = 0")

            avg = sum / len(self.nodes_to_visit)

            print ("avg = %s" % (avg,))

            for node in self.nodes_to_visit.values():
                if graph.get_cluster(node) not in self.clusters_to_visit.values():
                    continue
                p = graph.tau(curr_node, node) * math.pow(graph.etha(curr_node, node), self.Beta) 
                if p > avg:
                    print ("p = %s" % (p,))
                    max_node = node

            if max_node == -1: # bug??
                max_node = node
        '''
        if max_node < 0:
            raise Exception("max_node < 0")

        if max_node in self.nodes_to_visit.values():
            del self.nodes_to_visit[max_node]

        if graph.get_cluster(max_node) in self.clusters_to_visit.values():
            del self.clusters_to_visit[graph.get_cluster(max_node)]
        
        return max_node

    # phermone update rule for indiv ants
    def local_updating_rule(self, curr_node, next_node):
        val = (1 - self.Rho) * self.colony.graph.tau(curr_node, next_node) + (self.Rho * self.colony.graph.tau0)
        self.colony.graph.update_tau(curr_node, next_node, val)

