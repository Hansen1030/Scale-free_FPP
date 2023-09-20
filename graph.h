#pragma once

#include <random>
#include <vector>
#include <cstdlib>

using std::vector;

class Graph {
    public:
        Graph() = default;

        /*
        * Construct the graph
        * total_nodes: total number of nodes graph contains
        * alpha: parameter alpha
        * gamma: parameter gamma
        * distribution_type: 
        *       random_type               random_index_1              random_index_2
        *  0     uniform                  lower bound                 upper bound
        *  1      normal                      mean                        std
        *  2     expnential                   mean                         \
        *  3      binomial                   times                    probability
        */ 
        Graph(int total_nodes, int alpha, int gamma, int distribution_type, double random_index_1, double random_index_2);

        /*
        * find the shortest path and output to a file
        * start_node: start position
        * target_node: end position
        * 
        * output: 0 if run succesful, 1 otherwise
        */
        int find_shortest_path(int start_node, int target_node);

    private:
        /*
        * generate w_i
        * 
        * output: w_i's value as a double
        */
        double weight_generator();

        // this matrix record the edge weight
        vector<vector<double>> edge_matrix;

        // this array record the w of each node
        vector<double> node_array;
        
        int total_nodes;
        int alpha;
        int gamma;
        int distribution_type;
        double random_index_1;
        double random_index_2;
};