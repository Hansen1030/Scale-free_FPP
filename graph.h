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
        Graph(int total_nodes,int start_node, int target_node, double alpha, double gamma, int distribution_type, double random_index_1, double random_index_2, int t_);

        /*
        * find the shortest path and output to a file
        * start_node: start position
        * target_node: end position
        * 
        * output: 0 if run succesful, 1 otherwise
        */
        double find_shortest_path(bool road =false);
        double greedy_alg_poly(bool road =false);
        double bidirectional_dijkstra();
        // This function tranformed the node, and consited the new position of starting point and ending point
        void node_transform_defult();
        void node_transform_equation(int k);
        void node_transform_Omega(int n, int k);
        void write_in_file(double answer);
        void write_in_file_path(vector<int> path);
        void test();//test function

        ~Graph();
    private:
        /*
        * generate w_i
        * 
        * output: w_i's value as a double
        */
        double weight_generator();

        double random_num_gen(int random_type, double random_index_1, double random_index_2);

        // this matrix record the edge weight
        vector<vector<double>> edge_matrix;

        // this array record the w of each node
        vector<double> node_array;

        double findKthLargest(int t);
        
        int total_nodes;
        int start_node;
        int target_node;
        double alpha;
        double gamma;
        int distribution_type;
        double random_index_1;
        double random_index_2;
        int t;
};