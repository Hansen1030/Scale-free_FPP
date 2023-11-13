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
        Graph(int total_nodes,int start_node, int target_node, double alpha, double gamma, int distribution_type, double random_index_1, double random_index_2);

        /*
        * find the shortest path and output to a file
        * start_node: start position
        * target_node: end position
        * 
        * output: 0 if run succesful, 1 otherwise
        */
       // Control function
        void generate_output(int algorithm, bool path);

       // support function
        void test();//test function
        vector<double> line_scan(int parts, int nodes_num);

        // edge matrix setting: tranformed the node, and consited the new position of starting point and ending point
        void node_transform_defult();
        void node_transform_equation(int k);
        void node_transform_Omega(int n, int k);
        void generate_new_matrix(int parts);

       // shorest path
        vector<double> find_shortest_path(int start, int target, bool road =false, vector<vector<double>>* edge_matrix=nullptr);
        double greedy_alg_poly(bool road =false);
        double bidirectional_dijkstra();
        // write in file function
        void write_in_file(double answer);
        void write_in_file_path(vector<double> path);
        // overwrite
        void write_in_file_path(vector<double> path, vector<double> cost);

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
        vector<vector<double>> edge_matrix_1;
        vector<vector<double>> edge_matrix_2;

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
};