#include "graph.h"
#include <iostream>
#include <algorithm>

Graph::Graph(int total_nodes, int alpha, int gamma, int distribution_type, double random_index_1, double random_index_2) {
    total_nodes = total_nodes;
    alpha = alpha;
    gamma = gamma;
    distribution_type = distribution_type;
    random_index_1 = random_index_1;
    random_index_2 = random_index_2;

    node_array.resize(total_nodes);
    edge_matrix.resize(total_nodes);
    for (int i = 0; i < edge_matrix.size(); i++) {
        edge_matrix[i].resize(total_nodes);
    }

    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();
    }

    for (int i = 0; i < total_nodes; i++) {
        for (int j = i + 1; j < total_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
        }
    }
}

int Graph::find_shortest_path(int start_node, int target_node) {
    // TODO: find the shortest path and output to a file
    return 0;
}

double Graph::weight_generator() {
    // TODO: generate node weight by the given parameter
}