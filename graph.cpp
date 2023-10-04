#include "graph.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <utility> 
#include <limits>  
#include <fstream>
using namespace std;

Graph::Graph(int total_nodes_, double alpha_, double gamma_, int distribution_type_, double random_index_1_, double random_index_2_) {
    total_nodes = total_nodes_;
    alpha = alpha_;
    gamma = gamma_;
    distribution_type = distribution_type_;
    random_index_1 = random_index_1_;
    random_index_2 = random_index_2_;

    node_array.resize(total_nodes);
    edge_matrix.resize(total_nodes);
    for (int i = 0; i < edge_matrix.size(); i++) {
        edge_matrix[i].resize(total_nodes);
    }

    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();

    }

    for (int i = 0; i < total_nodes; i++) {
        vector<double> temp;
        for (int j = i + 1; j < total_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
            double Omega = random_num_gen(2, 1, 1);
            edge_matrix[i][j] = (double) std::pow((double) abs(i-j),alpha) * Omega / (node_array[i] * node_array[j]);
            edge_matrix[j][i] = edge_matrix[i][j];
        }
    }
}

int Graph::find_shortest_path(int start_node, int target_node) {
    int n = edge_matrix.size();
    vector<double> dist(n, numeric_limits<double>::max());
    vector<int> prev(n, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;

    dist[start_node] = 0;
    pq.push({0, start_node});

    while (!pq.empty()) {
        double current_dist = pq.top().first;
        int current_node = pq.top().second;
        pq.pop();

        if (current_dist > dist[current_node]) continue;

        for (int i = 0; i < n; ++i) {
            double weight = edge_matrix[current_node][i];
            if (weight < 0) continue; // Skip negative weights or non-edges

            double new_dist = current_dist + weight;
            if (new_dist < dist[i]) {
                dist[i] = new_dist;
                prev[i] = current_node;
                pq.push({new_dist, i});
            }
        }
    }
    // Reconstruct the shortest path from start to end
    vector<int> path;
    for (int at = target_node; at != -1; at = prev[at]) {
        path.push_back(at);
    }

    reverse(path.begin(), path.end());
    int row=path[0];
    double answer = 0;
    for(int i=1;i<path.size();i++){
        answer += edge_matrix[row][path[i]];
        row = path[i];
    }
    std::ofstream outfile;
    string filename = "./data/" + std::to_string(alpha) + "_" + std::to_string(gamma) + "_" + std::to_string(n);
    outfile.open(filename, std::ios::app);  // Open in append mode

    if (outfile.is_open()) {
        outfile << answer << "\n";  // Write the double value and a newline
        outfile.close();  // Close the file
    } else {
        std::cerr << "Unable to open the file: " << filename << std::endl;
    }
    return 0;
}

double Graph::weight_generator() {
    // TODO: generate node weight by the given parameter
    double rand = random_num_gen(distribution_type, random_index_1, random_index_2);
    return pow(rand, (-1.0/gamma));
}

double Graph::random_num_gen(int random_type, double random_index_1, double random_index_2)
{
        std::random_device rd;
        std::mt19937 gen(rd());
        if (random_type == 0)
        {
                // uniform distribution
                std::uniform_real_distribution<> distribution(random_index_1, random_index_2);
                return distribution(gen);
        }
        else if (random_type == 1)
        {
                // normal distribution
                std::normal_distribution<double> distribution(random_index_1, random_index_2);
                return distribution(gen);
        }
        else if (random_type == 2)
        {
                // exponential distribution
                std::exponential_distribution<> distribution(random_index_1);
                return distribution(gen);
        }
        else if (random_type == 3)
        {
                // binomial distribution
                std::binomial_distribution<> distribution(random_index_1, random_index_2);
                return distribution(gen);
        }
        return -2;
}

Graph::~Graph() {
    node_array.clear();
    for (auto v : edge_matrix) {
        v.clear();
    }
    edge_matrix.clear();
}