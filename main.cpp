#include "graph.h"
#include <iostream>
#include <chrono>

// rewrite a main function which use GPU 
int main(int argc, char* argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (argc != 10) {
        std::cout << "please type the input as int start_node, int target_node, int total_nodes, double alpha, double gamma, int distribution_type, double random_index_1, double random_index_2, int sample_size" << std::endl;
        return 1;
    }
    for (int i = 0; i < std::stoi(argv[9]); i++) {
        
        Graph* g = new Graph(std::stoi(argv[3]),std::stoi(argv[1]),std::stoi(argv[2]), std::stod(argv[4]), std::stod(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8]), std::stoi(argv[9]));
        g->find_shortest_path(true);

        delete g;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Elapsed time: " << elapsed_time << " milliseconds" << std::endl;
    return 0;
}