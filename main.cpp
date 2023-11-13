#include "graph.h"
#include <iostream>
#include <chrono>

// rewrite a main function which use GPU 
int main() {
    auto start_time = std::chrono::high_resolution_clock::now();
    int start_node = 2000;
    int target_node = 8000;
    int total_nodes = 10000;
    double alpha = 3;
    double gamma = 0.8;
    int distribution_type = 0;
    double random_index_1 = 0;
    double random_index_2 = 1;
    int sample_size = 1;
    /*
    choose the algorithm:
    0: default
    1: heightest weight
    2: divide by part
    */ 
   int algorithm = 2;

   //genetrate path or cost:
   bool path = true;

    for (int i = 0; i < sample_size; i++) {
        Graph* g = new Graph(total_nodes, start_node, target_node, alpha, gamma, distribution_type, random_index_1, random_index_2);
        g->generate_output(algorithm,path);
        delete g;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    std::cout << "Elapsed time: " << elapsed_time << " milliseconds" << std::endl;
    return 0;
}