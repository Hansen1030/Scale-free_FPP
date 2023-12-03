#include "graph.h"
#include <iostream>
#include <chrono>

int main() {

    auto start_time = std::chrono::high_resolution_clock::now();
    int total_nodes = 20000;
    int start_node = total_nodes*0.2;
    int target_node = total_nodes*0.8;
    double alpha = 2.2;
    double gamma = 0.8;
    int distribution_type = 0;
    double random_index_1 = 0;
    double random_index_2 = 1;
    int sample_size = 3;
    /*
    choose the algorithm:
    0: default
    1: heightest weight
    2: divide by part
    */ 
   int algorithm = 0;

   //genetrate path or cost:
   bool path = true;

    for (int i = 0; i < sample_size; i++) {
        Graph* g = new Graph(total_nodes, start_node, target_node, alpha, gamma, distribution_type, random_index_1, random_index_2);
        g->generate_output(algorithm,path);
        delete g;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()/1000.0;
    std::cout << "Elapsed time: " << elapsed_time << " second" << std::endl;
    return 0;
}