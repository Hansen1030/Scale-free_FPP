#include "graph.h"
#include <iostream>

int main(int argc, char* argv[]) {
    // vector<vector<double>> test;
    // for (auto i=0;i<6;i++){
    //     for(auto j=0;j<6;j++){
    //         if(i==j){
    //             test[i][j]=0;
    //             }
    //             else{
    //                 test[i][j]=(i+j)/2;
    //                 }
    //         if(i==4&&j!=4){
    //             test[i][j]=0;
    //         }
    //         if(i==2 &&j!=2){
    //             test[i][j]=1;
    //         }
    //     }
    // }
    // Graph v;
    //  vector<double> answer=v.find_shortest_path(1,5);
    // for(auto tem:answer){
    //     std::cout<<tem;
    // }
    if (argc != 9) {
        std::cout << "please type the input as int start_node, int target_node, int total_nodes, int alpha, int gamma, int distribution_type, double random_index_1, double random_index_2" << std::endl;
        return 1;
    }
    Graph* g = new Graph(std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stoi(argv[8]));
    g->find_shortest_path(std::stoi(argv[1]), std::stoi(argv[2]));

    return 0;
}