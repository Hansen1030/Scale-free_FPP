#include "graph.h"
#include <iostream>

int main(int argc, char* argv[]) {
    vector<vector<double>> test;
    for (auto i=0;i<6;i++){
        for(auto j=0;j<6;j++){
            if(i==j){
                test[i][j]=0;
                }
                else{
                    test[i][j]=(i+j)/2;
                    }
            if(i==4&&j!=4){
                test[i][j]=0;
            }
            if(i==2 &&j!=2){
                test[i][j]=1;
            }
        }
    }
     vector<double> answer=find_shortest_path(1,5);
    for(auto tem:answer){
        std::cout<<tem;
    }
    return 0;
}