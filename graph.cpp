#include "graph.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <utility> 
#include <limits>  
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
using namespace std;

Graph::Graph(int total_nodes_, int start, int target,double alpha_, double gamma_, int distribution_type_, double random_index_1_, double random_index_2_) {
    total_nodes = total_nodes_;
    start_node=start;
    target_node=target;
    alpha = alpha_;
    gamma = gamma_;
    distribution_type = distribution_type_;
    random_index_1 = random_index_1_;
    random_index_2 = random_index_2_;
    vector<double> output = line_scan(3, 3);
    for (int i = 0; i < 3; i++) {
        cout << output[i] << endl;
    }
     
}
void Graph::generate_output(int algorithm, bool path) {
    vector<double> answer;
    switch(algorithm){
        case 0:{// default
            node_transform_defult();
            answer=find_shortest_path(start_node,target_node,path,&edge_matrix);
        }
            break;
        case 1:{// height weight
            node_transform_equation(8);
            answer=find_shortest_path(start_node,target_node,path,&edge_matrix);
            break;
        }
        case 2:{// two direction
            int parts=3;
            int num_nodes=3;
            cout<<"generate_new_matrix?"<<endl;
            generate_new_matrix(parts);
            cout<<"line scan?"<<endl;
            vector<double> index= line_scan(parts,num_nodes);
            cout<<"seg?"<<endl;
            int gap = (parts-1)*(node_array.size()/parts);
            answer = find_shortest_path(index[0], int(index[0])+(node_array.size()/parts), path, &edge_matrix_1);
            vector<double> answer_rest = find_shortest_path(int(index[1])-gap, index[1], path, &edge_matrix_2);
            if(!path){
                answer[0]+=answer_rest[0]+index[2];
                break;
            }
            for (auto i : answer_rest) {
                i+=gap;
                answer.push_back(i);
            }
            if(!path) answer.push_back(index[2]);

            break;
        }
        default:
            node_transform_defult();
            break;
    }
    if(path){
        write_in_file_path(answer);
    }else{
        write_in_file(answer[0]);
    }
}


void Graph::node_transform_defult(){

    node_array.resize(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();
    }
    edge_matrix.resize(total_nodes);
    for (int i = 0; i < edge_matrix.size(); i++) {
        edge_matrix[i].resize(total_nodes);
    }
    double Omega=0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> d(1.0);
    for (int i = 0; i < total_nodes; i++) {
        for (int j = i + 1; j < total_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
            Omega = d(gen);
            edge_matrix[i][j] = (double) std::pow((double) abs(i-j),alpha) * Omega / (node_array[i] * node_array[j]);
            // cout<<edge_matrix[i][j]<<" ";
            edge_matrix[j][i] = edge_matrix[i][j];
        }
        // cout<<endl;
    }
}

void Graph::node_transform_equation(int k){
    vector<int> original_index;
    vector<double> temp;
    temp.resize(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        temp[i] = weight_generator();
    }
    double threshold = pow(k,1/gamma);
    // cout<<"node: "<<total_nodes<<" alpha:"<<alpha<<" gamma: "<<gamma<<endl;
    // cout<<"threshold:"<<threshold<<endl;
    for(auto i=0;i<temp.size();i++){
        if(i==start_node){
            original_index.push_back(i);
            node_array.push_back(temp[i]);
            start_node = node_array.size()-1;
        }
        if(i==target_node){
            original_index.push_back(i);
            node_array.push_back(temp[i]);
            target_node = node_array.size()-1;
        }
        if(temp[i]>threshold &&i!=start_node && i!=target_node){
            original_index.push_back(i);
            node_array.push_back(temp[i]);
        }
    }
    // cout<<"weight"<<endl;
    // for(auto i:node_array){
    //     cout<<i<<" ";
    // }
    // cout<<count<<" above threshold"<<endl;
    // cout<<"------------------------------------------------------"<<endl;
    // cout<<"new original index:";
    // for(auto i:original_index){
    //     cout<<i<<" ";
    // }
    // cout<<endl;
    cout<<"new total node: "<<node_array.size()<<endl;
    edge_matrix.resize(node_array.size());
    for (int i = 0; i < edge_matrix.size(); i++) {
        edge_matrix[i].resize(node_array.size());
    }
    for (int i = 0; i < node_array.size(); i++) {
        vector<double> temp;
        for (int j = i + 1; j < node_array.size(); j++) {
            // TODO: assign each edge a weight by the given node weight
            double Omega = random_num_gen(2, 1, 1);
            edge_matrix[i][j] = (double) std::pow((double) abs(original_index[i]-original_index[j]),alpha) * Omega / (node_array[i] * node_array[j]);
            edge_matrix[j][i] = edge_matrix[i][j];
        }
    }
}



void Graph::test(){
    double Omega=0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> d(1.0);
    for(auto i=0;i<total_nodes;i++){
        for(auto j=0;j<total_nodes;j++){
            Omega = d(gen);
        }
        cout<<Omega<<endl;
    }
    cout<<"test Omega:"<<Omega<<endl;
}

bool unique(vector<int> original_index, int index){
    for(auto i:original_index){
        if(i==index) return false;
    }
    return true;
}

void Graph::node_transform_Omega(int n, int k){
    double threshold = pow(k,1/gamma);
    double threshold_Omega = 4/(total_nodes*total_nodes);
    vector<int> original_index;
    node_array.resize(total_nodes); 
    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();
        // store the node index which the weight overed the threshold
        if(node_array[i]>threshold){
            original_index.push_back(i);
        }
    }
    if(unique(original_index,start_node)){
        original_index.push_back(start_node);
    }
    if(unique(original_index,target_node)){
        original_index.push_back(target_node);
    }

    cout<<"node after shrink: "<<original_index.size()<<endl;
    vector<vector<double>> temp_matrix;
    temp_matrix.resize(total_nodes);
    for (int i = 0; i < temp_matrix.size(); i++) {
        temp_matrix[i].resize(total_nodes);
    }
    // Priority queue to keep track of n smallest Omega values
    priority_queue<pair<double, pair<int, int>>> pq;

    for (int i = 0; i < total_nodes; i++) {
        vector<double> temp;
        for (int j = i + 1; j < total_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
            double Omega = random_num_gen(2, 1, 1);
            temp_matrix[i][j] = (double) std::pow((double) abs(i-j),alpha) * Omega / (node_array[i] * node_array[j]);
            temp_matrix[j][i] = temp_matrix[i][j];
            //add the smallest n Omega values to the priority queue
            if(Omega>threshold_Omega){
                continue;
            }
            pq.push({Omega, {i, j}});
            if (pq.size() > n) {
                pq.pop(); // Remove the largest
            }
        }
    }
    cout<<"Omega size: "<<pq.size()<<endl;
    // store the index which the Omega value is the smallest n
    while(pq.size()>0){
        pair<double, pair<int, int>> temp = pq.top();
        pq.pop();
        if(unique(original_index,temp.second.first)){
            original_index.push_back(temp.second.first);
        }
        if(unique(original_index,temp.second.second)){
            original_index.push_back(temp.second.second);
        }
    }
    cout<<"total node after added omega: "<<original_index.size()<<endl;
    std::sort(original_index.begin(), original_index.end());

    start_node=find(original_index.begin(), original_index.end(), start_node) - original_index.begin();
    target_node=find(original_index.begin(), original_index.end(), target_node) - original_index.begin();

    edge_matrix.resize(original_index.size());
    for (int i = 0; i < edge_matrix.size(); i++) {
        edge_matrix[i].resize(original_index.size());
    }
    for (int i = 0; i < original_index.size(); i++) {
        vector<double> temp;
        for (int j = i + 1; j < original_index.size(); j++) {
            edge_matrix[i][j]= temp_matrix[original_index[i]][original_index[j]];
            edge_matrix[j][i] = edge_matrix[i][j];
        }
    }
}


double Graph::greedy_alg_poly(bool road) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> d(1.0);
    int b = 3;
    node_array.resize(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();
    }
    int power = 0;
    int total_time = 0;
    int last = start_node;
    while(pow(b,power) < (total_nodes / b)){
        // std::cout << pow(b, power) << " " << (total_nodes / b) << endl;
        power++;
        int pos = 0;
        double largest = 0;
        for (int i = start_node; i < start_node + pow(b,power); i++) {
            if (node_array[i] > largest) {
                largest = node_array[i];
                pos = i;
            }
        }
        if (pos == last) continue;
        double Omega = d(gen);
        total_time += (double) std::pow((double) abs(pos - last),alpha) * Omega / (node_array[pos] * node_array[last]);
        last = pos;
        std::cout << "power: " << power << " " << "position: " << pos << " Time now: " << total_time << std::endl;
    }
    while(power--) {
        int pos = 0;
        double largest = 0;
        for (int i = (target_node - pow(b,power)); i <= target_node; i++) {
            if (node_array[i] > largest) {
                largest = node_array[i];
                pos = i;
            }
        }
        if (pos == last) continue;
        double Omega = d(gen);
        total_time += (double) std::pow((double) abs(pos - last),alpha) * Omega / (node_array[pos] * node_array[last]);
        last = pos;
        std::cout << "power: " << power << " " << "position: " << pos << " Time now: " << total_time << std::endl;
    }
    double Omega = d(gen);
    total_time += (double) std::pow((double) abs(last - target_node),alpha) * Omega / (node_array[target_node] * node_array[last]);
    return total_time;
}



vector<double> Graph::find_shortest_path(int start, int target, bool road, vector<vector<double>>* edge_matrix) {
    int n = edge_matrix->size();
    vector<double> dist(n, numeric_limits<double>::max());
    vector<int> prev(n, -1);
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    dist[start] = 0;
    pq.push({0, target});
    double current_dist=0;
    while (!pq.empty()) {
        current_dist = pq.top().first;
        int current_node = pq.top().second;
        pq.pop();
        if (current_node == target) {
            break;
        }
        if (current_dist > dist[current_node]) continue;

        for (int i = 0; i < n; ++i) {
            double weight = (*edge_matrix)[current_node][i];
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
    vector<double> path;

    for (int at = target; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());
    if(road){
        return path;   
    }
    return vector<double>{current_dist};
}

double Graph::bidirectional_dijkstra() {
    int n = edge_matrix.size();
    vector<double> dist_f(n, numeric_limits<double>::max()); // Forward distances
    vector<double> dist_b(n, numeric_limits<double>::max()); // Backward distances

    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq_f; // Forward queue
    priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq_b; // Backward queue

    dist_f[start_node] = 0;
    pq_f.push({0, start_node});

    dist_b[target_node] = 0;
    pq_b.push({0, target_node});

    double best_path = numeric_limits<double>::max();
    int meeting_node = -1;

    unordered_map<int, double> visited_f, visited_b;

    while (!pq_f.empty() && !pq_b.empty()) {
        // Forward search
        int current_f = pq_f.top().second;
        double current_dist_f = pq_f.top().first;
        pq_f.pop();

        if (visited_f.count(current_f) == 0) { // Skip already visited nodes
            visited_f[current_f] = current_dist_f;

            if (visited_b.count(current_f) > 0) {
                double potential_best_path = current_dist_f + visited_b[current_f];
                if (potential_best_path < best_path) {
                    best_path = potential_best_path;
                    meeting_node = current_f;
                }
            }

            for (int i = 0; i < n; ++i) {
                double weight = edge_matrix[current_f][i];
                if (weight >= 0) {
                    double new_dist = current_dist_f + weight;
                    if (new_dist < dist_f[i]) {
                        pq_f.push({new_dist, i});
                        dist_f[i] = new_dist;
                    }
                }
            }
        }

        // Backward search
        int current_b = pq_b.top().second;
        double current_dist_b = pq_b.top().first;
        pq_b.pop();

        if (visited_b.count(current_b) == 0) { // Skip already visited nodes
            visited_b[current_b] = current_dist_b;

            if (visited_f.count(current_b) > 0) {
                double potential_best_path = current_dist_b + visited_f[current_b];
                if (potential_best_path < best_path) {
                    best_path = potential_best_path;
                    meeting_node = current_b;
                }
            }

            for (int i = 0; i < n; ++i) {
                double weight = edge_matrix[i][current_b];
                if (weight >= 0) {
                    double new_dist = current_dist_b + weight;
                    if (new_dist < dist_b[i]) {
                        pq_b.push({new_dist, i});
                        dist_b[i] = new_dist;
                    }
                }
            }
        }
    }

    if (meeting_node == -1) {
        cout << "No path exists." << endl;
        return numeric_limits<double>::max();
    }

    return best_path;
}


void Graph::write_in_file(double answer){
    std::ofstream outfile;
    string filename = "./data_new/" + std::to_string(alpha) + "_" + std::to_string(gamma) + "_" + std::to_string(total_nodes);
    outfile.open(filename, std::ios::app);  // Open in append mode

    if (outfile.is_open()) {
        outfile << answer << "\n";  // Write the double value and a newline
        outfile.close();  // Close the file
    } else {
        std::cerr << "Unable to open the file: " << filename << std::endl;
    }
}
//-------------------------------------------------- spliting rule
void Graph::write_in_file_path(vector<double> path){
    std::ofstream outfile;
    string filename = "./path/" + std::to_string(alpha) + "_" + std::to_string(gamma) + "_" + std::to_string(total_nodes);
    outfile.open(filename, std::ios::app);  // Open in append mode

    if (outfile.is_open()) {
        for(auto i=1;i<path.size();i++){
            auto weight_l=node_array[i-1];
            auto weight_r=node_array[i];
            outfile << path[i-1] <<" "<<path[i]<<" "<<weight_l<<" "<<weight_r<<"\n";
        }
        outfile.close();
    } else {
        std::cerr << "Unable to open the file: " << filename << std::endl;
    }
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

vector<double> Graph::line_scan(int parts, int nodes_num) {
    node_array.resize(total_nodes);
    for (int i = 0; i<total_nodes;i++) {
        node_array[i] = weight_generator();
    }
    vector<double> output;
    output.resize(3);
    vector<pair<int, double>> begin;
    vector<pair<int, double>> end;
    begin.resize(nodes_num);
    end.resize(nodes_num);
    for (int i = 0; i < (node_array.size() / parts); i++) {
        pair<int, double> prev;
        for (int j = 0; j < nodes_num; j++) {
            if (begin[j].second <= node_array[i]) {
                prev = begin[j];
                begin[j].first = i;
                begin[j].second = node_array[i];
                for (int k = j + 1; k < nodes_num; k++) {
                    pair<int, double> mem = begin[k];
                    begin[k] = prev;
                    prev = mem;
                }
                break;
            }
        }
    }
    for (int i = (node_array.size() * (parts - 1) / parts); i < node_array.size(); i++) {
        pair<int, double> prev;
        for (int j = 0; j < nodes_num; j++) {
            if (end[j].second <= node_array[i]) {
                prev = end[j];
                end[j].first = i;
                end[j].second = node_array[i];
                for (int k = j + 1; k < nodes_num; k++) {
                    pair<int, double> mem = end[k];
                    end[k] = prev;
                    prev = mem;
                }
                break;
            }
        }
    }
    // for (int i = 0; i < parts; i++) {
    //     cout << begin[i].first << endl;
    // }
    // for (int i = 0; i < parts; i++) {
    //     cout << end[i].first << endl;
    // }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> d(1.0);
    output[2] = 9999.9;
    for (int i = 0; i < begin.size(); i++) {
        for (int j = 0; j < end.size(); j++) {
            double omega = d(gen);
            double length = (double) std::pow((double) abs(begin[i].first-end[j].first),alpha) * omega / (begin[i].second * end[j].second);
            if (length < output[2]) {
                output[0] = (double) begin[i].first;
                output[1] = (double) end[j].first;
                output[2] = length;
                // cout << output[2] << endl;
            }
        }
    }
    // for (int i = 0; i < 3; i++) {
    //     cout << output[0] << endl;
    // }
    return output;
}

void Graph::generate_new_matrix(int parts) {
    node_array.resize(total_nodes);
    for (int i = 0; i < total_nodes; i++) {
        node_array[i] = weight_generator();
    }
    int n_nodes = total_nodes / parts;
        edge_matrix_1.resize(n_nodes);
    for (int i = 0; i < edge_matrix_1.size(); i++) {
        edge_matrix_1[i].resize(n_nodes);
    }
    double Omega = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::exponential_distribution<> d(1.0);
    for (int i = 0; i < n_nodes; i++) {
        for (int j = i + 1; j < n_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
            Omega = d(gen);
            edge_matrix_1[i][j] = (double)std::pow((double)abs(i - j), alpha) * Omega / (node_array[i] * node_array[j]);
            // cout<<edge_matrix[i][j]<<" ";
            edge_matrix_1[j][i] = edge_matrix_1[i][j];
        }
    }

    edge_matrix_2.resize(n_nodes);
    for (int i = 0; i < edge_matrix_2.size(); i++) {
        edge_matrix_2[i].resize(n_nodes);
    }
    for (int i = 0; i < n_nodes; i++) {
        for (int j = i + 1; j < n_nodes; j++) {
            // TODO: assign each edge a weight by the given node weight
            Omega = d(gen);
            edge_matrix_2[i][j] = (double)std::pow((double)abs(i - j), alpha) * Omega / (node_array[i + (parts - 1) * n_nodes] * node_array[j + (parts - 1) * n_nodes]);
            // cout<<edge_matrix[i][j]<<" ";
            edge_matrix_2[j][i] = edge_matrix_2[i][j];
        }
    }
}
