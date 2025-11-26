//
// Created by Carina Z on 2025/11/19.
//

#ifndef RESTART_METATORUS_H
#define RESTART_METATORUS_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <cassert> // For assert
#include <random>  // For random number generation
#include <fstream> // For file I/O
#include <iomanip> // For stream manipulators (setw, etc.)
#include <sstream> // For string streams

struct Node {
    int id;
    std::vector<int> value;
    std::vector<int> neighbors;
    std::vector<int> link_status; // 1 for active, 0 for faulty
};

struct NeighborClasses {
    std::vector<int> preferred;
    std::vector<int> spare;
};

struct StrategicNeighborClasses2 {
    std::vector<int> shorter;
    std::vector<int> other;
};

struct StrategicNeighborClasses3 {
    std::vector<int> shorter;
    std::vector<int> same;
    std::vector<int> longer;
};


class metaTorus {
public:
    static const int DELIVERY_FAIL = -1;

    // --- Parameters ---
    int n; // Number of dimensions - 1
    int dimensions;
    int k; // Arity
    long long V; // Total vertices (k^n)
    int diameter; // Max distance in the network
    // Constructor
    metaTorus(int dims, int k_arity);
    std::vector<std::vector<std::vector<std::vector<long double>>>> directed_probabilities;


    int getId(const std::vector<int>& coord) const;
    std::string getCoordString(int id) const;

    void printNetwork() const;
    void setFaultyLinks(double fault_rate);
    std::vector<std::pair<int, int>> getFaultyLinks() const;

    int getDistance(int node1_id, int node2_id) const;
    int getHammingDistance(int node1_id, int node2_id) const;
    int getStrategicDistance(int node1_id, int node2_id) const;

    void calculateDirectedRoutingProbabilities();

    NeighborClasses classifyNeighbors(int current_node_id, int dest_node_id) const;

    StrategicNeighborClasses2 classifyNeighborsStrategic2(int current_node_id, int dest_node_id) const;
    StrategicNeighborClasses3 classifyNeighborsStrategic3(int current_node_id, int dest_node_id) const;

    std::vector<int> BFS(int start_node_id, int end_node_id) const;
    int findShortestPathBFS(int start_node_id, int end_node_id) const;

    //for testing
    int RecursiveStrategic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const;
    int findPathStrategic(int start_node_id, int end_node_id, bool verbose = false) const;

    int findPathBrute(int start_node_id, int end_node_id, bool verbose = false) const;
    int bruteRecursive(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const;

    int findPathBruteProbabilistic(int start_node_id, int end_node_id, bool verbose = false) const;
    int bruteRecursiveProbabilistic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const;

//    int RecursiveStrategic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) const;
//    int findPathStrategic(int start_node_id, int end_node_id) const;

    int findPathDirectedStrategic2(int start_id, int target_id, bool verbose) const;
    int directedRecursiveStrategic2(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const;

    int findPathDirectedStrategic3(int start_id, int target_id, bool verbose) const;
    int directedRecursiveStrategic3(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const;


private:
    std::mt19937 rand_generator; // Random number generator engine

    bool hasFaultyLink(int u, int v) const;

    std::vector<Node> nodes;
    void createMetaTorusNeighborsDOR();

    // --- Probability Helper Methods (from user) ---
    int getSpecificNeighbor(int nodeIndex, int neighborIndex) const;
    long double getDirectedProbability(int u_id, int b_id, int h, int d) const;
    void setDirectedProbability(int u_id, int b_id, int h, int d, long double prob);

};

#endif //RESTART_METATORUS_H
