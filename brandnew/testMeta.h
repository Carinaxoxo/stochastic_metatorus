#ifndef BRANDNEW_TESTMETA_H
#define BRANDNEW_TESTMETA_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

// A structure to hold information about a single node.
struct Node {
    int id;
    std::vector<int> neighbors; // Will contain the 4 Metatorus-like neighbors.
    std::map<int, int> neighbor_id_to_index; // Helper map for fast, O(1) access
};

class testMeta {
public:
    // Constructor
    testMeta(int n_dims, int k_arity);

    // Probability calculation and printing
    void calcDirectedRoutingProbabilities();
    void printProbabilitiesAsTable(const std::string& filename) const;

    // Pathfinding function
    std::vector<int> findFaultFreePathBFS(int start_node, int end_node) const;

    // --- Public Helper & Accessor Functions ---
    long long getV() const { return V; }
    void setRandomFaults(double fault_ratio, unsigned int seed);
    // MOVED to public for use in main.cpp
    std::string getCoordString(int id) const;

private:
    // --- Parameters ---
    int n; // Number of dimensions
    int k; // Arity
    long long V; // Total vertices (k^n)

    // --- Core Data Structures ---
    std::vector<Node> nodes;
    using ProbTableHD = std::vector<std::vector<long double>>;
    using ProbTableNeighbor = std::vector<ProbTableHD>;
    using ProbabilityStorage = std::vector<ProbTableNeighbor>;
    ProbabilityStorage directed_probabilities;
    std::unordered_map<long long, bool> faulty_links;

    // --- Private Helper Functions ---
    std::vector<int> getCoords(int id) const;
    int getId(const std::vector<int>& coords) const;
    void setDirectedProbability(int u_id, int b_id, int h, int d, long double prob);
    long double getDirectedProbability(int u_id, int b_id, int h, int d) const;
    bool hasFaultyLink(int u_id, int v_id) const;
};

#endif //BRANDNEW_TESTMETA_H