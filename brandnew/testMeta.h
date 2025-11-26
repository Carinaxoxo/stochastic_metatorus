#ifndef BRANDNEW_TESTMETA_H
#define BRANDNEW_TESTMETA_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

// A structure to hold information about a single node.
struct Node {
    int id;
    std::vector<int> neighbors;
    std::map<int, int> neighbor_id_to_index;
};

// Helper struct for classifying neighbors
struct NeighborClasses {
    std::vector<int> preferred;
    std::vector<int> spare;
};

class testMeta {
public:
    // A public constant to represent a failed delivery.
    static const int DELIVERY_FAIL = -1;

    // Constructor
    testMeta(int n_dims, int k_arity);

    // Probability calculation
    void calcDirectedRoutingProbabilities();
    void printProbabilitiesAsTable(const std::string& filename) const;

    // --- Routing Algorithms ---
    std::vector<int> findFaultFreePathBFS(int start_node, int end_node) const;
    // Public function to initiate the brute force adaptive routing
    int route_brute(int start_id, int target_id) const;


    // --- Public Helper & Accessor Functions ---
    long long getV() const { return V; }
    void setRandomFaults(double fault_ratio, unsigned int seed);
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
    void createMetaTorusNeighborsDOR(); // Use the Meta Torus specific logic
    std::vector<int> getCoords(int id) const;
    int getId(const std::vector<int>& coords) const;
    void setDirectedProbability(int u_id, int b_id, int h, int d, long double prob);
    long double getDirectedProbability(int u_id, int b_id, int h, int d) const;
    bool hasFaultyLink(int u_id, int v_id) const;

    // --- Helpers for Brute Routing Algorithm ---
    int distance(int node_a_id, int node_b_id) const;
    bool isPreferred(int neighbor_id, int current_id, int target_id) const;

    // --- NEW: Helpers for the final Brute Routing Algorithm ---
    NeighborClasses classifyNeighbors(int current_id, int target_id) const;
    int brute(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) const;
};

#endif //BRANDNEW_TESTMETA_H

