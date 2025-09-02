#pragma once
#include <vector>
#include <array>
#include <unordered_set>
#include <unordered_map>
#include <iosfwd>
#include <iostream>
#include <map>

// ---- Constants --------------------------------------------------------------
const int DELIVERY_FAIL = -1;

// ---- Basic node types -------------------------------------------------------

struct Node {
    int index{};
    std::vector<int> value;      // length n
    std::vector<int> neighbors;  // exactly NODE_DEGREE neighbors
};

struct DRPTable {
    int Hmax{0}, Dmax{0};
    std::vector<long double> data;
    DRPTable() = default;
    DRPTable(int H, int D) { reset(H, D); }
    void reset(int H, int D) { Hmax = H; Dmax = D; data.assign((size_t)(H+1)*(D+1), 0.0L); }
    inline long double& at(int h, int d) { return data[(size_t)h*(Dmax+1) + d]; }
    inline long double  get(int h, int d) const { return data[(size_t)h*(Dmax+1) + d]; }
};

// Classification result for neighbors of u wrt destination t (simple distance)
struct NeighborClasses {
    std::vector<int> preferred;
    std::vector<int> spare;
};

// Classification results for the new strategic algorithms
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
    static const int NODE_DEGREE = 4;

    int n, k, V, E, diameter;
    std::vector<Node> nodes;
    std::array<int,NODE_DEGREE> degs{{4,4,4,4}};
    std::vector<std::array<int,NODE_DEGREE>> rev_dir;
    std::unordered_set<unsigned long long> faulty_edges;
    std::vector<std::map<int, DRPTable>> P;

    explicit metaTorus(int n_, int k_);

    // DRP API
    void computeDirectedRoutingProbabilities();
    void addFaultEdge(int u, int v);
    void clearFaults();
    int  alpha(int u, int v) const;
    long double getP(int u, int b_id, int h, int d) const;

    // Routing and Pathfinding API
    int bfs(int start_node_id, int end_node_id);
    int route_brute(int start_id, int target_id);
    int route_directed(int start_id, int target_id);
    int route_strategic_two(int start_id, int target_id);
    int route_strategic_three(int start_id, int target_id);
    std::vector<int> route_strategic_two_with_path(int start_id, int target_id);
    std::vector<int> route_strategic_three_with_path(int start_id, int target_id);
    bool hasFaultyLink(int u_id, int v_id) const; // Moved to public

    // Distance helpers
    int distance(const Node* a, const Node* b) const;
    int distance(int u, int v) const { return distance(&nodes[u], &nodes[v]); }

    int hammingDistance(const Node* a, const Node* b) const;
    int hammingDistance(int u, int v) const { return hammingDistance(&nodes[u], &nodes[v]); }

    // Neighbor Classification
    NeighborClasses classifyNeighbors(int u, int t) const;
    StrategicNeighborClasses2 classify_strategic_two(int c_id, int t_id);
    StrategicNeighborClasses3 classify_strategic_three(int c_id, int t_id);

    // Faults + print helpers
    void setRandomFaults(double ratio, unsigned long long seed = 0);
    std::vector<int> findAllReachableNodes(int start_node_id) const;
    void printAllDRP(std::ostream& os) const;
    void printOld() const;

    // Coord/id helpers
    int getId(const std::vector<int>& coord) const;
    std::vector<int> idToCoord(int id) const;

private:
    // Private recursive worker functions
    int brute(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth);
    int directed(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth);
    int strategic_two(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth);
    int strategic_three(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth);
    std::vector<int> strategic_two_with_path_worker(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited);
    std::vector<int> strategic_three_with_path_worker(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited);

    // Helpers for strategic classification
    int ringDistance(int a, int b, int k) const;
    int acrossRing(int h_start, int h_target, const std::vector<int>& need_to_visit) const;

    static unsigned long long edgeKey(int u, int v);
    void computeCoordinatesFromIds();
    void createNeighbors4();
    void computeReverseDirs();
    void validateNeighbors() const;
};

