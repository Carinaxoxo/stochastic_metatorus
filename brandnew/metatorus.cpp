#include "metatorus.h"

#include <cassert>
#include <limits>
#include <algorithm>
#include <random>
#include <iomanip>
#include <cmath>
#include <queue> // Required for std::queue

using namespace std;

// --------- internal helpers ---------

static long long ipowll(long long base, int exp) {
    long long res = 1;
    while (exp) {
        if (exp & 1) res *= base;
        base *= base;
        exp >>= 1;
    }
    return res;
}


unsigned long long metaTorus::edgeKey(int u, int v) {
    unsigned a = (unsigned)min(u, v), b = (unsigned)max(u, v);
    return ( (unsigned long long)a << 32 ) | b;
}

// --------- ctor ---------

metaTorus::metaTorus(int n_, int k_) : n(n_), k(k_), V(0), E(0), diameter(0) {
    assert(n > 2 && "n must be > 2");
    assert(k > 2 && "k must be > 2");
    assert(k == (n - 1) && "This construction assumes k == n-1");

    long long vll = ipowll(k, n);
    assert(vll <= (long long)numeric_limits<int>::max());
    V = (int)vll;
    E = 2 * V;
    diameter = (int)((k / 2) * n);

    nodes.resize(V);
    for (int i = 0; i < V; ++i) {
        nodes[i].index = i;
        nodes[i].value.assign(n, 0);
    }

    computeCoordinatesFromIds();
    createNeighbors4();
    computeReverseDirs();
    validateNeighbors();

    int Hmax = n;
    int Dmax = diameter;
    P.resize(V);
    for (int u = 0; u < V; ++u)
        for (int b = 0; b < NODE_DEGREE; ++b)
            P[u][b] = DRPTable(Hmax, Dmax);
}

// --------- public API ---------

void metaTorus::addFaultEdge(int u, int v) {
    faulty_edges.insert(edgeKey(u, v));
}

void metaTorus::clearFaults() {
    faulty_edges.clear();
}

int metaTorus::alpha(int u, int v) const {
    return faulty_edges.count(edgeKey(u, v)) ? 0 : 1;
}


long double metaTorus::getP(int u, int b_id, int h, int d) const {
    // Check if h and d are within the valid range to prevent crashes.
    if (h < 0 || h > n || d < 0 || d > diameter) {
        return -1.0L;
    }
    // Find the DRP table associated with the incoming hop 'b_id'.
    auto it = P[u].find(b_id);
    if (it != P[u].end()) {
        return it->second.get(h, d);
    }
    return -1.0L; // Should not happen in a correct implementation
}

int metaTorus::getId(const vector<int>& coord) const {
    long long id = 0;
    for (int j = 0; j < n; ++j) {
        assert(0 <= coord[j] && coord[j] < k);
        id = id * k + coord[j];
    }
    assert(id >= 0 && id < V);
    return (int)id;
}

vector<int> metaTorus::idToCoord(int id) const {
    vector<int> coord(n, 0);
    int m = id;
    for (int j = n - 1; j > 0; --j) {
        int po = (int)ipowll(k, j);
        coord[n - j - 1] = m / po;
        m %= po;
    }
    coord[n - 1] = m;
    return coord;
}

// --------- DRP computation ---------

void metaTorus::computeDirectedRoutingProbabilities() {
    const int Dmax = diameter;
    const int Hcap = n;
    const int denom = 2 * n - 1;

    P.assign(V, map<int, DRPTable>());

    // Initialize DRP tables for all nodes and all incoming directions
    for (int u = 0; u < V; ++u) {
        for (int b_neighbor_id : nodes[u].neighbors) {
            P[u][b_neighbor_id] = DRPTable(Hcap, Dmax);
        }
        // Also need a dummy entry for the source node's initial call
        P[u][u] = DRPTable(Hcap, Dmax);
    }

    // Base case: (h, d) = (1, 1)
    for (int u = 0; u < V; ++u) {
        // For each possible incoming direction 'b'...
        for (int b_neighbor_id : nodes[u].neighbors) {
            long double numerator = 0.0L;
            // Sum over all *other* outgoing directions 'a'.
            for (int a_neighbor_id : nodes[u].neighbors) {
                if (a_neighbor_id == b_neighbor_id) continue;
                numerator += (long double)alpha(u, a_neighbor_id);
            }
            P[u][b_neighbor_id].at(1, 1) = numerator / (long double)denom;
        }
    }

    // Inductive steps: d = 2..Dmax, h = 1..min(d, n)
    for (int d = 2; d <= Dmax; ++d) {
        int hmax = min(d, Hcap);
        for (int h = 1; h <= hmax; ++h) {
            for (int u = 0; u < V; ++u) {
                // First, calculate the total sum S_h,d(u)
                long double Shd = 0.0L;
                for (int a_dir = 0; a_dir < NODE_DEGREE; ++a_dir) {
                    int a_id = nodes[u].neighbors[a_dir];
                    long double term_a = 0.0L;
                    // The probability at neighbor 'a' depends on the incoming direction being 'u'
                    if (h == 1) {
                        term_a = getP(a_id, u, 1, d - 1);
                    } else if (h == d) {
                        term_a = getP(a_id, u, h - 1, d - 1);
                    } else {
                        term_a = ((long double)(h - 1) / (d - 1)) * getP(a_id, u, h - 1, d - 1)
                                 + ((long double)(d - h) / (d - 1)) * getP(a_id, u, h, d - 1);
                    }
                    Shd += (long double)alpha(u, a_id) * term_a;
                }

                // Now, calculate P_h,d^b(u) for each incoming direction 'b'
                for (int b_neighbor_id : nodes[u].neighbors) {
                    long double term_b = 0.0L;
                    if (h == 1) {
                        term_b = getP(b_neighbor_id, u, 1, d - 1);
                    } else if (h == d) {
                        term_b = getP(b_neighbor_id, u, h - 1, d - 1);
                    } else {
                        term_b = ((long double)(h - 1) / (d - 1)) * getP(b_neighbor_id, u, h - 1, d - 1)
                                 + ((long double)(d - h) / (d - 1)) * getP(b_neighbor_id, u, h, d - 1);
                    }

                    long double num = Shd - (long double)alpha(u, b_neighbor_id) * term_b;
                    P[u][b_neighbor_id].at(h, d) = num / (long double)denom;
                }
            }
        }
    }
}

/**
 * @brief Calculates the directed routing probabilities for all nodes in the torus.
 *
 * This function is corrected to use the project's variable naming convention:
 * - n: number of dimensions
 * - k: arity (number of nodes per dimension)
 *
 * It implements the advanced algorithm from the paper "A Stochastic
 * Edge-fault-tolerant Routing Algorithm in Torus Networks" by Zhang and Kaneko.
 */



void metaTorus::calcDirectedRoutingProbabilities() {
    // --- Correct Parameterization using YOUR variable names ---
    // n: number of dimensions
    // k: arity (nodes per dimension)

    // Degree of each node in an n-dimensional torus is 2*n.
    const int degree = 2 * n;
    // Diameter is dimension * floor(arity / 2).
    const int diameter = n * (k / 2);

    // --- Step 1: Base Case Calculation (d=1, h=1) ---
    // This must be completed for ALL nodes before moving to the next step,
    // as the main loop depends on these initial values from neighbors.
    for (int i = 0; i < V; ++i) {
        Node* u = &nodes[i];

        // For each possible incoming direction from a neighbor 'b'...
        for (int b_idx = 0; b_idx < degree; ++b_idx) {
            Node* b = &nodes[u->neighbors[b_idx]];
            double sum_alpha = 0.0;

            // ...sum the link statuses of all OTHER outgoing neighbors 'a'.
            for (int a_idx = 0; a_idx < degree; ++a_idx) {
                Node* a = &nodes[u->neighbors[a_idx]];
                if (a->id == b->id) {
                    continue; // Exclude the incoming path
                }
                sum_alpha += !hasFaultyLink(u, a) ? 1.0 : 0.0;
            }

            // Probability is normalized by the number of valid outgoing paths (degree - 1).
            double p_1_1 = sum_alpha / (degree - 1.0);
            setDirectedProbability(u, b, 1, 1, p_1_1);
        }
    }

    // --- Step 2: Main Loop (Dynamic Programming) ---
    // Iterate through increasing path lengths (d) and Hamming distances (h).
    for (int d = 2; d <= diameter; ++d) {
        // --- CRITICAL FIX: Hamming distance 'h' is limited by the number of dimensions 'n'. ---
        for (int h = 1; h <= std::min(d, n); ++h) {

            // This loop completes for all nodes 'u' at a given (h,d) before the
            // next iteration, simulating synchronized message passing.
            for (int i = 0; i < V; ++i) {
                Node* u = &nodes[i];

                // Part A: Pre-calculate the total sum 'S_hd' over ALL neighbors 'a'.
                double S_hd = 0.0;
                for (int a_idx = 0; a_idx < degree; ++a_idx) {
                    Node* a = &nodes[u->neighbors[a_idx]];
                    double alpha = !hasFaultyLink(u, a) ? 1.0 : 0.0;
                    if (alpha == 0.0) continue;

                    // Apply the correct recurrence relation from the paper
                    if (h == 1) {
                        S_hd += alpha * getDirectedProbability(a, u, h, d - 1);
                    } else if (h == d) {
                        S_hd += alpha * getDirectedProbability(a, u, h - 1, d - 1);
                    } else { // General case
                        double p1 = getDirectedProbability(a, u, h - 1, d - 1);
                        double p2 = getDirectedProbability(a, u, h, d - 1);
                        if (d > 1) {
                            S_hd += alpha * (((h - 1.0) * p1 + (d - h) * p2) / (d - 1.0));
                        }
                    }
                }

                // Part B: Calculate the specific directed probability P(h,d,b)(u).
                // For each incoming direction 'b', subtract its term from the total sum 'S_hd'.
                for (int b_idx = 0; b_idx < degree; ++b_idx) {
                    Node* b = &nodes[u->neighbors[b_idx]];
                    double alpha_b = !hasFaultyLink(u, b) ? 1.0 : 0.0;
                    double term_for_b = 0.0;

                    if (alpha_b > 0.0) { // Calculate the term 'b' would have contributed
                        if (h == 1) {
                            term_for_b = alpha_b * getDirectedProbability(b, u, h, d - 1);
                        } else if (h == d) {
                            term_for_b = alpha_b * getDirectedProbability(b, u, h - 1, d - 1);
                        } else { // General case
                            double p1 = getDirectedProbability(b, u, h - 1, d - 1);
                            double p2 = getDirectedProbability(b, u, h, d - 1);
                            if (d > 1) {
                                term_for_b = alpha_b * (((h - 1.0) * p1 + (d - h) * p2) / (d - 1.0));
                            }
                        }
                    }

                    // The directed probability is the sum from all *other* paths, normalized.
                    double p_hd = (S_hd - term_for_b) / (degree - 1.0);
                    setDirectedProbability(u, b, h, d, p_hd);
                }
            }
        }
    }
}

// --------- Routing and Pathfinding Implementation ---------

bool metaTorus::hasFaultyLink(int u_id, int v_id) const {
    // A link is faulty if its key exists in the faulty_edges set.
    return faulty_edges.count(edgeKey(u_id, v_id)) > 0;
}

int metaTorus::bfs(int start_node_id, int end_node_id) {
    if (start_node_id == end_node_id) {
        return 0;
    }

    // Queue stores pairs of {node_id, distance_from_start}
    queue<pair<int, int>> q;
    // Visited set to prevent cycles and redundant explorations
    unordered_set<int> visited;

    q.push({start_node_id, 0});
    visited.insert(start_node_id);

    while (!q.empty()) {
        pair<int, int> current = q.front();
        q.pop();
        int current_id = current.first;
        int current_dist = current.second;

        if (current_id == end_node_id) {
            return current_dist; // Found the shortest path
        }

        // Explore neighbors
        for (int neighbor_id : nodes[current_id].neighbors) {
            // Check if neighbor has been visited AND if the link is fault-free
            if (visited.find(neighbor_id) == visited.end() && !hasFaultyLink(current_id, neighbor_id)) {
                visited.insert(neighbor_id);
                q.push({neighbor_id, current_dist + 1});
            }
        }
    }

    return -1; // No path found, return -1 to indicate failure
}


// --------- private topology builders ---------

void metaTorus::computeCoordinatesFromIds() {
    for (int i = 0; i < V; ++i) {
        int m = i;
        for (int j = n - 1; j > 0; --j) {
            int po = (int)ipowll(k, j);
            nodes[i].value[n - j - 1] = m / po;
            m %= po;
        }
        nodes[i].value[n - 1] = m;
    }
}

void metaTorus::createNeighbors4() {
    for (int i = 0; i < V; ++i) {
        const vector<int> current = nodes[i].value;
        nodes[i].neighbors.assign(NODE_DEGREE, -1); // Initialize with -1

        // A temporary structure to hold neighbor and the dimension it affects
        vector<pair<int, int>> temp_neighbors;

        int header = current[0];
        int x = header + 1;
        assert(1 <= x && x <= (n - 1));

        // Neighbor on header ring (dimension 0) +1
        vector<int> temp3_coords = current;
        temp3_coords[0] = (header + 1) % (n - 1);
        temp_neighbors.push_back({0, getId(temp3_coords)});

        // Neighbor on header ring (dimension 0) -1
        vector<int> temp4_coords = current;
        temp4_coords[0] = (header - 1 + (n - 1)) % (n - 1);
        temp_neighbors.push_back({0, getId(temp4_coords)});

        // Neighbor in dimension x, +1
        vector<int> temp1_coords = current;
        temp1_coords[x] = (temp1_coords[x] + 1) % k;
        temp_neighbors.push_back({x, getId(temp1_coords)});

        // Neighbor in dimension x, -1
        vector<int> temp2_coords = current;
        temp2_coords[x] = (temp2_coords[x] - 1 + k) % k;
        temp_neighbors.push_back({x, getId(temp2_coords)});

        // Sort the neighbors based on the dimension they affect
        std::sort(temp_neighbors.begin(), temp_neighbors.end());

        // Assign the sorted neighbor IDs to the node
        for(int j=0; j < NODE_DEGREE; ++j) {
            nodes[i].neighbors[j] = temp_neighbors[j].second;
        }
    }
}

void metaTorus::computeReverseDirs() {
    rev_dir.assign(V, array<int,NODE_DEGREE>{-1, -1, -1, -1});
    for (int u = 0; u < V; ++u) {
        for (int dir = 0; dir < NODE_DEGREE; ++dir) {
            int v = nodes[u].neighbors[dir];
            int r = -1;
            for (int t = 0; t < NODE_DEGREE; ++t) {
                if (nodes[v].neighbors[t] == u) { r = t; break; }
            }
            assert(r != -1 && "Neighbor relationship must be symmetric");
            rev_dir[u][dir] = r;
        }
    }
}

void metaTorus::validateNeighbors() const {
    for (int i = 0; i < V; ++i) {
        const auto& nbrs = nodes[i].neighbors;
        assert((int)nbrs.size() == NODE_DEGREE);
        for (int t = 0; t < NODE_DEGREE; ++t) {
            assert(0 <= nbrs[t] && nbrs[t] < V);
            assert(nbrs[t] != i);
            for (int u = t + 1; u < NODE_DEGREE; ++u) assert(nbrs[t] != nbrs[u]);
        }
    }
}

int metaTorus::distance(const Node* a, const Node* b) const {
    int d = 0;
    for (int i = 0; i < n; ++i) {
        int diff = std::abs(a->value[i] - b->value[i]);
        d += std::min(diff, k - diff);
    }
    assert(0 <= d && d <= diameter);
    return d;
}

int metaTorus::hammingDistance(const Node* a, const Node* b) const {
    int d = 0;
    for (int i = 0; i < n; ++i)
        if (a->value[i] != b->value[i]) ++d;
    assert(0 <= d && d <= n);
    return d;
}

NeighborClasses metaTorus::classifyNeighbors(int u, int t) const {
    NeighborClasses out;
    const int dist_ut = distance(u, t);
    for (int v : nodes[u].neighbors) {
        int dist_vt = distance(v, t);
        if (dist_vt < dist_ut) out.preferred.push_back(v);
        else                    out.spare.push_back(v);
    }
    return out;
}

void metaTorus::setRandomFaults(double ratio, unsigned long long seed) {
    if (ratio < 0.0) ratio = 0.0;
    if (ratio > 1.0) ratio = 1.0;
    clearFaults();
    std::vector<std::pair<int,int>> edges;
    edges.reserve(static_cast<size_t>(E));
    for (int u = 0; u < V; ++u) {
        for (int dir = 0; dir < NODE_DEGREE; ++dir) {
            int v = nodes[u].neighbors[dir];
            if (u < v) edges.emplace_back(u, v);
        }
    }
    const size_t total = edges.size();
    const size_t m = static_cast<size_t>(std::llround(ratio * total));
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::shuffle(edges.begin(), edges.end(), rng);
    for (size_t i = 0; i < m; ++i) {
        addFaultEdge(edges[i].first, edges[i].second);
    }
}

std::vector<int> metaTorus::findAllReachableNodes(int start_node_id) const {
    std::vector<int> reachable_nodes;
    std::queue<int> q;
    std::unordered_set<int> visited;

    q.push(start_node_id);
    visited.insert(start_node_id);

    while (!q.empty()) {
        int current_id = q.front();
        q.pop();
        reachable_nodes.push_back(current_id);

        // Explore neighbors
        for (int neighbor_id : nodes[current_id].neighbors) {
            // Check if neighbor has been visited AND if the link is fault-free
            if (visited.find(neighbor_id) == visited.end() && !hasFaultyLink(current_id, neighbor_id)) {
                visited.insert(neighbor_id);
                q.push(neighbor_id);
            }
        }
    }
    return reachable_nodes;
}

void metaTorus::printAllDRP(std::ostream& os) const {
    os.setf(std::ios::fixed);
    os << std::setprecision(6);
    os << "u,b,h,d,P\n";
    for (int u = 0; u < V; ++u) {
        for (int b = 0; b < NODE_DEGREE; ++b) {
            for (int d = 1; d <= diameter; ++d) {
                int hmax = std::min(d, n);
                for (int h = 1; h <= hmax; ++h) {
                    os << u << ',' << b << ',' << h << ',' << d << ','
                       << static_cast<double>(getP(u, b, h, d)) << '\n';
                }
            }
        }
    }
}

void metaTorus::printOld() const {
    std::cout << "Torus: n=" << n << ", k=" << k << std::endl;
    std::cout << "Number of vertices (V): " << V << std::endl;
    std::cout << "Number of edges (E): " << E << std::endl;
    std::cout << "Diameter: " << diameter << std::endl;
    std::cout << "\nNode coordinates and neighbors:" << std::endl;
    for (int i = 0; i < V; ++i) {
        std::cout << "Node " << i << " (";
        for (int j = 0; j < n; ++j) {
            std::cout << nodes[i].value[j];
            if (j < n - 1) std::cout << ",";
        }
        std::cout << "): " << std::endl;
        for (int j = 0; j < NODE_DEGREE; ++j) {
            int neighbor_id = nodes[i].neighbors[j];
            int faultiness = 1 - alpha(i, neighbor_id);
            std::cout << "  Neighbor " << j << " (" << neighbor_id << " -> ";
            for (int k = 0; k < n; ++k) {
                std::cout << nodes[neighbor_id].value[k];
                if (k < n - 1) std::cout << ",";
            }
            std::cout << ") Faultiness: " << faultiness << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "\n--- Directed Routing Probabilities (DRP) ---" << std::endl;
    for (int i = 0; i < V; i++) {
        for (int h = 1; h <= n; h++) {
            for (int d = h; d <= diameter; d++) {
                for (int index = 0; index < NODE_DEGREE; index++) {
                    long double prob = getP(i, index, h, d);
                    if (prob > 0.0L) {
                        std::cout << "nodes[" << i << "].d_P[" << h << "][" << d << "][" << index << "] = "
                                  << std::fixed << std::setprecision(6) << static_cast<double>(prob) << std::endl;
                    }
                }
            }
        }
    }
    std::cout << std::endl;
}

/**
 * @brief The recursive, greedy routing algorithm as you designed it.
 * It commits to the first valid path it finds.
 */
int metaTorus::route_brute(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    // Initial call: no previous node (-1), start at depth 0.
    return brute(-1, start_id, target_id, visited, 0);
}

/**
 * @brief The corrected brute algorithm. It is greedy, stateful (for cycle
 * detection), and uses deterministic dimension-order routing.
 */
int metaTorus::brute(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) {
    if (current_id == target_id) return depth;

    visited[current_id] = true;

    NeighborClasses classes = classifyNeighbors(current_id, target_id);

    // --- Step 1: Check preferred neighbors in dimension order ---
    // The neighbors are already sorted by dimension, so we just iterate.
    for (int neighbor_id : classes.preferred) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            // Check for a cycle before making the move.
            if (visited.count(neighbor_id)) {
                return DELIVERY_FAIL; // Fail immediately if a loop is detected
            }
            // The first valid neighbor found is the deterministic choice.
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    // --- Step 2: If no valid preferred, check spare neighbors in dimension order ---
    for (int neighbor_id : classes.spare) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            // Check for a cycle before making the move.
            if (visited.count(neighbor_id)) {
                return DELIVERY_FAIL; // Fail immediately if a loop is detected
            }
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    // --- Step 3: If no valid moves at all ---
    return DELIVERY_FAIL;
}

int metaTorus::route_directed(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    // Initial call: no previous node (-1), start at depth 0.
    return directed(-1, start_id, target_id, visited, 0);
}

/**
 * @brief Recursive DRP-based routing algorithm.
 * Finds the neighbor with the highest probability, checks its validity, and recurses.
 * Backtracks if a high-probability path leads to a dead end.
 */
int metaTorus::directed(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) {
    if (current_id == target_id) return depth;

    visited[current_id] = true;

    NeighborClasses classes = classifyNeighbors(current_id, target_id);
    int dist_ct = distance(current_id, target_id);

    // Find BEST preferred neighbor
    int best_preferred_id = -1;
    long double max_prob_p = -1.0L;
    for (int neighbor_id : classes.preferred) {
        if (neighbor_id == prev_id) continue;
        int h = hammingDistance(neighbor_id, target_id);
        // Probability at 'neighbor_id' given it came from 'current_id'
        long double prob = getP(neighbor_id, current_id, h, dist_ct - 1);
        if (prob > max_prob_p) {
            max_prob_p = prob;
            best_preferred_id = neighbor_id;
        }
    }

    if (best_preferred_id != -1 && !hasFaultyLink(current_id, best_preferred_id)) {
        if (visited.count(best_preferred_id) != 0) return DELIVERY_FAIL;
        return directed(current_id, best_preferred_id, target_id, visited, depth + 1);
    }

    // Find BEST spare neighbor
    int best_spare_id = -1;
    long double max_prob_s = -1.0L;
    for (int neighbor_id : classes.spare) {
        if (neighbor_id == prev_id) continue;
        int h = hammingDistance(neighbor_id, target_id);
        long double prob = getP(neighbor_id, current_id, h, std::min(dist_ct + 1, diameter));
        if (prob > max_prob_s) {
            max_prob_s = prob;
            best_spare_id = neighbor_id;
        }
    }

    if (best_spare_id != -1 && !hasFaultyLink(current_id, best_spare_id)) {
        if (visited.count(best_spare_id) != 0) return DELIVERY_FAIL;
        return directed(current_id, best_spare_id, target_id, visited, depth + 1);
    }

    return DELIVERY_FAIL;
}

// --------- Strategic Classification Helpers ---------

int metaTorus::ringDistance(int a, int b, int k) const {
    int diff = std::abs(a - b);
    return std::min(diff, k - diff);
}

int metaTorus::acrossRing(int h_start, int h_target, const std::vector<int>& need_to_visit) const {
    int num = need_to_visit.size();
    if (num == 0) return 0;
    std::vector<int> temp;
    temp.push_back(h_start);

    for (int i = 1; i <= (n-1); i++) {
        int pos = h_start + i;
        if (need_to_visit[pos % (n - 1)] == 1 || pos % (n - 1) == h_target)
            temp.push_back(pos);
    }
    temp.push_back(h_start + n-1);

    int side1 = 0, side2 = 0, half1 = 0, half2 = 0;
    bool switch_side = false;

    for (size_t i = 1; i < temp.size(); i++) {
        int delta = temp[i] - temp[i - 1];
        if (!switch_side){
            side1 = std::max(side1, delta);
        }
        else{
            side2 = std::max(side2, delta);
        }
        if (temp[i] % (n-1)  == h_target) {
            switch_side = true;
            half1 = temp[i] - h_start;
            half2 = n-1 - half1;
        }
    }

    if (half1 == 0 && half2 == 0 && h_start != h_target) {
        return ringDistance(h_start, h_target, n-1);
    }
    if (half1 + half2 != n-1) return ringDistance(h_start, h_target, n-1);


    return std::min(half1 + half2 * 2 - side2 * 2, half2 + half1 * 2 - side1 * 2);
}


// --------- Strategic Classification Implementations (Corrected) ---------

StrategicNeighborClasses2 metaTorus::classify_strategic_two(int c_id, int t_id) {
    const Node& c = nodes[c_id];
    const Node& t = nodes[t_id];
    StrategicNeighborClasses2 classes;

    const auto& standing = c.value;
    const auto& destination = t.value;

    // Iterate through each neighbor and dynamically determine its type and cost
    for (int neighbor_id : c.neighbors) {
        const auto& neighbor_coords = nodes[neighbor_id].value;
        int changed_dim = -1;
        // Find which dimension changed to create this neighbor
        for (int i = 0; i < n; ++i) {
            if (standing[i] != neighbor_coords[i]) {
                changed_dim = i;
                break;
            }
        }

        if (changed_dim == 0) { // It's a RING move
            std::vector<int> need_to_visit(n - 1, 0);
            for(int i = 0; i < (n - 1); i++) {
                if(standing[i + 1] != destination[i + 1]) need_to_visit[i] = 1;
            }
            int ring_cost = acrossRing(standing[0], destination[0], need_to_visit);
            int neighbor_ring_cost = acrossRing(neighbor_coords[0], destination[0], need_to_visit);

            if (neighbor_ring_cost < ring_cost) {
                classes.shorter.push_back(neighbor_id);
            } else {
                classes.other.push_back(neighbor_id);
            }
        } else if (changed_dim > 0) { // It's a DIMENSION move
            int d_dim = ringDistance(standing[changed_dim], destination[changed_dim], k);
            int d_neighbor_dim = ringDistance(neighbor_coords[changed_dim], destination[changed_dim], k);

            if (d_neighbor_dim < d_dim) {
                classes.shorter.push_back(neighbor_id);
            } else {
                classes.other.push_back(neighbor_id);
            }
        }
    }
    return classes;
}

StrategicNeighborClasses3 metaTorus::classify_strategic_three(int c_id, int t_id) {
    const Node& c = nodes[c_id];
    const Node& t = nodes[t_id];
    StrategicNeighborClasses3 classes;

    const auto& standing = c.value;
    const auto& destination = t.value;

    for (int neighbor_id : c.neighbors) {
        const auto& neighbor_coords = nodes[neighbor_id].value;
        int changed_dim = -1;
        for (int i = 0; i < n; ++i) {
            if (standing[i] != neighbor_coords[i]) {
                changed_dim = i;
                break;
            }
        }

        if (changed_dim == 0) { // RING move
            std::vector<int> need_to_visit(n - 1, 0);
            for(int i = 0; i < (n - 1); i++) {
                if(standing[i + 1] != destination[i + 1]) need_to_visit[i] = 1;
            }
            int ring_cost = acrossRing(standing[0], destination[0], need_to_visit);
            int neighbor_ring_cost = acrossRing(neighbor_coords[0], destination[0], need_to_visit);

            if (neighbor_ring_cost < ring_cost) classes.shorter.push_back(neighbor_id);
            else if (neighbor_ring_cost == ring_cost) classes.same.push_back(neighbor_id);
            else classes.longer.push_back(neighbor_id);

        } else if (changed_dim > 0) { // DIMENSION move
            int d_dim = ringDistance(standing[changed_dim], destination[changed_dim], k);
            int d_neighbor_dim = ringDistance(neighbor_coords[changed_dim], destination[changed_dim], k);

            if (d_neighbor_dim < d_dim) classes.shorter.push_back(neighbor_id);
            else if (d_neighbor_dim == d_dim) classes.same.push_back(neighbor_id);
            else classes.longer.push_back(neighbor_id);
        }
    }
    return classes;
}

// --------- Strategic Routing Implementations ---------

int metaTorus::route_strategic_two(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    return strategic_two(-1, start_id, target_id, visited, 0);
}

int metaTorus::strategic_two(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) {
    if (current_id == target_id) return depth;
    visited[current_id] = true;

    StrategicNeighborClasses2 classes = classify_strategic_two(current_id, target_id);

    // --- Step 1: Try best neighbor from 'shorter' class ---
    int best_shorter_id = -1;
    long double max_prob_shorter = -1.0L;
    for (int neighbor_id : classes.shorter) {
        if (neighbor_id == prev_id) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        // [FIXED] Use 'current_id' as the incoming node for the DRP lookup.
        long double prob = getP(neighbor_id, current_id, h, d);

        if (prob > max_prob_shorter) {
            max_prob_shorter = prob;
            best_shorter_id = neighbor_id;
        }
    }
    if (best_shorter_id != -1 && !hasFaultyLink(current_id, best_shorter_id)) {
        if (visited.count(best_shorter_id) > 0) return DELIVERY_FAIL;
        return strategic_two(current_id, best_shorter_id, target_id, visited, depth + 1);
    }

    // --- Step 2: Try best neighbor from 'other' class ---
    int best_other_id = -1;
    long double max_prob_other = -1.0L;
    for (int neighbor_id : classes.other) {
        if (neighbor_id == prev_id) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        // [FIXED] Use 'current_id' as the incoming node for the DRP lookup.
        long double prob = getP(neighbor_id, current_id, h, d);

        if (prob > max_prob_other) {
            max_prob_other = prob;
            best_other_id = neighbor_id;
        }
    }
    if (best_other_id != -1 && !hasFaultyLink(current_id, best_other_id)) {
        if (visited.count(best_other_id) > 0) return DELIVERY_FAIL;
        return strategic_two(current_id, best_other_id, target_id, visited, depth + 1);
    }

    visited.erase(current_id);
    return DELIVERY_FAIL;
}


int metaTorus::route_strategic_three(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    return strategic_three(-1, start_id, target_id, visited, 0);
}

int metaTorus::strategic_three(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) {
    if (current_id == target_id) return depth;
    visited[current_id] = true;

    StrategicNeighborClasses3 classes = classify_strategic_three(current_id, target_id);

    // --- Step 1: Try 'shorter' class ---
    int best_id = -1; long double max_prob = -1.0L;
    for (int neighbor_id : classes.shorter) {
        if (neighbor_id == prev_id) continue;
        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        // [FIXED] Use 'current_id' as the incoming node for the DRP lookup.
        long double prob = getP(neighbor_id, current_id, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id)) {
        if (visited.count(best_id) > 0) return DELIVERY_FAIL;
        return strategic_three(current_id, best_id, target_id, visited, depth + 1);
    }

    // --- Step 2: Try 'same' class ---
    best_id = -1; max_prob = -1.0L;
    for (int neighbor_id : classes.same) {
        if (neighbor_id == prev_id) continue;
        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        // [FIXED] Use 'current_id' as the incoming node for the DRP lookup.
        long double prob = getP(neighbor_id, current_id, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id)) {
        if (visited.count(best_id) > 0) return DELIVERY_FAIL;
        return strategic_three(current_id, best_id, target_id, visited, depth + 1);
    }

    // --- Step 3: Try 'longer' class ---
    best_id = -1; max_prob = -1.0L;
    for (int neighbor_id : classes.longer) {
        if (neighbor_id == prev_id) continue;
        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        // [FIXED] Use 'current_id' as the incoming node for the DRP lookup.
        long double prob = getP(neighbor_id, current_id, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id)) {
        if (visited.count(best_id) > 0) return DELIVERY_FAIL;
        return strategic_three(current_id, best_id, target_id, visited, depth + 1);
    }

    return DELIVERY_FAIL;
}

std::vector<int> metaTorus::route_strategic_three_with_path(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    return strategic_three_with_path_worker(-1, start_id, target_id, visited);
}

std::vector<int> metaTorus::strategic_three_with_path_worker(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited) {
    if (current_id == target_id) {
        return {current_id};
    }
    visited[current_id] = true;

    StrategicNeighborClasses3 classes = classify_strategic_three(current_id, target_id);

    auto try_candidate_class = [&](const std::vector<int>& candidates) -> std::vector<int> {
        int best_id = -1;
        long double max_prob = -1.0L;

        for (int neighbor_id : candidates) {
            if (neighbor_id == prev_id) continue;
            int h = hammingDistance(neighbor_id, target_id);
            int d = distance(neighbor_id, target_id);
            long double prob = getP(neighbor_id, current_id, h, d);
            if (prob > max_prob) {
                max_prob = prob;
                best_id = neighbor_id;
            }
        }

        if (best_id != -1 && !hasFaultyLink(current_id, best_id) && visited.count(best_id) == 0) {
            std::vector<int> result_path = strategic_three_with_path_worker(current_id, best_id, target_id, visited);
            if (!result_path.empty()) {
                result_path.insert(result_path.begin(), current_id);
                return result_path;
            }
        }
        return {}; // Return empty vector if this class yielded no path
    };

    std::vector<int> path = try_candidate_class(classes.shorter);
    if (!path.empty()) return path;

    path = try_candidate_class(classes.same);
    if (!path.empty()) return path;

    path = try_candidate_class(classes.longer);
    if (!path.empty()) return path;

    return {}; // Return empty vector if no path found
}

std::vector<int> metaTorus::route_strategic_two_with_path(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    return strategic_two_with_path_worker(-1, start_id, target_id, visited);
}

std::vector<int> metaTorus::strategic_two_with_path_worker(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited) {
    if (current_id == target_id) {
        return {current_id};
    }
    visited[current_id] = true;

    StrategicNeighborClasses2 classes = classify_strategic_two(current_id, target_id);

    auto try_candidate_class = [&](const std::vector<int>& candidates) -> std::vector<int> {
        int best_id = -1;
        long double max_prob = -1.0L;

        for (int neighbor_id : candidates) {
            if (neighbor_id == prev_id) continue;
            int h = hammingDistance(neighbor_id, target_id);
            int d = distance(neighbor_id, target_id);
            long double prob = getP(neighbor_id, current_id, h, d);
            if (prob > max_prob) {
                max_prob = prob;
                best_id = neighbor_id;
            }
        }

        if (best_id != -1 && !hasFaultyLink(current_id, best_id) && visited.count(best_id) == 0) {
            std::vector<int> result_path = strategic_two_with_path_worker(current_id, best_id, target_id, visited);
            if (!result_path.empty()) {
                result_path.insert(result_path.begin(), current_id);
                return result_path;
            }
        }
        return {};
    };

    std::vector<int> path = try_candidate_class(classes.shorter);
    if (!path.empty()) return path;

    path = try_candidate_class(classes.other);
    if (!path.empty()) return path;

    return {};
}