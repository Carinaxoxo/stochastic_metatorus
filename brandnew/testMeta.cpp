#include "testMeta.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <cmath>
#include <random>
#include <utility>
#include <algorithm>
#include <queue> // Needed for BFS

// Helper function for safe integer exponentiation.
long long ipow(int base, int exp) {
    long long res = 1;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

// --- Part 1: Metatorus Construction ---
testMeta::testMeta(int n_dims, int k_arity) : n(n_dims), k(k_arity) {
    if (n < 2 || k < 2) {
        throw std::invalid_argument("Metatorus requires n>=2 and k>=2.");
    }
    V = ipow(k, n);
    nodes.resize(V);

// Generate neighbors using the Meta Torus specific logic
    createMetaTorusNeighborsDOR();

    // Initialize probability storage
    int degree = 4; // Metatorus-like degree
    int diameter = n * (k / 2); // Torus-based diameter for calculation
    directed_probabilities.resize(V);
    for (int i = 0; i < V; ++i) {
        directed_probabilities[i].resize(degree);
        for (int j = 0; j < degree; ++j) {
            directed_probabilities[i][j].resize(n + 1);
            for (int h = 0; h <= n; ++h) {
                directed_probabilities[i][j][h].resize(diameter + 1, -1.0L);
            }
        }
    }
}
/**
 * @brief Creates the 4 Meta Torus neighbors for each node, ensuring they are
 * stored in dimension order to support deterministic routing.
 * This implementation is based on your createNeighbors4() logic.
 */
void testMeta::createMetaTorusNeighborsDOR() {
    for (int i = 0; i < V; ++i) {
        std::vector<int> current_coords = getCoords(i);
        nodes[i].neighbors.assign(4, -1);

        // A temporary structure to hold {dimension, neighbor_id} pairs.
        std::vector<std::pair<int, int>> temp_neighbors;

        int header = current_coords[0];
        // The header selects another dimension. We'll use modulo to keep it in bounds.
        int x = (header + 1) % n;
        // Ensure x is not dimension 0, as per the implied logic.
        if (x == 0) x = 1;

        // 1. Neighbor on header ring (dimension 0) in positive direction
        std::vector<int> temp_coords_1 = current_coords;
        temp_coords_1[0] = (header + 1) % (n - 1); // Use k for arity of the ring
        temp_neighbors.push_back({0, getId(temp_coords_1)});

        // 2. Neighbor on header ring (dimension 0) in negative direction
        std::vector<int> temp_coords_2 = current_coords;
        temp_coords_2[0] = (header - 1 + (n - 1)) % (n - 1);
        temp_neighbors.push_back({0, getId(temp_coords_2)});

        // 3. Neighbor in selected dimension 'x' in positive direction
        std::vector<int> temp_coords_3 = current_coords;
        temp_coords_3[x] = (current_coords[x] + 1) % k;
        temp_neighbors.push_back({x, getId(temp_coords_3)});

        // 4. Neighbor in selected dimension 'x' in negative direction
        std::vector<int> temp_coords_4 = current_coords;
        temp_coords_4[x] = (current_coords[x] - 1 + k) % k;
        temp_neighbors.push_back({x, getId(temp_coords_4)});

        // CRITICAL STEP: Sort the pairs. This sorts by dimension number first,
        // which enforces the dimension-order needed for deterministic routing.
        std::sort(temp_neighbors.begin(), temp_neighbors.end());

        // Assign the sorted neighbor IDs to the final neighbors list.
        for (int j = 0; j < 4; ++j) {
            nodes[i].neighbors[j] = temp_neighbors[j].second;
        }
    }
}
// --- Part 2: Directed Routing Probability Calculation (Hybrid Model) ---
void testMeta::calcDirectedRoutingProbabilities() {
    const int full_torus_degree = 2 * n;
    const int diameter = n * (k / 2);
    const double denominator = (full_torus_degree > 1) ? (full_torus_degree - 1.0) : 1.0;

    // Base Case (d=1, h=1)
    for (int i = 0; i < V; ++i) {
        for (int b_id : nodes[i].neighbors) {
            double sum_alpha = 0.0;
            for (int a_id : nodes[i].neighbors) {
                if (a_id == b_id) continue;
                sum_alpha += !hasFaultyLink(i, a_id) ? 1.0 : 0.0;
            }
            double p_1_1 = sum_alpha / denominator;
            setDirectedProbability(i, b_id, 1, 1, p_1_1);
        }
    }

    // Main Dynamic Programming Loop
    for (int d = 2; d <= diameter; ++d) {
        for (int h = 1; h <= std::min(d, n); ++h) {
            for (int i = 0; i < V; ++i) {
                double S_hd = 0.0;
                for (int a_id : nodes[i].neighbors) {
                    if (hasFaultyLink(i, a_id)) continue;

                    long double p1 = getDirectedProbability(a_id, i, h - 1, d - 1);
                    if(p1 < 0) p1 = 0.0L;
                    long double p2 = getDirectedProbability(a_id, i, h, d - 1);
                    if(p2 < 0) p2 = 0.0L;

                    if (h == 1) S_hd += p2;
                    else if (h == d) S_hd += p1;
                    else S_hd += ((h - 1.0) * p1 + (d - h) * p2) / (d - 1.0);
                }

                for (int b_id : nodes[i].neighbors) {
                    double term_for_b = 0.0;
                    if (!hasFaultyLink(i, b_id)) {
                        long double p1 = getDirectedProbability(b_id, i, h - 1, d - 1);
                        if(p1 < 0) p1 = 0.0L;
                        long double p2 = getDirectedProbability(b_id, i, h, d - 1);
                        if(p2 < 0) p2 = 0.0L;

                        if (h == 1) term_for_b = p2;
                        else if (h == d) term_for_b = p1;
                        else term_for_b = ((h - 1.0) * p1 + (d - h) * p2) / (d - 1.0);
                    }
                    double p_hd = (S_hd - term_for_b) / denominator;
                    setDirectedProbability(i, b_id, h, d, p_hd);
                }
            }
        }
    }
}

// --- Part 3: Printing the Results Table ---
void testMeta::printProbabilitiesAsTable(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file '" << filename << "'." << std::endl;
        return;
    }

    const int diameter = n * (k / 2);
    std::vector<std::pair<int, int>> hd_pairs;
    for (int d = 1; d <= diameter; ++d) {
        for (int h = 1; h <= std::min(d, n); ++h) {
            hd_pairs.push_back({h, d});
        }
    }

    for (int node_id = 0; node_id < V; ++node_id) {
        outFile << "\n\n--- Probability Table for Node u = " << node_id << " " << getCoordString(node_id) << " ---\n";
        outFile << std::left << std::setw(25) << "Incoming Neighbor b";
        for (const auto& p : hd_pairs) {
            outFile << "| P(" << p.first << "," << p.second << ")";
        }
        outFile << "\n" << std::string(25, '-');
        for (size_t i = 0; i < hd_pairs.size(); ++i) outFile << "+---------";
        outFile << "\n";

        for (int neighbor_id : nodes[node_id].neighbors) {
            std::string neighbor_str = std::to_string(neighbor_id) + " " + getCoordString(neighbor_id);
            outFile << std::left << std::setw(25) << neighbor_str;
            for (const auto& p : hd_pairs) {
                long double prob = getDirectedProbability(node_id, neighbor_id, p.first, p.second);
                outFile << "| " << std::fixed << std::setprecision(5) << std::setw(7) << (prob < 0 ? 0.0 : prob);
            }
            outFile << "\n";
        }
    }
    outFile.close();
}

// --- Helper Function Implementations ---
std::vector<int> testMeta::getCoords(int id) const {
    std::vector<int> coords(n);
    long long temp_id = id;
    for (int i = n - 1; i >= 0; --i) {
        coords[i] = temp_id % k;
        temp_id /= k;
    }
    return coords;
}

int testMeta::getId(const std::vector<int>& coords) const {
    int id = 0;
    int multiplier = 1;
    for (int i = n - 1; i >= 0; --i) {
        id += coords[i] * multiplier;
        multiplier *= k;
    }
    return id;
}

std::string testMeta::getCoordString(int id) const {
    if (id < 0 || id >= V) return "(invalid)";
    std::vector<int> coords = getCoords(id);
    std::string s = "(";
    for(size_t i = 0; i < coords.size(); ++i) {
        s += std::to_string(coords[i]) + (i == coords.size() - 1 ? "" : ",");
    }
    s += ")";
    return s;
}

void testMeta::setDirectedProbability(int u_id, int b_id, int h, int d, long double prob) {
    if (nodes[u_id].neighbor_id_to_index.count(b_id)) {
        int b_local_index = nodes[u_id].neighbor_id_to_index.at(b_id);
        directed_probabilities[u_id][b_local_index][h][d] = prob;
    }
}

long double testMeta::getDirectedProbability(int u_id, int b_id, int h, int d) const {
    if (nodes[u_id].neighbor_id_to_index.count(b_id)) {
        int b_local_index = nodes[u_id].neighbor_id_to_index.at(b_id);
        if (h >= 0 && (size_t)h < directed_probabilities[u_id][b_local_index].size() &&
            d >= 0 && (size_t)d < directed_probabilities[u_id][b_local_index][h].size()) {
            return directed_probabilities[u_id][b_local_index][h][d];
        }
    }
    return -1.0L;
}

void testMeta::setRandomFaults(double fault_ratio, unsigned int seed) {
    faulty_links.clear();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    for (int i = 0; i < V; ++i) {
        for (int neighbor_id : nodes[i].neighbors) {
            if (i < neighbor_id) {
                if (dist(rng) < fault_ratio) {
                    long long key = static_cast<long long>(i) * V + neighbor_id;
                    faulty_links[key] = true;
                }
            }
        }
    }
}

bool testMeta::hasFaultyLink(int u_id, int v_id) const {
    if (u_id > v_id) std::swap(u_id, v_id);
    long long key = static_cast<long long>(u_id) * V + v_id;
    return faulty_links.count(key);
}

// --- Part 4: Fault-Free Pathfinding using BFS ---
std::vector<int> testMeta::findFaultFreePathBFS(int start_node, int end_node) const {
    // Handle the edge case where start and end are the same.
    if (start_node == end_node) {
        return {start_node};
    }

    // A queue for the nodes to visit.
    std::queue<int> q;
    // A vector to track visited nodes.
    std::vector<bool> visited(V, false);
    // A vector to reconstruct the path by storing parent pointers.
    std::vector<int> parent(V, -1);

    // Initialize BFS
    q.push(start_node);
    visited[start_node] = true;

    int current_node = -1;
    bool path_found = false;

    while (!q.empty()) {
        current_node = q.front();
        q.pop();

        // If we reached the destination, we can stop.
        if (current_node == end_node) {
            path_found = true;
            break;
        }

        // Explore neighbors
        for (int neighbor_id : nodes[current_node].neighbors) {
            // CRITICAL CHECK: Ensure the link is not faulty and the neighbor hasn't been visited.
            if (!hasFaultyLink(current_node, neighbor_id) && !visited[neighbor_id]) {
                visited[neighbor_id] = true;
                parent[neighbor_id] = current_node; // Set parent for path reconstruction
                q.push(neighbor_id);
            }
        }
    }

    // --- Path Reconstruction ---
    std::vector<int> path;
    if (path_found) {
        int at = end_node;
        // Backtrack from the end node to the start node using the parent map.
        while (at != -1) {
            path.push_back(at);
            at = parent[at];
        }
        // The path is constructed backwards, so reverse it.
        std::reverse(path.begin(), path.end());
    }

    // If a path was found, return it. Otherwise, return an empty vector.
    return path;
}

// --- Part 5: Brute Force "Dummy" Adaptive Routing Algorithm ---

/**
 * Calculates the Manhattan distance between two nodes on the torus.
 */
int testMeta::distance(int node_a_id, int node_b_id) const {
    std::vector<int> coords_a = getCoords(node_a_id);
    std::vector<int> coords_b = getCoords(node_b_id);
    int d = 0;
    for (int i = 0; i < n; ++i) {
        int diff = std::abs(coords_a[i] - coords_b[i]);
        d += std::min(diff, k - diff);
    }
    return d;
}

/**
 * Checks if a neighbor is a "preferred" neighbor (i.e., gets one step closer).
 */
bool testMeta::isPreferred(int neighbor_id, int current_id, int target_id) const {
    return distance(neighbor_id, target_id) == distance(current_id, target_id) - 1;
}

NeighborClasses testMeta::classifyNeighbors(int current_id, int target_id) const {
    NeighborClasses classes;
    for (int neighbor_id : nodes.at(current_id).neighbors) {
        if (isPreferred(neighbor_id, current_id, target_id)) {
            classes.preferred.push_back(neighbor_id);
        } else {
            classes.spare.push_back(neighbor_id);
        }
    }
    return classes;
}


/**
 * Public wrapper to start the brute force routing.
 */
int testMeta::route_brute(int start_id, int target_id) const {
    std::unordered_map<int, bool> visited;
    return brute(-1, start_id, target_id, visited, 0);
}

int testMeta::brute(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) const {
    if (current_id == target_id) return depth;
    visited[current_id] = true;

    NeighborClasses classes = classifyNeighbors(current_id, target_id);

    for (int neighbor_id : classes.preferred) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    for (int neighbor_id : classes.spare) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    return DELIVERY_FAIL;
}