//
// Created by Carina Z on 2025/11/19.
//

#include "metatorus.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <random>
#include <set>

long long ipow(int base, int exp) {
    long long res = 1;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

metaTorus::metaTorus(int dims, int k_arity) : dimensions(dims), k(k_arity) {
    n = dimensions - 1;

    if (n < 3 || k < 2) {
        throw std::invalid_argument("Metatorus requires n>=3 and k>=2.");
    }
    V = n * ipow(k, n);
    //dimensions = n + 1
//    V = ipow(k, n) - ipow(k, n-1); // 0 <= h <= n-1
    nodes.resize(V);

    for (long long i = 0; i < V; ++i) {
        nodes[i].id = (int)i;
        nodes[i].value.resize(dimensions); // Correctly resizes to 4
        long long temp_id = i;

        // Handle data dimensions (v_0, v_1, v_2), which use arity k
        for (int j = dimensions - 1; j >= 1; --j) { // j = 3, 2, 1
            nodes[i].value[j] = temp_id % k;
            temp_id /= k;
        }

        // Handle header (h), which uses arity (n-1)
        int header_arity = n;
        nodes[i].value[0] = temp_id % header_arity; // j = 0
    }

    createMetaTorusNeighborsDOR();

    // Seed the random number generator
    rand_generator.seed(std::random_device()());
    // After creating neighbors, initialize all links to be active (status = 1)
    for (long long i = 0; i < V; ++i) {
        nodes[i].link_status.assign(nodes[i].neighbors.size(), 1);
    }

    // Initialize probability storage
    int degree = 4;
    diameter = static_cast<int>(floor(k / 2.0)) * dimensions;
    directed_probabilities.resize(V);
    for (int i = 0; i < V; ++i) {
        directed_probabilities[i].resize(degree);
        for (int j = 0; j < degree; ++j) {
            directed_probabilities[i][j].resize(dimensions + 1);
            for (int h = 0; h <= dimensions; ++h) {
                directed_probabilities[i][j][h].resize(diameter + 1, -1.0L); // Initialize with -1
            }
        }
    }
}

int metaTorus::getId(const std::vector<int>& coord) const {
    long long id = 0;

// Calculate ID
    id = coord[0]; // Start with header
    for (int j = 1; j < dimensions; ++j) {
        id = id * k + coord[j]; // Add data dimensions
    }

    assert(id >= 0 && id < V); // This will now pass (0 <= id < 81)
    return (int)id;
}

std::string metaTorus::getCoordString(int id) const {
    if (id < 0 || id >= V) return "Invalid";
    std::stringstream ss;
    ss << "(";
    for (size_t i = 0; i < nodes[id].value.size(); ++i) {
        ss << nodes[id].value[i];
        if (i < nodes[id].value.size() - 1) ss << ",";
    }
    ss << ")";
    return ss.str();
}

void metaTorus::createMetaTorusNeighborsDOR() {
    for (int i = 0; i < V; ++i) {
        std::vector<int> current_coords = nodes[i].value;
        nodes[i].neighbors.assign(4, -1);

        int header = current_coords[0];

        // x is the specific dimension this node is currently "pointing" to
        int x = header + 1;

        std::vector<int> temp_coords_1 = current_coords;
        temp_coords_1[0] = (header + 1) % n;
        nodes[i].neighbors[0] = getId(temp_coords_1);
//        temp_neighbors.push_back({0, getId(temp_coords_1)});
        std::vector<int> temp_coords_2 = current_coords;
        temp_coords_2[0] = (header - 1 + n) % n;
        nodes[i].neighbors[1] = getId(temp_coords_2);
//        temp_neighbors.push_back({0, getId(temp_coords_2)});
        std::vector<int> temp_coords_3 = current_coords;
        temp_coords_3[x] = (current_coords[x] + 1) % k;
        nodes[i].neighbors[2] = getId(temp_coords_3);
//        temp_neighbors.push_back({x, getId(temp_coords_3)});
        std::vector<int> temp_coords_4 = current_coords;
        temp_coords_4[x] = (current_coords[x] - 1 + k) % k;
        nodes[i].neighbors[3] = getId(temp_coords_4);
//        temp_neighbors.push_back({x, getId(temp_coords_4)});
    }
}

int metaTorus::getDistance(int node1_id, int node2_id) const {
    if (node1_id < 0 || node1_id >= V || node2_id < 0 || node2_id >= V) {
        return -1; // Invalid ID
    }

    const auto& coords1 = nodes[node1_id].value;
    const auto& coords2 = nodes[node2_id].value;
    int distance = 0;

    // 1. Handle the header dimension (i=0) with its arity (n-1)
    int header_arity = n;
    int diff_h = std::abs(coords1[0] - coords2[0]);
    distance += std::min(diff_h, header_arity - diff_h);

    // 2. Handle the data dimensions (i=1 to n-1) with arity k
    for (int i = 1; i < dimensions; ++i) {
        int diff = std::abs(coords1[i] - coords2[i]);
        // In a torus, the distance is the shorter of the two ways around
        distance += std::min(diff, k - diff);
    }
    return distance;
}

int metaTorus::getHammingDistance(int node1_id, int node2_id) const {
    if (node1_id < 0 || node1_id >= V || node2_id < 0 || node2_id >= V) {
        return -1; // Invalid ID
    }
    const auto& coords1 = nodes[node1_id].value;
    const auto& coords2 = nodes[node2_id].value;
    int h_dist = 0;
    for(int i = 0; i < dimensions; ++i) {
        if (coords1[i] != coords2[i]) {
            h_dist++;
        }
    }
    return h_dist;
}


void metaTorus::setFaultyLinks(double fault_rate) {
    assert(fault_rate >= 0.0 && fault_rate <= 1.0);

    // 1. Create a list of all unique links
    std::vector<std::pair<int, int>> all_links;
    for (long long i = 0; i < V; ++i) {
        for (int neighbor_id : nodes[i].neighbors) {
            if (i < neighbor_id) {
                all_links.push_back({(int)i, neighbor_id});
            }
        }
    }

    // 2. Shuffle the list of links randomly
    std::shuffle(all_links.begin(), all_links.end(), rand_generator);

    // 3. Calculate the exact number of links to make faulty
    int total_links = all_links.size();
    int faults_to_create = static_cast<int>(round(total_links * fault_rate));

    // 4. Mark the first 'faults_to_create' links from the shuffled list as faulty
    for (int i = 0; i < faults_to_create; ++i) {
        int node1_id = all_links[i].first;
        int node2_id = all_links[i].second;

        // Find and mark the forward link
        for (size_t j = 0; j < nodes[node1_id].neighbors.size(); ++j) {
            if (nodes[node1_id].neighbors[j] == node2_id) {
                nodes[node1_id].link_status[j] = 0;
                break;
            }
        }
        // Find and mark the reverse link
        for (size_t k = 0; k < nodes[node2_id].neighbors.size(); ++k) {
            if (nodes[node2_id].neighbors[k] == node1_id) {
                nodes[node2_id].link_status[k] = 0;
                break;
            }
        }
    }
}

std::vector<std::pair<int, int>> metaTorus::getFaultyLinks() const {
    std::vector<std::pair<int, int>> faulty_list;
    for (const auto& node : nodes) {
        for (size_t i = 0; i < node.neighbors.size(); ++i) {
            if (node.link_status[i] == 0) {
                faulty_list.push_back({node.id, node.neighbors[i]});
            }
        }
    }
    return faulty_list;
}

bool metaTorus::hasFaultyLink(int u, int v) const {
    if (u < 0 || u >= V) return true;

    const Node& node = nodes[u];
    for (size_t i = 0; i < node.neighbors.size(); ++i) {
        if (node.neighbors[i] == v) {
            // Return true if status is 0 (faulty), false if 1 (active)
            return (node.link_status[i] == 0);
        }
    }
    return true; // Assume faulty if they aren't neighbors
}

void metaTorus::printNetwork() const {
    std::cout << "\n--- Metatorus Network Details (dimensions= " << dimensions << ", n=" << n << ", k=" << k << ", V=" << V << ") ---\n";
    for (const auto& node : nodes) {
        std::cout << "Node " << std::setw(4) << node.id << " Coords(";
        for (size_t i = 0; i < node.value.size(); ++i) {
            std::cout << node.value[i] << (i == node.value.size() - 1 ? "" : ",");
        }
        std::cout << ") \tNeighbors: [";

        for (size_t i = 0; i < node.neighbors.size(); ++i) {
            int neighbor_id = node.neighbors[i];
            int status = node.link_status[i];

            // Check the status of the link
            if (status == 1) {
                const Node& neighbor_node = nodes[neighbor_id];
                std::cout << std::setw(4) << neighbor_node.id << " (";
                for (size_t j = 0; j < neighbor_node.value.size(); ++j) {
                    std::cout << neighbor_node.value[j] << (j == neighbor_node.value.size() - 1 ? "" : ",");
                }
                std::cout << ")";
            } else {
                // Print a clear indicator for a faulty link (status == 0)
                std::cout << "   [FAULTY]   ";
            }
            if (i < node.neighbors.size() - 1) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "---------------------------------------------------------------------------\n" << std::endl;
}

int metaTorus::getSpecificNeighbor(int nodeIndex, int neighborIndex) const{
    // Check if the neighborIndex is valid
    assert(0 <= neighborIndex && neighborIndex < V);
    // Iterate through the neighbors of node `a` to find the index of the specific neighbor
    for (int i = 0; i < 4; i++) {
        if (nodes[nodeIndex].neighbors[i] == neighborIndex) {
            return i; // Return the index if the neighbor is found
        }
    }

    // If no matching neighbor is found, return -1
    return -1;
}

long double metaTorus::getDirectedProbability(int u_id, int b_id, int h, int d) const {
    assert(0 <= h && h <= dimensions);
    assert(0 <= d && d <= diameter);
//
//    // Base case from your logic, destination is reached.
    if (h == 0) {
        assert(d == 0);
        return 1.0;// Base case: self-loop probability
    }
//
//    // Find the local index for the incoming neighbor 'b_id'
    int b_local_index = getSpecificNeighbor(u_id, b_id);
//    assert(b_local_index >= 0 && b_local_index < 4);
//
//    assert(directed_probabilities[u_id][b_local_index][h][d] != -1);
//
//    return directed_probabilities[u_id][b_local_index][h][d];

    // Check if the value was calculated
    long double prob = directed_probabilities[u_id][b_local_index][h][d];
    if (prob == -1.0L) {
        std::cerr << "!! DEBUG WARNING [getDirectedProbability]: Probability P(h=" << h << ", d=" << d << ") for u=" << u_id
                  << ", b=" << b_id << " (local_idx=" << b_local_index << ") was NOT CALCULATED (-1.0)." << std::endl;
        std::cerr << "  ->   Node (u_id): " << u_id << " (Coords: " << getCoordString(u_id) << ")" << std::endl;
        std::cerr << "  ->   From (b_id): " << b_id << " (Coords: " << getCoordString(b_id) << ")" << std::endl;
        std::cerr << "  ->   Passed h = " << h << ", d = " << d << std::endl;
    }

    return prob;
}

void metaTorus::setDirectedProbability(int u_id, int b_id, int h, int d, long double prob) {
    int b_local_index = getSpecificNeighbor(u_id, b_id);
    assert(b_local_index >= 0 && b_local_index < 4);

    assert(0 <= h && h <= dimensions);
    assert(0 <= d && d <= diameter);
    assert(0.0 <= prob && prob <= 1.0);

    assert(directed_probabilities[u_id][b_local_index][h][d] == -1.0L);
//    if (directed_probabilities[u_id][b_local_index][h][d] != -1.0L) {
//        // Keep the first computed value — it’s identical in all valid paths.
//        return;
//    }
    directed_probabilities[u_id][b_local_index][h][d] = prob;
}

void metaTorus::calculateDirectedRoutingProbabilities() {
//    std::cout << "[DEBUG] === Starting calculateDirectedRoutingProbabilities ===" << std::endl;
    // This denominator is based on the theoretical degree of a standard n-dimensional torus.
//    const double denominator = 2.0 * n - 1.0;
    const double denominator = 2.0 * dimensions - 1.0;

    // Base case: d=1
//    std::cout << "[DEBUG:d=1] --- Calculating Base Case (d=1) ---" << std::endl;
    for (int u = 0; u < V; ++u) { // For each node u (by ID 'i')
//         std::cout << "[DEBUG]   u = " << u << std::endl;
        for (int b_id : nodes[u].neighbors) { // For each incoming neighbor b
//             std::cout << "[DEBUG]       b_id = " << b_id << std::endl;
            double sum_alpha = 0.0;
            for (int a_id : nodes[u].neighbors) { // For each outgoing neighbor a
                if (a_id == b_id) continue;
                double alpha = 0.0;
                alpha = !hasFaultyLink(u, a_id) ? 1.0L : 0.0;
//                 std::cout << "[DEBUG]       a_id = " << a_id << ", alpha = " << alpha << std::endl;
                sum_alpha += alpha;
            }
            double p_1_1 = sum_alpha / denominator;
//             std::cout << "[DEBUG]     -> sum_alpha = " << sum_alpha << ", p_1_1 = " << p_1_1 << std::endl;
            setDirectedProbability(u, b_id, 1, 1, p_1_1);
        }
    }
//    std::cout << "[DEBUG] --- Base Case (d=1) Complete ---" << std::endl;
    //base case all good

    // Dynamic programming step for d > 1
//    std::cout << "[DEBUG:d>1] --- Calculating Case (d > 1) ---" << std::endl;
    for (int d = 2; d <= diameter; ++d) {
//        std::cout << "[DEBUG] --- Starting DP Step (d=" << d << ") ---" << std::endl;
        for (int h = 1; h <= std::min(dimensions, d); ++h) {
//            std::cout << "[DEBUG]     --- h=" << h << " ---" << std::endl;
            for (int  u = 0; u < V; ++u) { // For each node u (by ID 'i')
                // std::cout << "[DEBUG]       u = " << u << std::endl;
                long double sum = 0.0;

//                 std::cout << "[DEBUG]         --- Sum loop (over 'a') ---" << std::endl;
                for (int a_id : nodes[u].neighbors) { // Sum over all potential next hops 'a'
                    double alpha = !hasFaultyLink(u, a_id) ? 1.0L : 0.0;
//                    std::cout << "[DEBUG]           a_id = " << a_id << ", alpha = " << alpha << std::endl;
//                    std::cout << "[DEBUG]           -> h - 1 = " << h - 1 << ", d - 1 = " << d - 1 << std::endl;
//                    long double p1 = getDirectedProbability(a_id, u, h - 1, d - 1);
//                    std::cout << "[DEBUG]           -> h = " << h << ", d - 1 = " << d - 1 << std::endl;
//                    long double p2 = getDirectedProbability(a_id, u, h, d - 1);
                    long double term_to_add = 0.0;
                    if (h == 1) term_to_add = alpha * getDirectedProbability(a_id, u, h, d - 1);
                    else if (h == d) term_to_add = alpha * getDirectedProbability(a_id, u, h - 1, d - 1);
                    else {
                        long double p1 = getDirectedProbability(a_id, u, h - 1, d - 1);
                        long double p2 = getDirectedProbability(a_id, u, h, d - 1);
//                        std::cout << "[DEBUG]           -> p1(h-1, d-1) = " << p1 << ", p2(h, d-1) = " << p2 << std::endl;
                        term_to_add = alpha * (((h - 1.0L) * p1 + (d - h) * p2) / (d - 1.0L));
//                        std::cout << "[DEBUG]           -> adding to sum: " << term_to_add << std::endl;
                    }
                    sum += term_to_add;
//                     if (h == 1) sum += alpha * p2;
//                    else if (h == d) sum += alpha * p1;
//                    else sum += alpha * (((h - 1.0L) * p1 + (d - h) * p2) / (d - 1.0L));
                }
//                 std::cout << "[DEBUG]         --- Total sum = " << sum << " ---" << std::endl;

//                 std::cout << "[DEBUG]         --- Set loop (over 'b') ---" << std::endl;
                for (int b_id : nodes[u].neighbors) { // For each incoming neighbor 'b'
//                     std::cout << "[DEBUG]           b_id = " << b_id << std::endl;
                    long double prob_b = 0.0;
                    double alpha = !hasFaultyLink(u, b_id) ? 1.0L : 0.0;
//                     std::cout << "[DEBUG]           -> alpha(b) = " << alpha << std::endl;

                    long double subtraction_term = 0.0;
                    if (h == 1) {
//                        std::cout << "[DEBUG]           -> h = " << h << ", d - 1 = " << d - 1 << std::endl;
                        subtraction_term = alpha * getDirectedProbability(b_id, u, h, d - 1); // p2
//                        std::cout << "[DEBUG]           -> (h=1) subtraction_term = " << subtraction_term << std::endl;
                    } else if (h == d) {
//                        std::cout << "[DEBUG]           -> h - 1 = " << h - 1 << ", d - 1 = " << d - 1 << std::endl;
                        subtraction_term = alpha * getDirectedProbability(b_id, u, h - 1, d - 1); // p1
//                        std::cout << "[DEBUG]           -> (h=d) subtraction_term = " << subtraction_term << std::endl;
                    } else {
//                        std::cout << "[DEBUG]           -> h - 1 = " << h - 1 << ", d - 1 = " << d - 1 << std::endl;
                        long double p1_b = getDirectedProbability(b_id, u, h - 1, d - 1);
//                        std::cout << "[DEBUG]           -> h = " << h << ", d - 1 = " << d - 1 << std::endl;
                        long double p2_b = getDirectedProbability(b_id, u, h, d - 1);
//                        std::cout << "[DEBUG]           -> (else) p1_b = " << p1_b << ", p2_b = " << p2_b << std::endl;
                        long double inner = ((h - 1.0) * p1_b + (d - h) * p2_b) / (d - 1.0);
                        subtraction_term = alpha * inner;
//                         std::cout << "[DEBUG]           -> (else) inner = " << inner << ", subtraction_term = " << subtraction_term << std::endl;
                    }
                    prob_b = (sum - subtraction_term) / denominator;
//                    if (h == 1) prob_b = (sum - alpha * getDirectedProbability(b_id, u, h , d - 1)) / denominator;
//                    else if (h == d) prob_b = (sum - alpha * getDirectedProbability(b_id, u, h - 1, d - 1)) / denominator;
//                    else {
//                        long double inner = ((h - 1.0) * getDirectedProbability(b_id, u, h - 1, d - 1)
//                                             + (d - h) * getDirectedProbability(b_id, u, h, d - 1))
//                                            / (d - 1.0);
//                        prob_b = (sum - alpha * inner) / denominator;
//                    }
                    setDirectedProbability(u, b_id, h, d, prob_b);
                }
            }
        }
//        std::cout << "[DEBUG] --- Finished DP Step (d=" << d << ") ---" << std::endl;
    }
//    std::cout << "[DEBUG] === calculateDirectedRoutingProbabilities Complete ===" << std::endl;
}

// Helper: Robust Ring Sweep for Line-like traversal on a Ring
int calculateRingSweep(int n, int start_header, int dest_header, const std::vector<int>& dimensions_to_visit) {
    std::set<int> points_set;
    points_set.insert(start_header);
    points_set.insert(dest_header);

    // Add headers required for dimension changes
    // Note: To change dim 'd', we need header '(d-1)'.
    for (int dim : dimensions_to_visit) {
        points_set.insert((dim - 1 + n) % n);
    }

    std::vector<int> p(points_set.begin(), points_set.end());

    // Edge Case: Only 1 point (Start == Dest and no stops)
    if (p.size() <= 1) return 0;

    // 1. Sort points to process segments sequentially
    std::sort(p.begin(), p.end());

    int min_total_dist = 2147483647; // Max int

    // 2. Iterate over every possible "Gap" in the ring.
    // The "Gap" is the segment between p[i] and p[i+1] that we DO NOT traverse.
    // By removing one gap, we unroll the ring into a linear line starting at u, ending at v.
    for (size_t i = 0; i < p.size(); ++i) {

        // The line starts at the point AFTER the gap (u) and goes CW to the point BEFORE the gap (v).
        int u = p[(i + 1) % p.size()];
        int v = p[i];

        // Length of this linear segment (covering all required points)
        int L = (v - u + n) % n;

        // Map Start and Dest to relative positions on this line [0, L]
        int s_pos = (start_header - u + n) % n;
        int d_pos = (dest_header - u + n) % n;

        // We are now on a straight line [0, L].
        // We start at s_pos, must visit 0 and L (the extremities), and end at d_pos.

        // Strategy 1: Go Left to 0, then Right to L, then to Dest
        // Path: s_pos -> 0 -> L -> d_pos
        // Dist: (s_pos - 0) + (L - 0) + |L - d_pos|  (Note: d_pos is <= L)
        int cost1 = s_pos + L + (L - d_pos);

        // Strategy 2: Go Right to L, then Left to 0, then to Dest
        // Path: s_pos -> L -> 0 -> d_pos
        // Dist: (L - s_pos) + (L - 0) + |0 - d_pos|
        int cost2 = (L - s_pos) + L + d_pos;

        int current_cost = std::min(cost1, cost2);

        if (current_cost < min_total_dist) {
            min_total_dist = current_cost;
        }
    }

    return min_total_dist;
}

// --- Public Calculation: Strategic Distance ---

int metaTorus::getStrategicDistance(int node1_id, int node2_id) const {
    if (node1_id == node2_id) return 0;

    const std::vector<int>& c1 = nodes[node1_id].value;
    const std::vector<int>& c2 = nodes[node2_id].value;

    // 1. Data Cost: Simple ring distance for data dimensions (1 to n)
    // See Page 37 of PDF
    int data_cost = 0;
    std::vector<int> diff_dims;

    for (int i = 1; i < dimensions; ++i) {
        if (c1[i] != c2[i]) {
            int diff = std::abs(c1[i] - c2[i]);
            // min( |a-b|, k-|a-b| )
            data_cost += std::min(diff, k - diff);
            diff_dims.push_back(i);
        }
    }

    // 2. Header Cost: Use the Sweep Algorithm
    int header_cost = calculateRingSweep(n, c1[0], c2[0], diff_dims);

    return data_cost + header_cost;
}

NeighborClasses metaTorus::classifyNeighbors(int current_node_id, int dest_node_id) const {
    NeighborClasses result;
    if (current_node_id < 0 || current_node_id >= V || dest_node_id < 0 || dest_node_id >= V) {
        return result; // Return empty result for invalid IDs
    }

    int dist_to_dest = getDistance(current_node_id, dest_node_id);
    const Node& current_node = nodes[current_node_id];

    for (int neighbor_id : current_node.neighbors) {
        if (getDistance(neighbor_id, dest_node_id) < dist_to_dest) {
            result.preferred.push_back(neighbor_id);
        } else {
            result.spare.push_back(neighbor_id);
        }
    }
    return result;
}

StrategicNeighborClasses2 metaTorus::classifyNeighborsStrategic2(int current_node_id, int dest_node_id) const {
    StrategicNeighborClasses2 result;
    if (current_node_id < 0 || current_node_id >= V || dest_node_id < 0 || dest_node_id >= V) {
        return result;
    }

    int dist_to_dest = getStrategicDistance(current_node_id, dest_node_id);
    const Node& current_node = nodes[current_node_id];

    for (int neighbor_id : current_node.neighbors) {
        int neighbor_dist = getStrategicDistance(neighbor_id, dest_node_id);

        // Strictly shorter neighbors go to 'shorter'
        // Everything else (same distance or longer) goes to 'other'
        if (neighbor_dist < dist_to_dest) {
            result.shorter.push_back(neighbor_id);
        } else {
            result.other.push_back(neighbor_id);
        }
    }
    return result;
}

StrategicNeighborClasses3 metaTorus::classifyNeighborsStrategic3(int current_node_id, int dest_node_id) const {
    StrategicNeighborClasses3 result;
    if (current_node_id < 0 || current_node_id >= V || dest_node_id < 0 || dest_node_id >= V) {
        return result; // Return empty result for invalid IDs
    }

    int dist_to_dest = getStrategicDistance(current_node_id, dest_node_id);
    const Node& current_node = nodes[current_node_id];

    for (int neighbor_id : current_node.neighbors) {
        int neighbor_dist = getStrategicDistance(neighbor_id, dest_node_id);
        if (neighbor_dist < dist_to_dest) {
            result.shorter.push_back(neighbor_id);
        } else if (neighbor_dist == dist_to_dest) {
            result.same.push_back(neighbor_id);
        } else {
            result.longer.push_back(neighbor_id);
        }
    }
    return result;
}


