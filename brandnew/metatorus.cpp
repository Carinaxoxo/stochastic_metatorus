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
        for (int b = 0; b < 4; ++b)
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

long double metaTorus::getP(int u, int bdir, int h, int d) const {
    return P[u][bdir].get(h, d);
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

    for (int u = 0; u < V; ++u) {
        long double total_sum_over_all_neighbors = 0.0L;
        for (int dir = 0; dir < 4; ++dir) {
            int a = nodes[u].neighbors[dir];
            total_sum_over_all_neighbors += (long double)alpha(u, a);
        }
        for (int bdir = 0; bdir < 4; ++bdir) {
            int b = nodes[u].neighbors[bdir];
            long double numerator = total_sum_over_all_neighbors - (long double)alpha(u, b);
            P[u][bdir].at(1, 1) = numerator / (long double)denom;
        }
    }

    for (int d = 2; d <= Dmax; ++d) {
        int hmax = min(d, Hcap);
        for (int h = 1; h <= hmax; ++h) {
            for (int u = 0; u < V; ++u) {
                long double Shd = 0.0L;
                for (int dir = 0; dir < 4; ++dir) {
                    int a = nodes[u].neighbors[dir];
                    int rev = rev_dir[u][dir];
                    long double term = 0.0L;
                    if (h == 1) term = P[a][rev].get(1, d - 1);
                    else if (h == d) term = P[a][rev].get(h - 1, d - 1);
                    else term = ((long double)(h - 1) / (d - 1)) * P[a][rev].get(h - 1, d - 1)
                                + ((long double)(d - h) / (d - 1)) * P[a][rev].get(h, d - 1);
                    Shd += (long double)alpha(u, a) * term;
                }
                for (int bdir = 0; bdir < 4; ++bdir) {
                    int b = nodes[u].neighbors[bdir];
                    int rev = rev_dir[u][bdir];
                    long double term_b = 0.0L;
                    if (h == 1) term_b = P[b][rev].get(1, d - 1);
                    else if (h == d) term_b = P[b][rev].get(h - 1, d - 1);
                    else term_b = ((long double)(h - 1) / (d - 1)) * P[b][rev].get(h - 1, d - 1)
                                  + ((long double)(d - h) / (d - 1)) * P[b][rev].get(h, d - 1);
                    long double num = Shd - (long double)alpha(u, b) * term_b;
                    P[u][bdir].at(h, d) = num / (long double)denom;
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
        nodes[i].neighbors.clear();
        nodes[i].neighbors.reserve(4);
        int header = current[0] % (n - 1);
        int x = header + 1;
        assert(1 <= x && x <= (n - 1));
        {
            vector<int> temp = current;
            temp[x] = (temp[x] + 1) % k;
            nodes[i].neighbors.push_back(getId(temp));
        }
        {
            vector<int> temp = current;
            temp[x] = (temp[x] - 1 + k) % k;
            nodes[i].neighbors.push_back(getId(temp));
        }
        {
            vector<int> temp = current;
            temp[0] = (header + 1) % (n - 1);
            nodes[i].neighbors.push_back(getId(temp));
        }
        {
            vector<int> temp = current;
            temp[0] = (header - 1 + (n - 1)) % (n - 1);
            nodes[i].neighbors.push_back(getId(temp));
        }
    }
}

void metaTorus::computeReverseDirs() {
    rev_dir.assign(V, array<int,4>{-1, -1, -1, -1});
    for (int u = 0; u < V; ++u) {
        for (int dir = 0; dir < 4; ++dir) {
            int v = nodes[u].neighbors[dir];
            int r = -1;
            for (int t = 0; t < 4; ++t) {
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
        assert((int)nbrs.size() == 4);
        for (int t = 0; t < 4; ++t) {
            assert(0 <= nbrs[t] && nbrs[t] < V);
            assert(nbrs[t] != i);
            for (int u = t + 1; u < 4; ++u) assert(nbrs[t] != nbrs[u]);
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
        for (int dir = 0; dir < 4; ++dir) {
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

void metaTorus::printAllDRP(std::ostream& os) const {
    os.setf(std::ios::fixed);
    os << std::setprecision(6);
    os << "u,b,h,d,P\n";
    for (int u = 0; u < V; ++u) {
        for (int b = 0; b < 4; ++b) {
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
        for (int j = 0; j < 4; ++j) {
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
                for (int index = 0; index < 4; index++) {
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

// --------- Routing and Pathfinding Implementation ---------
/**
 * @brief Public wrapper to start the greedy routing algorithm.
 */
int metaTorus::route_brute(int start_id, int target_id) {
    std::unordered_map<int, bool> visited;
    // Initial call: no previous node (-1), start at depth 0.
    return brute(-1, start_id, target_id, visited, 0);
}

/**
 * @brief The recursive, greedy routing algorithm as you designed it.
 * It commits to the first valid path it finds.
 */
int metaTorus::brute(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth)
{
    if (current_id == target_id) return depth;

    visited[current_id] = true;

    // Get the classification of neighbors (preferred vs. spare)
    NeighborClasses neighbors = classifyNeighbors(current_id, target_id);

    // --- Step 1: Try to route to any available preferred neighbor ---
    for (int neighbor_id : neighbors.preferred) {
        // Don't immediately go back to the node we just came from
        if (neighbor_id == prev_id) continue;

        // Check if the link is fault-free
        if (!hasFaultyLink(current_id, neighbor_id)) {
            // If the neighbor has already been visited in this path, it's a cycle.
            if (visited.count(neighbor_id) > 0) continue; // Just skip, don't fail the whole route

            // Greedily commit to this path and return the result immediately.
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    // --- Step 2: If no preferred neighbors worked, try spare neighbors ---
    for (int neighbor_id : neighbors.spare) {
        if (neighbor_id == prev_id) continue;

        if (!hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id) > 0) continue;

            // Greedily commit to this path and return the result immediately.
            return brute(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }

    // If no preferred or spare neighbors could be taken, this path is a dead end.
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
    const Node& current_node = nodes[current_id];
    NeighborClasses neighbors = classifyNeighbors(current_id, target_id);

    // --- Step 1: Find the BEST candidate among all preferred neighbors ---
    int best_preferred_id = -1;
    long double max_prob_p = -1.0L;

    for (int neighbor_id : neighbors.preferred) {
        if (neighbor_id == prev_id) continue;

        int fwd_dir = -1;
        for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);

        if (prob > max_prob_p) {
            max_prob_p = prob;
            best_preferred_id = neighbor_id;
        }
    }

    // --- Step 2: If a best preferred neighbor was found, try to route through it ---
    if (best_preferred_id != -1 && !hasFaultyLink(current_id, best_preferred_id) && visited.count(best_preferred_id) == 0) {
        int result = directed(current_id, best_preferred_id, target_id, visited, depth + 1);
        // If this path was successful, return the result.
        if (result != DELIVERY_FAIL) return result;
    }

    // --- Step 3: If the preferred path failed or was invalid, find the BEST spare neighbor ---
    int best_spare_id = -1;
    long double max_prob_s = -1.0L;

    for (int neighbor_id : neighbors.spare) {
        if (neighbor_id == prev_id) continue;

        int fwd_dir = -1;
        for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);

        if (prob > max_prob_s) {
            max_prob_s = prob;
            best_spare_id = neighbor_id;
        }
    }

    // --- Step 4: If a best spare neighbor was found, try to route through it ---
    if (best_spare_id != -1 && !hasFaultyLink(current_id, best_spare_id) && visited.count(best_spare_id) == 0) {
        int result = directed(current_id, best_spare_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
    }

    // --- Step 5: If all options have been exhausted and failed, backtrack ---
    // Un-mark the current node so other paths can use it
    visited.erase(current_id);
    return DELIVERY_FAIL;
}


// --------- Strategic Classification Helpers ---------

int metaTorus::ringDistance(int p1, int p2, int size) const {
    int diff = std::abs(p1 - p2);
    return std::min(diff, size - diff);
}

// CORRECTED: This function now uses a simpler, more robust heuristic.
int metaTorus::acrossRing(int h_start, int h_target, const std::vector<int>& need_to_visit) const {
    int ring_size = n - 1;
    if (ring_size <= 0) return 0;

    // Find the min and max dimension that needs to be visited, including start and end
    int min_visit = h_start;
    int max_visit = h_start;
    bool needs_travel = false;

    for (int i = 0; i < ring_size; ++i) {
        if (need_to_visit[i] == 1) {
            min_visit = std::min(min_visit, i);
            max_visit = std::max(max_visit, i);
            needs_travel = true;
        }
    }
    min_visit = std::min(min_visit, h_target);
    max_visit = std::max(max_visit, h_target);

    if (!needs_travel && h_start == h_target) return 0;

    // Calculate the cost of traversing the segment that contains all required nodes
    int forward_dist = max_visit - min_visit;
    int backward_dist = ring_size - forward_dist;

    // The total travel is the shorter of the two paths, plus the distance to cover the segment itself.
    return std::min(ringDistance(h_start, min_visit, ring_size) + forward_dist,
                    ringDistance(h_start, max_visit, ring_size) + forward_dist);
}


// --------- Strategic Classification Implementations ---------

StrategicNeighborClasses2 metaTorus::classify_strategic_two(int c_id, int t_id) {
    const Node& c = nodes[c_id];
    const Node& t = nodes[t_id];
    StrategicNeighborClasses2 classes;

    const auto& standing = c.value;
    const auto& destination = t.value;
    // CORRECTED: The header index must be within the bounds of the coordinate vector.
    int header = (standing[0] + 1) % n;
    if (header == 0) header = 1; // Header must be a dimension from 1 to n-1

    // Dimension move cost
    int d_dim = ringDistance(standing[header], destination[header], k);
    int d1_dim = ringDistance(nodes[c.neighbors[0]].value[header], destination[header], k);
    int d2_dim = ringDistance(nodes[c.neighbors[1]].value[header], destination[header], k);

    if (d1_dim < d_dim) classes.shorter.push_back(c.neighbors[0]);
    else classes.other.push_back(c.neighbors[0]);

    if (d2_dim < d_dim) classes.shorter.push_back(c.neighbors[1]);
    else classes.other.push_back(c.neighbors[1]);

    // Ring move cost
    std::vector<int> need_to_visit(n - 1, 0);
    for(int i = 0; i < (n - 1); i++) {
        if(standing[i + 1] != destination[i + 1]) need_to_visit[i] = 1;
    }

    int ring_cost = acrossRing(standing[0], destination[0], need_to_visit);
    int ring_cost3 = acrossRing(nodes[c.neighbors[2]].value[0], destination[0], need_to_visit);
    int ring_cost4 = acrossRing(nodes[c.neighbors[3]].value[0], destination[0], need_to_visit);

    if (ring_cost3 < ring_cost) classes.shorter.push_back(c.neighbors[2]);
    else classes.other.push_back(c.neighbors[2]);

    if (ring_cost4 < ring_cost) classes.shorter.push_back(c.neighbors[3]);
    else classes.other.push_back(c.neighbors[3]);

    return classes;
}

StrategicNeighborClasses3 metaTorus::classify_strategic_three(int c_id, int t_id) {
    const Node& c = nodes[c_id];
    const Node& t = nodes[t_id];
    StrategicNeighborClasses3 classes;

    const auto& standing = c.value;
    const auto& destination = t.value;
    // CORRECTED: The header index must be within the bounds of the coordinate vector.
    int header = (standing[0] + 1) % n;
    if (header == 0) header = 1; // Header must be a dimension from 1 to n-1

    // Dimension move cost
    int d_dim = ringDistance(standing[header], destination[header], k);
    int d1_dim = ringDistance(nodes[c.neighbors[0]].value[header], destination[header], k);
    int d2_dim = ringDistance(nodes[c.neighbors[1]].value[header], destination[header], k);

    if (d1_dim < d_dim) classes.shorter.push_back(c.neighbors[0]);
    else if (d1_dim == d_dim) classes.same.push_back(c.neighbors[0]);
    else classes.longer.push_back(c.neighbors[0]);

    if (d2_dim < d_dim) classes.shorter.push_back(c.neighbors[1]);
    else if (d2_dim == d_dim) classes.same.push_back(c.neighbors[1]);
    else classes.longer.push_back(c.neighbors[1]);

    // Ring move cost
    std::vector<int> need_to_visit(n - 1, 0);
    for(int i = 0; i < (n - 1); i++) {
        if(standing[i + 1] != destination[i + 1]) need_to_visit[i] = 1;
    }

    int ring_cost = acrossRing(standing[0], destination[0], need_to_visit);
    int ring_cost3 = acrossRing(nodes[c.neighbors[2]].value[0], destination[0], need_to_visit);
    int ring_cost4 = acrossRing(nodes[c.neighbors[3]].value[0], destination[0], need_to_visit);

    if (ring_cost3 < ring_cost) classes.shorter.push_back(c.neighbors[2]);
    else if (ring_cost3 == ring_cost) classes.same.push_back(c.neighbors[2]);
    else classes.longer.push_back(c.neighbors[2]);

    if (ring_cost4 < ring_cost) classes.shorter.push_back(c.neighbors[3]);
    else if (ring_cost4 == ring_cost) classes.same.push_back(c.neighbors[3]);
    else classes.longer.push_back(c.neighbors[3]);

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

    const Node& current_node = nodes[current_id];
    StrategicNeighborClasses2 classes = classify_strategic_two(current_id, target_id);

    // --- Step 1: Try best neighbor from 'shorter' class ---
    int best_shorter_id = -1;
    long double max_prob_shorter = -1.0L;
    for (int neighbor_id : classes.shorter) {
        if (neighbor_id == prev_id) continue;
        int fwd_dir = -1; for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);

        if (prob > max_prob_shorter) {
            max_prob_shorter = prob;
            best_shorter_id = neighbor_id;
        }
    }
    if (best_shorter_id != -1 && !hasFaultyLink(current_id, best_shorter_id) && visited.count(best_shorter_id) == 0) {
        int result = strategic_two(current_id, best_shorter_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
    }

    // --- Step 2: Try best neighbor from 'other' class ---
    int best_other_id = -1;
    long double max_prob_other = -1.0L;
    for (int neighbor_id : classes.other) {
        if (neighbor_id == prev_id) continue;
        int fwd_dir = -1; for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;

        int h = hammingDistance(neighbor_id, target_id);
        int d = distance(neighbor_id, target_id);
        int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);

        if (prob > max_prob_other) {
            max_prob_other = prob;
            best_other_id = neighbor_id;
        }
    }
    if (best_other_id != -1 && !hasFaultyLink(current_id, best_other_id) && visited.count(best_other_id) == 0) {
        int result = strategic_two(current_id, best_other_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
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

    const Node& current_node = nodes[current_id];
    StrategicNeighborClasses3 classes = classify_strategic_three(current_id, target_id);

    // --- Step 1: Try 'shorter' class ---
    int best_id = -1; long double max_prob = -1.0L;
    for (int neighbor_id : classes.shorter) {
        if (neighbor_id == prev_id) continue;
        int fwd_dir = -1; for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;
        int h = hammingDistance(neighbor_id, target_id); int d = distance(neighbor_id, target_id); int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id) && visited.count(best_id) == 0) {
        int result = strategic_three(current_id, best_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
    }

    // --- Step 2: Try 'same' class ---
    best_id = -1; max_prob = -1.0L;
    for (int neighbor_id : classes.same) {
        if (neighbor_id == prev_id) continue;
        int fwd_dir = -1; for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;
        int h = hammingDistance(neighbor_id, target_id); int d = distance(neighbor_id, target_id); int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id) && visited.count(best_id) == 0) {
        int result = strategic_three(current_id, best_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
    }

    // --- Step 3: Try 'longer' class ---
    best_id = -1; max_prob = -1.0L;
    for (int neighbor_id : classes.longer) {
        if (neighbor_id == prev_id) continue;
        int fwd_dir = -1; for(int j=0; j<4; ++j) if(current_node.neighbors[j] == neighbor_id) { fwd_dir = j; break; }
        if (fwd_dir == -1) continue;
        int h = hammingDistance(neighbor_id, target_id); int d = distance(neighbor_id, target_id); int rev = rev_dir[current_id][fwd_dir];
        long double prob = getP(neighbor_id, rev, h, d);
        if (prob > max_prob) { max_prob = prob; best_id = neighbor_id; }
    }
    if (best_id != -1 && !hasFaultyLink(current_id, best_id) && visited.count(best_id) == 0) {
        int result = strategic_three(current_id, best_id, target_id, visited, depth + 1);
        if (result != DELIVERY_FAIL) return result;
    }

    visited.erase(current_id);
    return DELIVERY_FAIL;
}