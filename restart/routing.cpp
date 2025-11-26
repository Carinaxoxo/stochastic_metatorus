//
// Created by Carina Z on 2025/11/19.
//

#include "metatorus.h"
#include <iostream>
#include <queue>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <cassert> // For assert
#include <random>  // For random number generation
#include <fstream> // For file I/O
#include <iomanip> // For stream manipulators (setw, etc.)
#include <sstream> // For string streams

std::vector<int> metaTorus::BFS(int start_node_id, int end_node_id) const {
    // If start and end are the same, the path is just the node itself.
    if (start_node_id == end_node_id) {
        return {start_node_id};
    }

    // A queue for the nodes to visit.
    std::queue<int> q;
    // A vector to store the parent of each node in the path, to reconstruct it later.
    std::vector<int> parent(V, -1);
    // A vector to keep track of visited nodes to avoid cycles.
    std::vector<bool> visited(V, false);

    // Start the search from the start_node_id.
    q.push(start_node_id);
    visited[start_node_id] = true;

    while (!q.empty()) {
        int current_id = q.front();
        q.pop();

        // If we've reached the destination, we can stop the search.
        if (current_id == end_node_id) {
            break;
        }

        const Node& current_node = nodes[current_id];
        // Check all neighbors of the current node.
        for (size_t i = 0; i < current_node.neighbors.size(); ++i) {
            int neighbor_id = current_node.neighbors[i];
            int link_is_active = current_node.link_status[i];

            // CRITICAL CHECK: Only proceed if the link is NOT faulty and the neighbor hasn't been visited.
            if (link_is_active == 1 && !visited[neighbor_id]) {
                visited[neighbor_id] = true;
                parent[neighbor_id] = current_id;
                q.push(neighbor_id);
            }
        }
    }

    // --- Path Reconstruction ---
    std::vector<int> path;
    // If a path was found (i.e., the end node was visited).
    if (visited[end_node_id]) {
        int current = end_node_id;
        // Backtrack from the end node to the start node using the parent pointers.
        while (current != -1) {
            path.push_back(current);
            current = parent[current];
        }
        // The path is constructed backwards, so we need to reverse it.
        std::reverse(path.begin(), path.end());
    }

    // If no path was found, 'path' will be empty. Otherwise, it holds the shortest path.
    return path;
}

int metaTorus::findShortestPathBFS(int start_node_id, int end_node_id) const {
    // This function now returns the path length (size - 1)
    std::vector<int> path_vector = BFS(start_node_id, end_node_id);
    if(path_vector.empty()){
        return DELIVERY_FAIL;
    }
    return path_vector.size() -1;
}

/*
int metaTorus::findPathStrategic(int start_node_id, int end_node_id) const {
    std::unordered_map<int, bool> visited;
    return RecursiveStrategic(-1, start_node_id, end_node_id, visited, 0);
}

int metaTorus::RecursiveStrategic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth) const {
    if (current_id == target_id) return depth;

    visited[current_id] = true;
    StrategicNeighborClasses3 classes = classifyNeighborsStrategic3(current_id, target_id);

    for (int neighbor_id : classes.shorter) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }
    for (int neighbor_id : classes.same) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }
    for (int neighbor_id : classes.longer) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1);
        }
    }
    return DELIVERY_FAIL;
}
 */

int metaTorus::findPathStrategic(int start_node_id, int end_node_id, bool verbose) const {
    std::unordered_map<int, bool> visited;
    if (verbose) {
        std::cout << "\n>>> STRATEGIC TRACE START >>>" << std::endl;
        std::cout << "Route: " << getCoordString(start_node_id) << " --> " << getCoordString(end_node_id) << std::endl;
    }
    return RecursiveStrategic(-1, start_node_id, end_node_id, visited, 0, verbose);
}

int metaTorus::RecursiveStrategic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const {
    // 1. Base Case: Reached Target
    if (current_id == target_id) {
        if (verbose) std::cout << "    [Goal Reached at Depth " << depth << "]" << std::endl;
        return depth;
    }

    visited[current_id] = true;

    // 2. Classify Neighbors
    StrategicNeighborClasses3 classes = classifyNeighborsStrategic3(current_id, target_id);

    // --- VERBOSE LOGGING BLOCK ---
    if (verbose) {
        int current_dist = getStrategicDistance(current_id, target_id);
        std::cout << "\n[Depth " << depth << "] At " << std::setw(12) << std::left << getCoordString(current_id)
                  << " (Dist to Dest: " << current_dist << ")" << std::endl;

        // Helper lambda to print a list of neighbors
        auto print_list = [&](const std::string& label, const std::vector<int>& list) {
            std::cout << "  " << std::setw(8) << label << ": ";
            if (list.empty()) { std::cout << "None"; }
            for (int nid : list) {
                int d = getStrategicDistance(nid, target_id);
                std::string faulty = hasFaultyLink(current_id, nid) ? "[X]" : "[safe]";
                bool visited_already = visited.count(nid);
                std::string status = visited_already ? "(Vis)" : "";

                std::cout << getCoordString(nid) << "[d:" << d << "]" << faulty << status << " ";
            }
            std::cout << std::endl;
        };

        print_list("Shorter", classes.shorter);
        print_list("Same", classes.same);
        print_list("Longer", classes.longer);
    }
    // -----------------------------

    // 3. Iterate Shorter
    for (int neighbor_id : classes.shorter) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (!visited.count((neighbor_id))) {
                if (verbose) std::cout << "    -> Moving to Shorter: " << getCoordString(neighbor_id) << std::endl;
                return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1, verbose);
            } else {
                return DELIVERY_FAIL;
            }
        }
    }

    // 4. Iterate Same (Misrouting)
    for (int neighbor_id : classes.same) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (!visited.count((neighbor_id))) {
                if (verbose) std::cout << "    -> Moving to Same: " << getCoordString(neighbor_id) << std::endl;
                return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1, verbose);
            } else {
                return DELIVERY_FAIL;
            }
        }
    }

    // 5. Iterate Longer (Backtracking / Detour)
    for (int neighbor_id : classes.longer) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (!visited.count((neighbor_id))) {
                if (verbose) std::cout << "    -> Moving to Longer: " << getCoordString(neighbor_id) << std::endl;
                return RecursiveStrategic(current_id, neighbor_id, target_id, visited, depth + 1, verbose);
            } else {
                return DELIVERY_FAIL;
            }
        }
    }

    return DELIVERY_FAIL;
}

int metaTorus::findPathDirectedStrategic2(int start_id, int target_id, bool verbose) const {
    std::unordered_map<int, bool> visited;
    if (verbose) {
        std::cout << "\n>>> STRATEGIC TRACE (Probabilistic - 2 Categories) >>> "
                  << getCoordString(start_id) << " --> " << getCoordString(target_id) << std::endl;
    }
    return directedRecursiveStrategic2(-1, start_id, target_id, visited, 0, verbose);
}

// --- Recursive Routing Logic (2 Categories) ---
int metaTorus::directedRecursiveStrategic2(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const {
    // 1. Base Case: Reached Target
    if (current_id == target_id) {
        if (verbose) std::cout << "  >> Destination reached at depth " << depth << "!" << std::endl;
        return depth;
    }

    visited[current_id] = true;

    // 2. Verbose Header
    if (verbose) {
        std::cout << "\n  [Step " << depth << "] At " << std::setw(14) << std::left << getCoordString(current_id)
                  << " (Dist to Dest: " << getStrategicDistance(current_id, target_id) << ")" << std::endl;
    }

    // 3. Classify Neighbors (2 Categories)
    StrategicNeighborClasses2 classes = classifyNeighborsStrategic2(current_id, target_id);

    // Helper lambda (Same as before, just reused)
    auto pick_best = [&](const std::vector<int>& neighbors, const std::string& label) -> int {
        int best_id = -1;
        long double max_prob = -1.0L;

        if (verbose && !neighbors.empty()) {
            std::cout << "    --- Checking " << label << " Neighbors ---" << std::endl;
            std::cout << "    "
                      << std::setw(16) << "Neighbor"
                      << std::setw(10) << "Link"
                      << std::setw(8)  << "Dist"
                      << std::setw(12) << "Prob(%)"
                      << "Notes" << std::endl;
        }

        for (int neighbor_id : neighbors) {
            bool is_faulty = hasFaultyLink(current_id, neighbor_id);
            bool is_prev = (neighbor_id == prev_id);
            bool is_visited = visited.count(neighbor_id);

            int d_neighbor = getStrategicDistance(neighbor_id, target_id);
            int h_neighbor = getHammingDistance(neighbor_id, target_id);
            int d_lookup = std::min(d_neighbor, diameter);

            long double prob = -1.0L;
            if (!is_faulty) {
                prob = getDirectedProbability(neighbor_id, current_id, h_neighbor, d_lookup);
            }

            if (verbose) {
                std::cout << "    " << std::setw(16) << getCoordString(neighbor_id);
                if (is_faulty) std::cout << std::setw(10) << "[FAULTY]";
                else           std::cout << std::setw(10) << "[ OK ]";
                std::cout << std::setw(8) << d_neighbor;
                if (is_faulty) std::cout << std::setw(12) << "---";
                else           std::cout << std::setw(12) << std::fixed << std::setprecision(4) << prob;
                std::string note = "";
                if (is_prev) note += "(Prev) ";
                if (is_visited) note += "(Vis) ";
                std::cout << note << std::endl;
            }

            if (is_faulty) continue;
            if (is_prev) continue;

            if (prob > max_prob) {
                max_prob = prob;
                best_id = neighbor_id;
            }
        }
        return best_id;
    };

    // 4. Try Shorter Neighbors
    int next_node = pick_best(classes.shorter, "Shorter");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Shorter: " << getCoordString(next_node) << std::endl;
        return directedRecursiveStrategic2(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    // 5. Try Other Neighbors (Same + Longer combined)
    next_node = pick_best(classes.other, "Other");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Other: " << getCoordString(next_node) << std::endl;
        return directedRecursiveStrategic2(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    if (verbose) std::cout << "    >> DEAD END (No valid moves found)." << std::endl;
    return DELIVERY_FAIL;
}

int metaTorus::findPathDirectedStrategic3(int start_id, int target_id, bool verbose) const {
    std::unordered_map<int, bool> visited;
    if (verbose) {
        std::cout << "\n>>> STRATEGIC TRACE (Probabilistic) >>> " << getCoordString(start_id) << " --> " << getCoordString(target_id) << std::endl;
    }
    return directedRecursiveStrategic3(-1, start_id, target_id, visited, 0, verbose);
}

int metaTorus::directedRecursiveStrategic3(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const {
    if (current_id == target_id) {
        if (verbose) std::cout << "  >> Destination reached at depth " << depth << "!" << std::endl;
        return depth;
    }

    visited[current_id] = true;

    if (verbose) {
        std::cout << "\n  [Step " << depth << "] At " << getCoordString(current_id)
                  << " (Dist: " << getStrategicDistance(current_id, target_id) << ")" << std::endl;
    }

    StrategicNeighborClasses3 classes = classifyNeighborsStrategic3(current_id, target_id);
    int dist_ct = getDistance(current_id, target_id);

    // Helper lambda to find best neighbor in a list based on probability
    auto pick_best = [&](const std::vector<int>& neighbors, const std::string& label) -> int {
        int best_id = -1;
        long double max_prob = -1.0L;

//        if (verbose && !neighbors.empty()) std::cout << "    --- Checking " << label << " ---" << std::endl;
        if (verbose && !neighbors.empty()) {
            std::cout << "    --- Checking " << label << " Neighbors ---" << std::endl;
            // Print Table Header
            std::cout << "    "
                      << std::setw(16) << "Neighbor"
                      << std::setw(10) << "Link"
                      << std::setw(8)  << "Dist"
                      << std::setw(12) << "Prob(%)"
                      << "Notes" << std::endl;
        }

        for (int neighbor_id : neighbors) {
            // 1. Gather Information
            bool is_faulty = hasFaultyLink(current_id, neighbor_id);
            bool is_prev = (neighbor_id == prev_id);
            bool is_visited = visited.count(neighbor_id);
            int d_neighbor = getStrategicDistance(neighbor_id, target_id);
            int h_neighbor = getHammingDistance(neighbor_id, target_id);
            // Clamp 'd' to diameter to prevent index out of bounds in lookup
            int d_lookup = std::min(d_neighbor, diameter);

            bool faulty = hasFaultyLink(current_id, neighbor_id);

            // Calculate params for probability lookup
            int h = getHammingDistance(neighbor_id, target_id);
            int d = std::min(getStrategicDistance(neighbor_id, target_id), diameter);

            // If link is faulty, effectively prob is 0.0, but we skip strictly
            long double prob = getDirectedProbability(neighbor_id, current_id, h, d);

            if (verbose) {
                std::cout << "    " << std::setw(16) << getCoordString(neighbor_id);

                // Link Status
                if (is_faulty) std::cout << std::setw(10) << "[FAULTY]";
                else           std::cout << std::setw(10) << "[ OK ]";

                // Distance
                std::cout << std::setw(8) << d_neighbor;

                // Probability
                if (is_faulty) std::cout << std::setw(12) << "---";
                else           std::cout << std::setw(12) << std::fixed << std::setprecision(4) << prob;

                // Notes (Prev / Visited)
                std::string note = "";
                if (is_prev) note += "(Prev) ";
                if (is_visited) note += "(Vis) ";
                std::cout << note << std::endl;
            }

            if (neighbor_id == prev_id) continue; // Don't go back immediately
            if (faulty) continue;

            if (prob > max_prob) {
                max_prob = prob;
                best_id = neighbor_id;
            }
        }
        return best_id;
    };


    // 1. Try Shorter
    int next_node = pick_best(classes.shorter, "Shorter");
    if (next_node != -1) {
        if (visited.count(next_node)){
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Selected Shorter: " << getCoordString(next_node) << std::endl;
        return directedRecursiveStrategic3(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    // 2. Try Same
    next_node = pick_best(classes.same, "Same");
    if (next_node != -1) {
        if (visited.count(next_node)){
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Selected Same: " << getCoordString(next_node) << std::endl;
        return directedRecursiveStrategic3(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    // 3. Try Longer
    next_node = pick_best(classes.longer, "Longer");
    if (next_node != -1) {
        if (visited.count(next_node)){
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Selected Longer: " << getCoordString(next_node) << std::endl;
        return directedRecursiveStrategic3(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    if (verbose) std::cout << "    >> DEAD END." << std::endl;
    return DELIVERY_FAIL;

    /*
    int best_shorter_id = -1;
    long double max_prob_shorter = -1.0L;

    for (int neighbor_id : classes.shorter) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            int h = getHammingDistance(neighbor_id, target_id);
            long double prob = getDirectedProbability(neighbor_id, current_id, h, std::min(getStrategicDistance(neighbor_id, target_id), diameter));
            if (prob > max_prob_shorter) {
                max_prob_shorter = prob;
                best_shorter_id = neighbor_id;
            }
        }
    }

    if (best_shorter_id != -1 && max_prob_shorter > 0) {
        if (visited.count(best_shorter_id) != 0) {
            return DELIVERY_FAIL;
        }
        return directedRecursiveStrategic3(current_id, best_shorter_id, target_id, visited, depth + 1);
    }

    int best_same_id = -1;
    long double max_prob_same = -1.0L;

    for (int neighbor_id : classes.same) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            int h = getHammingDistance(neighbor_id, target_id);
            long double prob = getDirectedProbability(neighbor_id, current_id, h, std::min(getStrategicDistance(neighbor_id, target_id), diameter));
            if (prob > max_prob_same) {
                max_prob_same = prob;
                best_same_id = neighbor_id;
            }
        }
    }

    if (best_same_id != -1) {
        if (visited.count(best_same_id) != 0) {
            return DELIVERY_FAIL;
        }
        return directedRecursiveStrategic3(current_id, best_same_id, target_id, visited, depth + 1);
    }


    int best_longer_id = -1;
    long double max_prob_longer = -1.0L;

    for (int neighbor_id : classes.longer) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            int h = getHammingDistance(neighbor_id, target_id);
            long double prob = getDirectedProbability(neighbor_id, current_id, h, std::min(getStrategicDistance(neighbor_id, target_id), diameter));
            if (prob > max_prob_longer) {
                max_prob_longer = prob;
                best_longer_id = neighbor_id;
            }
        }
    }

    if (best_longer_id != -1) {
        if (visited.count(best_longer_id) != 0) {
            return DELIVERY_FAIL;
        }
        return directedRecursiveStrategic3(current_id, best_longer_id, target_id, visited, depth + 1);
    }

    return DELIVERY_FAIL;
    */
}

// In routing.cpp

int metaTorus::findPathBrute(int start_node_id, int end_node_id, bool verbose) const {
    std::unordered_map<int, bool> visited;
    if (verbose) {
        std::cout << "\n>>> BRUTE FORCE TRACE (Preferred/Spare) >>> "
                  << getCoordString(start_node_id) << " --> " << getCoordString(end_node_id) << std::endl;
    }
    return bruteRecursive(-1, start_node_id, end_node_id, visited, 0, verbose);
}

int metaTorus::bruteRecursive(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const {
    // 1. Base Case
    if (current_id == target_id) {
        if (verbose) std::cout << "  >> Destination reached at depth " << depth << "!" << std::endl;
        return depth;
    }

    visited[current_id] = true;

    // 2. Verbose Header
    if (verbose) {
        std::cout << "\n  [Step " << depth << "] At " << std::setw(14) << std::left << getCoordString(current_id)
                  << " (Manhattan Dist: " << getDistance(current_id, target_id) << ")" << std::endl;
    }

    // 3. Classify (Preferred vs Spare)
    NeighborClasses classes = classifyNeighbors(current_id, target_id);

    for (int neighbor_id : classes.preferred) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return bruteRecursive(current_id, neighbor_id, target_id, visited, depth + 1, verbose);
        }
    }
    for (int neighbor_id : classes.spare) {
        if (neighbor_id != prev_id && !hasFaultyLink(current_id, neighbor_id)) {
            if (visited.count(neighbor_id)) return DELIVERY_FAIL;
            return bruteRecursive(current_id, neighbor_id, target_id, visited, depth + 1, verbose);
        }
    }
    return DELIVERY_FAIL;

    // Helper lambda to pick the FIRST valid neighbor
    /*
    auto pick_first_valid = [&](const std::vector<int>& neighbors, const std::string& label) -> int {
        if (verbose && !neighbors.empty()) {
            std::cout << "    --- Checking " << label << " Neighbors ---" << std::endl;
            std::cout << "    "
                      << std::setw(16) << "Neighbor"
                      << std::setw(10) << "Link"
                      << std::setw(8)  << "Dist"
                      << "Notes" << std::endl;
        }

        int selected_id = -1;

        for (int neighbor_id : neighbors) {
            bool is_faulty = hasFaultyLink(current_id, neighbor_id);
            bool is_prev = (neighbor_id == prev_id);
            bool is_visited = visited.count(neighbor_id);
            int d_neighbor = getDistance(neighbor_id, target_id); // using Manhattan Distance

            if (verbose) {
                std::cout << "    " << std::setw(16) << getCoordString(neighbor_id);

                if (is_faulty) std::cout << std::setw(10) << "[FAULTY]";
                else           std::cout << std::setw(10) << "[ OK ]";

                std::cout << std::setw(8) << d_neighbor;

                std::string note = "";
                if (is_prev) note += "(Prev) ";
                if (is_visited) note += "(Vis) ";

                if (selected_id == -1 && !is_faulty && !is_prev) {
                    note += " <-- SELECTED";
                }
                std::cout << note << std::endl;
            }

            // Greedy selection: Take the first valid one
            if (selected_id == -1 && !is_faulty && !is_prev) {
                selected_id = neighbor_id;
                if (!verbose) return selected_id;
            }
        }
        return selected_id;
    };

    // 4. Try Preferred
    int next_node = pick_first_valid(classes.preferred, "Preferred");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Selected " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Preferred: " << getCoordString(next_node) << std::endl;
        return bruteRecursive(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    // 5. Try Spare
    next_node = pick_first_valid(classes.spare, "Spare");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Selected " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Spare: " << getCoordString(next_node) << std::endl;
        return bruteRecursive(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    if (verbose) std::cout << "    >> DEAD END (No valid moves found)." << std::endl;
    return DELIVERY_FAIL;
     */
}

int metaTorus::findPathBruteProbabilistic(int start_node_id, int end_node_id, bool verbose) const {
    std::unordered_map<int, bool> visited;
    if (verbose) {
        std::cout << "\n>>> BRUTE PROBABILISTIC TRACE (Preferred/Spare + Max Prob) >>> "
                  << getCoordString(start_node_id) << " --> " << getCoordString(end_node_id) << std::endl;
    }
    return bruteRecursiveProbabilistic(-1, start_node_id, end_node_id, visited, 0, verbose);
}

int metaTorus::bruteRecursiveProbabilistic(int prev_id, int current_id, int target_id, std::unordered_map<int, bool>& visited, int depth, bool verbose) const {
    // 1. Base Case
    if (current_id == target_id) {
        if (verbose) std::cout << "  >> Destination reached at depth " << depth << "!" << std::endl;
        return depth;
    }

    visited[current_id] = true;

    // 2. Verbose Header
    if (verbose) {
        std::cout << "\n  [Step " << depth << "] At " << std::setw(14) << std::left << getCoordString(current_id)
                  << " (Manhattan Dist: " << getDistance(current_id, target_id) << ")" << std::endl;
    }

    // 3. Classify (Using Standard Brute Force Logic: Preferred vs Spare)
    NeighborClasses classes = classifyNeighbors(current_id, target_id);
    int dist_ct = getDistance(current_id, target_id);


    // Helper lambda to pick the neighbor with HIGHEST PROBABILITY
    auto pick_best_prob = [&](const std::vector<int>& neighbors, const std::string& label) -> int {
        int best_id = -1;
        long double max_prob = -1.0L;

        if (verbose && !neighbors.empty()) {
            std::cout << "    --- Checking " << label << " Neighbors ---" << std::endl;
            std::cout << "    "
                      << std::setw(16) << "Neighbor"
                      << std::setw(10) << "Link"
                      << std::setw(8)  << "Dist"
                      << std::setw(12) << "Prob(%)"
                      << "Notes" << std::endl;
        }

        for (int neighbor_id : neighbors) {
            // 1. Gather Info
            bool is_faulty = hasFaultyLink(current_id, neighbor_id);
            bool is_prev = (neighbor_id == prev_id);
            bool is_visited = visited.count(neighbor_id);

            // We must calculate Strategic/Hamming distance specifically for the probability lookup
            // even though we used Manhattan for classification.
            int d_strat = getDistance(neighbor_id, target_id);
            int h_val = getHammingDistance(neighbor_id, target_id);
            int d_lookup = std::min(d_strat, diameter);

            long double prob = -1.0L;
            if (!is_faulty) {
                prob = getDirectedProbability(neighbor_id, current_id, h_val, d_lookup);
            }

            // 2. Verbose Output
            if (verbose) {
                std::cout << "    " << std::setw(16) << getCoordString(neighbor_id);
                if (is_faulty) std::cout << std::setw(10) << "[FAULTY]";
                else           std::cout << std::setw(10) << "[ OK ]";

                std::cout << std::setw(8) << d_strat; // Showing strategic dist here as it relates to prob

                if (is_faulty) std::cout << std::setw(12) << "---";
                else           std::cout << std::setw(12) << std::fixed << std::setprecision(4) << prob;

                std::string note = "";
                if (is_prev) note += "(Prev) ";
                if (is_visited) note += "(Vis) ";
                std::cout << note << std::endl;
            }

            // 3. Selection Logic
            if (is_faulty) continue;
            if (is_prev) continue;

            if (prob > max_prob) {
                max_prob = prob;
                best_id = neighbor_id;
            }
        }
        return best_id;
    };

    // 4. Try Preferred (Pick highest prob)
    int next_node = pick_best_prob(classes.preferred, "Preferred");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Preferred: " << getCoordString(next_node) << std::endl;
        return bruteRecursiveProbabilistic(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    // 5. Try Spare (Pick highest prob)
    next_node = pick_best_prob(classes.spare, "Spare");
    if (next_node != -1) {
        if (visited.count(next_node)) {
            if (verbose) std::cout << "    >> Best option " << getCoordString(next_node) << " is visited. STOP." << std::endl;
            return DELIVERY_FAIL;
        }
        if (verbose) std::cout << "    >> Moving to Spare: " << getCoordString(next_node) << std::endl;
        return bruteRecursiveProbabilistic(current_id, next_node, target_id, visited, depth + 1, verbose);
    }

    if (verbose) std::cout << "    >> DEAD END (No valid moves found)." << std::endl;
    return DELIVERY_FAIL;
}
