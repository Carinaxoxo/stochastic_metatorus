#include "testMeta.h"
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iomanip>

// Helper function to run a test for a Metatorus configuration
void constructionTest(int dimensions, int arity, const std::string& output_filename) {
    std::cout << "\n========================================================\n";
    std::cout << "--- Metatorus Probability Test for n=" << dimensions << ", k=" << arity << " ---" << std::endl;
    std::cout << "========================================================\n";

    try {
        std::cout << "Constructing Metatorus-like network..." << std::endl;
        testMeta network(dimensions, arity);
        std::cout << "Network constructed with " << network.getV() << " nodes." << std::endl;

        std::cout << "\nCalculating directed routing probabilities..." << std::endl;
        network.calcDirectedRoutingProbabilities();
        std::cout << "Calculation complete." << std::endl;

        network.printProbabilitiesAsTable(output_filename);
        std::cout << "Probability tables have been saved to '" << output_filename << "'." << std::endl;

    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << "Skipping test for n=" << dimensions << ", k=" << arity << "." << std::endl;
    }
}

// Helper to run a single pathfinding test and print results
void runPathfindingTest(int dimensions, int arity, double fault_ratio, unsigned int seed) {
    std::cout << "\n--- Test: n=" << dimensions << ", k=" << arity
              << ", Fault Ratio=" << std::fixed << std::setprecision(2) << fault_ratio
              << " ---\n";

    try {
        testMeta network(dimensions, arity);
        network.setRandomFaults(fault_ratio, seed);

        // Generate a random pair of nodes
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> dist(0, network.getV() - 1);
        int start_node = dist(rng);
        int end_node = dist(rng);

        std::cout << "Finding path from " << network.getCoordString(start_node)
                  << " to " << network.getCoordString(end_node) << "...\n";

        // Find and report the path
        std::vector<int> path = network.findFaultFreePathBFS(start_node, end_node);

        if (path.empty()) {
            std::cout << "Result: No fault-free path found. 😔\n";
        } else {
            std::cout << "Result: Path found! ✅\n  Path: ";
            for (size_t i = 0; i < path.size(); ++i) {
                std::cout << network.getCoordString(path[i]) << (i == path.size() - 1 ? "" : " -> ");
            }
            std::cout << "\n  Length: " << path.size() - 1 << " hops.\n";
        }

    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main() {
    // Define the network configurations to test
    std::vector<std::pair<int, int>> configs = {{4, 3}, {5, 4}};
    unsigned int base_seed = 42;

    for (const auto& config : configs) {
        int n = config.first;
        int k = config.second;

        std::cout << "\n========================================================\n";
        std::cout << "   STARTING SIMULATION FOR n=" << n << ", k=" << k;
        std::cout << "\n========================================================\n";

        // Loop through fault ratios from 0.0 to 0.5
        for (double ratio = 0.0; ratio <= 0.5; ratio += 0.05) {
            // Use a consistent seed for each test run for reproducibility
            runPathfindingTest(n, k, ratio, base_seed);
        }

        base_seed++; // Use a different seed for the next configuration
    }

    std::cout << "\n--- All simulations complete. ---" << std::endl;

    return 0;
}
