#include "metatorus.h"
#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <string>
#include <iomanip>

// A simple struct to hold the statistics for each algorithm
struct AlgoStats {
    int success_count = 0;
    long long total_path_length = 0;
    long long total_bfs_length = 0; // To store the optimal path length for comparison

    void recordSuccess(int length, int bfs_length) {
        if (length != DELIVERY_FAIL) {
            success_count++;
            total_path_length += length;
            total_bfs_length += bfs_length;
        }
    }

    double getSuccessRate(int total_trials) const {
        if (total_trials == 0) return 0.0;
        return static_cast<double>(success_count) / total_trials;
    }

    double getAveragePathLength() const {
        if (success_count == 0) return 0.0;
        return static_cast<double>(total_path_length) / success_count;
    }

    double getAverageBfsLength() const {
        if (success_count == 0) return 0.0;
        return static_cast<double>(total_bfs_length) / success_count;
    }
};

// The main experiment function
void runExperiment() {
    const int TRIALS_PER_RATIO = 1000;
    metaTorus torus(4, 3);

    // Data structure to hold all results
    std::map<double, std::map<std::string, AlgoStats>> all_results;

    // Setup random number generation
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, torus.V - 1);

    // --- Outer Loop: Iterate through fault ratios ---
    for (double fault_ratio = 0.0; fault_ratio <= 0.5; fault_ratio += 0.05) {
        std::cout << "\n--- Starting Experiment for Fault Ratio: " << std::fixed << std::setprecision(2) << fault_ratio << " ---" << std::endl;

        // --- Inner Loop: Run 1000 trials for the current fault ratio ---
        for (int i = 0; i < TRIALS_PER_RATIO; ++i) {
            int start_node, target_node, bfs_len;

            // Loop until we find a pair of nodes that are connected
            do {
                torus.setRandomFaults(fault_ratio, rd());
                start_node = uni(rng);
                do {
                    target_node = uni(rng);
                } while (start_node == target_node);
                bfs_len = torus.bfs(start_node, target_node);
            } while (bfs_len == DELIVERY_FAIL);

            torus.computeDirectedRoutingProbabilities();

            // Run each algorithm and record the results, including the BFS length
            int brute_len = torus.route_brute(start_node, target_node);
            all_results[fault_ratio]["brute"].recordSuccess(brute_len, bfs_len);

            int directed_len = torus.route_directed(start_node, target_node);
            all_results[fault_ratio]["directed"].recordSuccess(directed_len, bfs_len);

            int strat2_len = torus.route_strategic_two(start_node, target_node);
            all_results[fault_ratio]["strategic_two"].recordSuccess(strat2_len, bfs_len);

            int strat3_len = torus.route_strategic_three(start_node, target_node);
            all_results[fault_ratio]["strategic_three"].recordSuccess(strat3_len, bfs_len);

            if ((i + 1) % 100 == 0) {
                std::cout << "Completed trial " << (i + 1) << "/" << TRIALS_PER_RATIO << std::endl;
            }
        }
    }

    // --- Print Final Results Summary Table ---
    std::cout << "\n\n--- Experiment Complete: Final Results ---" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(10) << "Fault %"
              << "| " << std::setw(25) << "Brute"
              << "| " << std::setw(25) << "Directed (Simple)"
              << "| " << std::setw(25) << "Strategic (2-Class)"
              << "| " << std::setw(25) << "Strategic (3-Class)" << std::endl;
    std::cout << std::setw(10) << ""
              << "| " << std::setw(8) << "Succ. %" << std::setw(9) << "Avg. Len" << std::setw(8) << "BFS Len"
              << "| " << std::setw(8) << "Succ. %" << std::setw(9) << "Avg. Len" << std::setw(8) << "BFS Len"
              << "| " << std::setw(8) << "Succ. %" << std::setw(9) << "Avg. Len" << std::setw(8) << "BFS Len"
              << "| " << std::setw(8) << "Succ. %" << std::setw(9) << "Avg. Len" << std::setw(8) << "BFS Len" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;

    for (auto const& [ratio, results] : all_results) {
        const auto& brute_stats = results.at("brute");
        const auto& directed_stats = results.at("directed");
        const auto& strat2_stats = results.at("strategic_two");
        const auto& strat3_stats = results.at("strategic_three");

        std::cout << std::fixed << std::setprecision(1) << std::left;
        std::cout << std::setw(9) << (ratio * 100.0) << "%" << "| ";

        std::cout << std::setw(8) << brute_stats.getSuccessRate(TRIALS_PER_RATIO) * 100.0
                  << std::setw(9) << brute_stats.getAveragePathLength()
                  << std::setw(8) << brute_stats.getAverageBfsLength() << "| ";

        std::cout << std::setw(8) << directed_stats.getSuccessRate(TRIALS_PER_RATIO) * 100.0
                  << std::setw(9) << directed_stats.getAveragePathLength()
                  << std::setw(8) << directed_stats.getAverageBfsLength() << "| ";

        std::cout << std::setw(8) << strat2_stats.getSuccessRate(TRIALS_PER_RATIO) * 100.0
                  << std::setw(9) << strat2_stats.getAveragePathLength()
                  << std::setw(8) << strat2_stats.getAverageBfsLength() << "| ";

        std::cout << std::setw(8) << strat3_stats.getSuccessRate(TRIALS_PER_RATIO) * 100.0
                  << std::setw(9) << strat3_stats.getAveragePathLength()
                  << std::setw(8) << strat3_stats.getAverageBfsLength() << std::endl;
    }
    std::cout << "------------------------------------------------------------------------------------------------------------------------" << std::endl;
}


int main() {
    runExperiment();
    return 0;
}
