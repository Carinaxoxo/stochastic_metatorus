#include "metatorus.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <map>
#include <string>
#include <iomanip>
#include <sstream>

using namespace std;

// A simple struct to hold the statistics for each algorithm
struct AlgoStats {
    int success_count = 0;
    long long total_path_length = 0;

    void recordSuccess(int length) {
        if (length != DELIVERY_FAIL) {
            success_count++;
            total_path_length += length;
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
};


// The main experiment function
void runExperiment(int n, int k, int trials_per_ratio, std::ostream& out) {
    out << "\n\n--- Starting Experiment for Torus n=" << n << ", k=" << k << " ---" << std::endl;
    out << "--- Trials per Fault Ratio: " << trials_per_ratio << " ---" << std::endl;

    metaTorus torus(n, k);

    std::map<double, std::map<std::string, AlgoStats>> all_results;
    std::map<double, double> avg_bfs_lengths;

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uni(0, torus.V - 1);

    for (double fault_ratio = 0.0; fault_ratio <= 0.5; fault_ratio += 0.05) {
        std::cout << "Running n=" << n << ", k=" << k << " at Fault Ratio: " << std::fixed << std::setprecision(2) << fault_ratio << "..." << std::endl;

        long long total_bfs_for_ratio = 0;
        int successful_trials_for_ratio = 0; // Keep track of valid trials

        for (int i = 0; i < trials_per_ratio; ++i) {
            int start_node, target_node, bfs_len;

            while (true) {
                // 1. Set random faults
                torus.setRandomFaults(fault_ratio, rd());

                // 2. Pick random start and target nodes
                start_node = uni(rng);
                target_node = uni(rng);

                // 3. Try to find a path using BFS
                bfs_len = torus.bfs(start_node, target_node);

                // 4. If a path exists, break the loop and use this trial
                if (bfs_len != DELIVERY_FAIL) {
                    cout << "found path of length " << bfs_len << " from " << start_node <<
                    " to " << target_node << " with fault ratio " << fault_ratio << "";
                    break;
                }
                // If no path, the while-loop repeats, re-rolling everything
            }

            successful_trials_for_ratio++;

            // The pair is now guaranteed to be connected.
//            bfs_len = torus.bfs(start_node, target_node);
            total_bfs_for_ratio += bfs_len;

            torus.computeDirectedRoutingProbabilities();

            all_results[fault_ratio]["brute"].recordSuccess(torus.route_brute(start_node, target_node));
            all_results[fault_ratio]["directed"].recordSuccess(torus.route_directed(start_node, target_node));
            all_results[fault_ratio]["strategic_two"].recordSuccess(torus.route_strategic_two(start_node, target_node));
            all_results[fault_ratio]["strategic_three"].recordSuccess(torus.route_strategic_three(start_node, target_node));

            if ((i + 1) % (trials_per_ratio / 10) == 0 && trials_per_ratio >= 10) {
                std::cout << "  ...completed trial " << (i + 1) << "/" << trials_per_ratio << std::endl;
            }
        }

        if (successful_trials_for_ratio > 0) {
            avg_bfs_lengths[fault_ratio] = static_cast<double>(total_bfs_for_ratio) / successful_trials_for_ratio;
        } else {
            avg_bfs_lengths[fault_ratio] = 0;
        }

    }

    // --- Print Final Results Summary Table to the file ---
    out << "\n--- Final Results for n=" << n << ", k=" << k << " ---" << std::endl;
    out << "------------------------------------------------------------------------------------------------------------" << std::endl;
    out << std::left
        << std::setw(10) << "Fault %"
        << "| " << std::setw(12) << "Avg BFS Len"
        << "| " << std::setw(21) << "Brute"
        << "| " << std::setw(21) << "Directed (Simple)"
        << "| " << std::setw(21) << "Strategic (2-Class)"
        << "| " << std::setw(21) << "Strategic (3-Class)" << std::endl;
    out << std::setw(10) << ""
        << "| " << std::setw(12) << ""
        << "| " << std::setw(10) << "Succ. %" << std::setw(11) << "Avg. Len"
        << "| " << std::setw(10) << "Succ. %" << std::setw(11) << "Avg. Len"
        << "| " << std::setw(10) << "Succ. %" << std::setw(11) << "Avg. Len"
        << "| " << std::setw(10) << "Succ. %" << std::setw(11) << "Avg. Len" << std::endl;
    out << "------------------------------------------------------------------------------------------------------------" << std::endl;

    for (auto const& [ratio, results] : all_results) {
        const auto& brute_stats = results.at("brute");
        const auto& directed_stats = results.at("directed");
        const auto& strat2_stats = results.at("strategic_two");
        const auto& strat3_stats = results.at("strategic_three");

        // Find the total number of successful trials for this ratio to calculate success rate correctly
        int total_valid_trials = brute_stats.success_count + directed_stats.success_count + strat2_stats.success_count + strat3_stats.success_count;
        // This is not quite right, let's find a better way. The number of trials for each algo is the same.
        // Let's find the number of trials that were actually run for this ratio.
        int trials_run_for_ratio = 0;
        if (results.count("brute")) { // Check if any stats were recorded
            // This is tricky. Let's assume the number of trials is the denominator for success rate.
            // A better way is to count the successful_trials_for_ratio for each ratio.
            // For now, we will use the trials_per_ratio as a simplification.
            // The logic inside the loop needs to be adjusted to pass this count.
        }


        std::stringstream ratio_ss;
        ratio_ss << std::fixed << std::setprecision(1) << (ratio * 100.0) << " %";

        out << std::fixed << std::setprecision(1) << std::left;
        out << std::setw(10) << ratio_ss.str()
            << "| " << std::setw(12) << avg_bfs_lengths.at(ratio) << "| ";

        out << std::setw(10) << brute_stats.getSuccessRate(trials_per_ratio) * 100.0
            << std::setw(11) << brute_stats.getAveragePathLength() << "| ";

        out << std::setw(10) << directed_stats.getSuccessRate(trials_per_ratio) * 100.0
            << std::setw(11) << directed_stats.getAveragePathLength() << "| ";

        out << std::setw(10) << strat2_stats.getSuccessRate(trials_per_ratio) * 100.0
            << std::setw(11) << strat2_stats.getAveragePathLength() << "| ";

        out << std::setw(10) << strat3_stats.getSuccessRate(trials_per_ratio) * 100.0
            << std::setw(11) << strat3_stats.getAveragePathLength() << std::endl;
    }
    out << "------------------------------------------------------------------------------------------------------------" << std::endl;
}


int main() {
    std::ofstream results_file("experiment_results.txt");
    if (!results_file.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }

    const int TRIALS = 100;

    runExperiment(4, 3, TRIALS, results_file);
    runExperiment(5, 4, TRIALS, results_file);

    results_file.close();
    std::cout << "\nExperiment complete. Results have been saved to experiment_results.txt" << std::endl;

    return 0;
}
