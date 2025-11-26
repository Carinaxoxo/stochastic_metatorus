#include "testMeta.h"
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <iomanip>

struct ExperimentStats {
    int brute_successes = 0;
    long long total_bfs_hops = 0;
    long long total_brute_hops = 0;
};

void run_large_scale_experiment(int n, int k, int num_iterations, double fault_ratio) {
    std::mt19937 master_rng(std::random_device{}());
    ExperimentStats stats;
    int successful_trials = 0;
    long long attempts = 0;
    testMeta network(n, k);

    while (successful_trials < num_iterations) {
        attempts++;
        unsigned int trial_seed = master_rng();
        network.setRandomFaults(fault_ratio, trial_seed);

        std::mt19937 trial_rng(trial_seed);
        std::uniform_int_distribution<int> node_dist(0, network.getV() - 1);
        int start_node = node_dist(trial_rng);
        int end_node = node_dist(trial_rng);

        std::vector<int> bfs_path = network.findFaultFreePathBFS(start_node, end_node);

        if (!bfs_path.empty()) {
            successful_trials++;

            int brute_length = network.route_brute(start_node, end_node);

            stats.total_bfs_hops += (bfs_path.size() - 1);

            if (brute_length != testMeta::DELIVERY_FAIL) {
                stats.brute_successes++;
                stats.total_brute_hops += brute_length;
            }
        }

        if (attempts > num_iterations * 1000 && successful_trials < 10) {
            std::cout << " (Warning: Could not find enough connectable pairs.)";
            break;
        }
    }

    double brute_success_rate = (successful_trials > 0) ? (static_cast<double>(stats.brute_successes) / successful_trials) * 100.0 : 0.0;
    double avg_bfs_len = (successful_trials > 0) ? static_cast<double>(stats.total_bfs_hops) / successful_trials : 0.0;
    double avg_brute_len = (stats.brute_successes > 0) ? static_cast<double>(stats.total_brute_hops) / stats.brute_successes : 0.0;

    std::cout << std::fixed << std::setprecision(2)
              << std::setw(12) << std::left << fault_ratio
              << std::setw(20) << std::right << brute_success_rate << "%"
              << std::setw(18) << avg_bfs_len
              << std::setw(18) << avg_brute_len << std::endl;
}

int main() {
    const int NUM_ITERATIONS = 1000;
    std::vector<std::pair<int, int>> configs = {{4, 3}, {5, 4}};

    for (const auto& config : configs) {
        int n = config.first;
        int k = config.second;

        std::cout << "\n\n====================================================================\n";
        std::cout << "   STARTING EXPERIMENT: " << NUM_ITERATIONS << " iterations for n=" << n << ", k=" << k;
        std::cout << "\n====================================================================\n";
        std::cout << std::left
                  << std::setw(12) << "Fault Ratio"
                  << std::right
                  << std::setw(20) << "Brute Success Rate"
                  << std::setw(18) << "Avg. BFS Length"
                  << std::setw(18) << "Avg. Brute Length" << std::endl;
        std::cout << "--------------------------------------------------------------------\n";

        for (double ratio = 0.0; ratio <= 0.5; ratio += 0.05) {
            run_large_scale_experiment(n, k, NUM_ITERATIONS, ratio);
        }
    }
    std::cout << "\n--- All experiments complete. ---" << std::endl;
    return 0;
}

