#include "metatorus.h"
#include "experiments.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <future> // Required for std::async
#include <mutex>  // Required for std::mutex

// Global mutex to prevent jumbled CSV output when running in parallel
std::mutex print_mutex;

// [Add this to experiments.cpp]

void test_recursive_strategic_vs_bfs(int n, int k, int trials, double fault_rate) {
    std::cout << "\n\n============================================================" << std::endl;
    std::cout << "   COMPARISON: BFS vs RecursiveStrategic" << std::endl;
    std::cout << "   (n=" << n << ", k=" << k << ", Faults=" << fault_rate << ", Trials=" << trials << ")" << std::endl;
    std::cout << "============================================================" << std::endl;

    // 1. Setup Network (Single topology for throughput testing)
    metaTorus network(n, k);
    network.setFaultyLinks(fault_rate);
    // Note: We do NOT need calculateDirectedRoutingProbabilities for pure RecursiveStrategic

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> node_dist(0, (int)network.V - 1);

    long long total_bfs_hops = 0;
    long long total_strat_hops = 0;
    int common_successes = 0;
    int bfs_found_pairs = 0;

    // 2. Loop until we have enough valid starting pairs
    while (bfs_found_pairs < trials) {
        int start = node_dist(rng);
        int end = node_dist(rng);

        if (start == end) continue;

        // A. Run BFS (Baseline)
        int bfs_len = network.findShortestPathBFS(start, end);

        if (bfs_len != metaTorus::DELIVERY_FAIL) {
            bfs_found_pairs++;

            // B. Run Recursive Strategic
            // false = turn off verbose printing so console doesn't flood
            int strat_len = network.findPathStrategic(start, end, false);

            if (strat_len != metaTorus::DELIVERY_FAIL) {
                // Only compare lengths if BOTH succeeded
                common_successes++;
                total_bfs_hops += bfs_len;
                total_strat_hops += strat_len;
            }
        }
    }

    // 3. Calculate Stats
    double avg_bfs = (common_successes > 0) ? (double)total_bfs_hops / common_successes : 0.0;
    double avg_strat = (common_successes > 0) ? (double)total_strat_hops / common_successes : 0.0;
    double success_rate = (100.0 * common_successes) / trials;

    // 4. Output Results
    std::cout << "\n---------------- Results ----------------" << std::endl;
    std::cout << "Total Valid Pairs (BFS Found): " << trials << std::endl;
    std::cout << "Strategic Successes:           " << common_successes << std::endl;
    std::cout << "Strategic Success Rate:        " << std::fixed << std::setprecision(2) << success_rate << "%" << std::endl;
    std::cout << "\n--- Path Length Comparison (on " << common_successes << " common successes) ---" << std::endl;
    std::cout << "Total BFS Hops:       " << total_bfs_hops << std::endl;
    std::cout << "Total Strategic Hops: " << total_strat_hops << std::endl;
    std::cout << "Avg BFS Length:       " << avg_bfs << std::endl;
    std::cout << "Avg Strategic Length: " << avg_strat << std::endl;
    std::cout << "Efficiency Ratio:     " << (avg_strat / avg_bfs) << "x BFS length" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
}

void run_visual_trace_experiment(int dimension, int k, int num_traces) {
    std::cout << "\n\n============================================================\n";
    std::cout << "   VISUAL TRACE EXPERIMENT (n=" << dimension << ", k=" << k << ")";
    std::cout << "\n============================================================\n";

    std::mt19937 rng(12345); // Fixed seed for visual traces
    metaTorus network(dimension, k);

    network.setFaultyLinks(0.25);
    network.calculateDirectedRoutingProbabilities();

    int completed = 0;
    while (completed < num_traces) {
        std::uniform_int_distribution<int> node_dist(0, (int)network.V - 1);
        int start = node_dist(rng);
        int end = node_dist(rng);

        if (start == end) continue;

        std::vector<int> bfs_path = network.BFS(start, end);

        if (!bfs_path.empty()) {
            std::cout << "\n------------------------------------------------------------" << std::endl;
            std::cout << "Trial " << (completed + 1) << ": "
                      << network.getCoordString(start) << " to " << network.getCoordString(end) << std::endl;

            // Verbose = true prints the trace
            int hops = network.findPathStrategic(start, end, true);
            std::cout << "\n>>> Result: " << (hops == metaTorus::DELIVERY_FAIL ? "FAILED" : "SUCCESS")
                      << " (Strategy Hops: " << hops << ")" << std::endl;
            completed++;
        }
    }
}

void run_strategy_comparison(int dimension, int k, int num_traces, unsigned int seed) {
    std::cout << "\n\n============================================================\n";
    std::cout << "   STRATEGY COMPARISON (n=" << dimension << ", k=" << k << ")";
    std::cout << "\n============================================================\n";

    std::mt19937 rng(seed);
    metaTorus network(dimension, k);

    double fault_rate = 0.4;
    network.setFaultyLinks(fault_rate);
    network.calculateDirectedRoutingProbabilities();

    int completed = 0;
    while (completed < num_traces) {
        std::uniform_int_distribution<int> node_dist(0, (int)network.V - 1);
        int start = node_dist(rng);
        int end = node_dist(rng);

        if (start == end) continue;

        if (network.findShortestPathBFS(start, end) != metaTorus::DELIVERY_FAIL) {
            std::cout << "\n--- Trial " << (completed + 1) << " ---" << std::endl;

            // 1. Brute Force
            int hops_brute = network.findPathBrute(start, end, true);
            std::cout << "Brute: " << (hops_brute == -1 ? "FAIL" : "OK") << std::endl;

            // 2. Strat 2
            int hops_strat2 = network.findPathDirectedStrategic2(start, end, true);
            std::cout << "Strat2: " << (hops_strat2 == -1 ? "FAIL" : "OK") << std::endl;

            // 3. Strat 3
            int hops_strat3 = network.findPathDirectedStrategic3(start, end, true);
            std::cout << "Strat3: " << (hops_strat3 == -1 ? "FAIL" : "OK") << std::endl;

            completed++;
        }
    }
}

/*

void run_large_scale_experiment(int dimension, int k, int num_iterations, double fault_ratio, unsigned int seed) {
    ExperimentStats stats;
    int successful_trials = 0;
    long long total_attempts = 0;

    std::mt19937 rng(seed);

    // Loop until we get 'num_iterations' of VALID paths (path exists in BFS)
    while (successful_trials < num_iterations) {
        // 1. New Topology per trial (Statistically Robust)
        metaTorus network(dimension, k);
        network.setFaultyLinks(fault_ratio);

        total_attempts++;
        std::uniform_int_distribution<int> node_dist(0, (int)network.V - 1);
        int start_node = node_dist(rng);
        int end_node = node_dist(rng);

//        if (start_node == end_node) continue;

        // 2. Baseline Check (BFS)
        int bfs_length = network.findShortestPathBFS(start_node, end_node);

        if (bfs_length != metaTorus::DELIVERY_FAIL) {
            // Only calculate DP probabilities if the path actually exists
            network.calculateDirectedRoutingProbabilities();

            successful_trials++;
            stats.total_bfs_hops += bfs_length;

            // A. Brute Force
            int brute_ps_length = network.findPathBrute(start_node, end_node, false);
            if (brute_ps_length != metaTorus::DELIVERY_FAIL) {
                stats.brute_ps_successes++;
                stats.total_brute_ps_hops += brute_ps_length;
            }

            // --- NEW: STRATEGIC STANDARD ---
            int strat_std_len = network.findPathStrategic(start_node, end_node, false);
            if (strat_std_len != metaTorus::DELIVERY_FAIL) {
                stats.strategic_std_successes++;
                stats.total_strategic_std_hops += strat_std_len;
            }

            // B. Brute Force Probabilistic
            int brute_so_length = network.findPathBruteProbabilistic(start_node, end_node, false);
            if (brute_so_length != metaTorus::DELIVERY_FAIL) {
                stats.brute_so_successes++;
                stats.total_brute_so_hops += brute_so_length;
            }

            // C. Directed Strategic 2
            int directed_s2_length = network.findPathDirectedStrategic2(start_node, end_node, false);
            if (directed_s2_length != metaTorus::DELIVERY_FAIL) {
                stats.directed_s2_successes++;
                stats.total_directed_s2_hops += directed_s2_length;
            }

            // D. Directed Strategic 3
            int directed_s3_length = network.findPathDirectedStrategic3(start_node, end_node, false);
            if (directed_s3_length != metaTorus::DELIVERY_FAIL) {
                stats.directed_s3_successes++;
                stats.total_directed_s3_hops += directed_s3_length;
            }
        }

        // Safety break for extremely fragmented networks
        if (total_attempts > num_iterations * 500 && successful_trials < 1) {
            break; // Abort this ratio
        }
    }

    // --- Statistics ---
    auto get_sr = [&](int success) { return (successful_trials > 0) ? (100.0 * success / successful_trials) : 0.0; };
    auto get_avg = [&](long long total, int success) { return (success > 0) ? ((double)total / success) : 0.0; };

    // --- CRITICAL: LOCK OUTPUT FOR PARALLEL PRINTING ---
    {
        std::lock_guard<std::mutex> lock(print_mutex);

        std::cout << std::fixed << std::setprecision(2)
                  << std::left  << std::setw(10) << fault_ratio
                  << std::right
                  << std::setw(10) << get_avg(stats.total_bfs_hops, successful_trials)
                  // Brute
                  << std::setw(12) << get_sr(stats.brute_ps_successes) << std::setw(10) << get_avg(stats.total_brute_ps_hops, stats.brute_ps_successes)
                  // Brute Prob
                  << std::setw(12) << get_sr(stats.brute_so_successes) << std::setw(10) << get_avg(stats.total_brute_so_hops, stats.brute_so_successes)
                  // NEW: Strategic Standard
                  << std::setw(12) << get_sr(stats.strategic_std_successes) << std::setw(10) << get_avg(stats.total_strategic_std_hops, stats.strategic_std_successes)
                  // Strat 2
                  << std::setw(12) << get_sr(stats.directed_s2_successes) << std::setw(10) << get_avg(stats.total_directed_s2_hops, stats.directed_s2_successes)
                  // Strat 3
                  << std::setw(12) << get_sr(stats.directed_s3_successes) << std::setw(10) << get_avg(stats.total_directed_s3_hops, stats.directed_s3_successes)
                  << std::endl;
    }
}

void runExperiments() {
    const int NUM_ITERATIONS = 100;
    // You can adjust configs here
    std::vector<std::pair<int, int>> configs = {{5, 3}, {5, 4}, {6, 4}, {4, 12}, {5, 12}};

    std::random_device rd;
    unsigned int master_seed = rd();

    std::cout << "MASTER_SEED: " << master_seed << std::endl;

    for (const auto& config : configs) {
        int n = config.first;
        int k = config.second;

        std::cout << "\n--- CONFIG: n=" << n-1 << ", k=" << k << " ---\n";
        // CSV Header
        // --- PRINT HEADER WITH MATCHING WIDTHS ---
        std::cout << std::left  << std::setw(10) << "F_Ratio"
                  << std::right
                  << std::setw(10) << "BFS_L"
                  << std::setw(12) << "Brute_SR" << std::setw(10) << "Brute_L"
                  << std::setw(12) << "B_Prob_SR" << std::setw(10) << "B_Prob_L"
                  << std::setw(12) << "StratStd_SR" << std::setw(10) << "StratStd_L" // <--- NEW HEADER
                  << std::setw(12) << "Strat2_SR" << std::setw(10) << "Strat2_L"
                  << std::setw(12) << "Strat3_SR" << std::setw(10) << "Strat3_L"
                  << std::endl;
        std::cout << std::string(132, '-') << std::endl; // Increased separator length to fit new column

        std::vector<std::future<void>> futures;

        // Launch Parallel Threads
        for (double ratio = 0.0; ratio <= 0.5; ratio += 0.05) {
            unsigned int run_seed = master_seed + (n * 1000) + (k * 100) + static_cast<unsigned int>(ratio * 100);

            // Async Launch
            futures.push_back(std::async(std::launch::async,
                                         run_large_scale_experiment,
                                         n, k, NUM_ITERATIONS, ratio, run_seed));
        }

        // Wait for all threads to finish this config
        for (auto& f : futures) {
            f.get();
        }
    }
    std::cout << "\n--- All experiments complete. ---" << std::endl;
}

*/


// 2. CHANGE RETURN TYPE: void -> std::string
std::string run_large_scale_experiment_task(int dimension, int k, int num_iterations, double fault_ratio, unsigned int seed) {
    ExperimentStats stats;
    int successful_trials = 0;
    long long total_attempts = 0;

    std::mt19937 rng(seed);

    while (successful_trials < num_iterations) {
        metaTorus network(dimension, k);
        network.setFaultyLinks(fault_ratio);

        total_attempts++;
        std::uniform_int_distribution<int> node_dist(0, (int)network.V - 1);
        int start_node = node_dist(rng);
        int end_node = node_dist(rng);

        if (start_node == end_node) continue;

        int bfs_length = network.findShortestPathBFS(start_node, end_node);

        if (bfs_length != metaTorus::DELIVERY_FAIL) {
            network.calculateDirectedRoutingProbabilities();
            successful_trials++;
            stats.total_bfs_hops += bfs_length;

            // A. Brute Force
            int brute_ps_length = network.findPathBrute(start_node, end_node, false);
            if (brute_ps_length != metaTorus::DELIVERY_FAIL) {
                stats.brute_ps_successes++;
                stats.total_brute_ps_hops += brute_ps_length;
            }

            // B. Brute Probabilistic
            int brute_so_length = network.findPathBruteProbabilistic(start_node, end_node, false);
            if (brute_so_length != metaTorus::DELIVERY_FAIL) {
                stats.brute_so_successes++;
                stats.total_brute_so_hops += brute_so_length;
            }

            // C. Strategic Standard (The new column you wanted)
            int strat_std_length = network.findPathStrategic(start_node, end_node, false);
            if (strat_std_length != metaTorus::DELIVERY_FAIL) {
                stats.strategic_std_successes++;
                stats.total_strategic_std_hops += strat_std_length;
            }

            // D. Strat 2
            int directed_s2_length = network.findPathDirectedStrategic2(start_node, end_node, false);
            if (directed_s2_length != metaTorus::DELIVERY_FAIL) {
                stats.directed_s2_successes++;
                stats.total_directed_s2_hops += directed_s2_length;
            }

            // E. Strat 3
            int directed_s3_length = network.findPathDirectedStrategic3(start_node, end_node, false);
            if (directed_s3_length != metaTorus::DELIVERY_FAIL) {
                stats.directed_s3_successes++;
                stats.total_directed_s3_hops += directed_s3_length;
            }
        }

        if (total_attempts > num_iterations * 500 && successful_trials < 1) break;
    }

    auto get_sr = [&](int success) { return (successful_trials > 0) ? (100.0 * success / successful_trials) : 0.0; };
    auto get_avg = [&](long long total, int success) { return (success > 0) ? ((double)total / success) : 0.0; };

    // 3. FORMAT INTO STRINGSTREAM instead of printing
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2)
       << std::left  << std::setw(10) << fault_ratio
       << std::right
       << std::setw(10) << get_avg(stats.total_bfs_hops, successful_trials)
       << std::setw(12) << get_sr(stats.brute_ps_successes) << std::setw(10) << get_avg(stats.total_brute_ps_hops, stats.brute_ps_successes)
       << std::setw(12) << get_sr(stats.brute_so_successes) << std::setw(10) << get_avg(stats.total_brute_so_hops, stats.brute_so_successes)
       << std::setw(12) << get_sr(stats.strategic_std_successes) << std::setw(10) << get_avg(stats.total_strategic_std_hops, stats.strategic_std_successes)
       << std::setw(12) << get_sr(stats.directed_s2_successes) << std::setw(10) << get_avg(stats.total_directed_s2_hops, stats.directed_s2_successes)
       << std::setw(12) << get_sr(stats.directed_s3_successes) << std::setw(10) << get_avg(stats.total_directed_s3_hops, stats.directed_s3_successes)
       << std::endl;

    return ss.str(); // Return the formatted string
}

void runExperiments() {
    const int NUM_ITERATIONS = 10000;
    std::vector<std::pair<int, int>> configs = {{6, 5}, {4, 12}, {5, 12}};

    std::random_device rd;
    unsigned int master_seed = rd();

    std::cout << "MASTER_SEED: " << master_seed << std::endl;

    for (const auto& config : configs) {
        int n = config.first;
        int k = config.second;

        std::cout << "\n--- CONFIG: n=" << n-1 << ", k=" << k << " ---\n";
        std::cout << std::left  << std::setw(10) << "F_Ratio"
                  << std::right
                  << std::setw(10) << "BFS_L"
                  << std::setw(12) << "Brute_SR" << std::setw(10) << "Brute_L"
                  << std::setw(12) << "B_Prob_SR" << std::setw(10) << "B_Prob_L"
                  << std::setw(12) << "StrStd_SR" << std::setw(10) << "StrStd_L" // Added Header
                  << std::setw(12) << "Strat2_SR" << std::setw(10) << "Strat2_L"
                  << std::setw(12) << "Strat3_SR" << std::setw(10) << "Strat3_L"
                  << std::endl;
        std::cout << std::string(132, '-') << std::endl;

        // 4. STORE FUTURES IN A VECTOR
        std::vector<std::future<std::string>> futures;

        // Launch them all (Still Parallel!)
        for (double ratio = 0.0; ratio <= 0.5; ratio += 0.05) {
            unsigned int run_seed = master_seed + (n * 1000) + (k * 100) + static_cast<unsigned int>(ratio * 100);

            // Note: We use 'run_large_scale_experiment_task' now
            futures.push_back(std::async(std::launch::async,
                                         run_large_scale_experiment_task,
                                         n, k, NUM_ITERATIONS, ratio, run_seed));
        }

        // 5. PRINT THEM IN ORDER
        // The vector 'futures' has them pushed in order (0.0, 0.05, 0.10...).
        // calling f.get() blocks until THAT specific one is ready.
        for (auto& f : futures) {
            std::string result_row = f.get();
            std::cout << result_row; // This guarantees order
        }
    }
    std::cout << "\n--- All experiments complete. ---" << std::endl;
}
