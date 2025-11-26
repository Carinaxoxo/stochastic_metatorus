//
// Created by Carina Z on 2025/11/20.
//

#ifndef RESTART_EXPERIMENTS_H
#define RESTART_EXPERIMENTS_H

#include <iostream>
#include <vector>
#include <iomanip>
#include <random>
#include "metatorus.h"

struct ExperimentStats {
    int brute_ps_successes = 0;
    int brute_so_successes = 0;
    int strategic_std_successes = 0; // <--- NEW
    int directed_s2_successes = 0;
    int directed_s3_successes = 0;

    long long total_bfs_hops = 0;
    long long total_brute_ps_hops = 0;
    long long total_brute_so_hops = 0;
    long long total_strategic_std_hops = 0; // <--- NEW
    long long total_directed_s2_hops = 0;
    long long total_directed_s3_hops = 0;
};

void run_visual_trace_experiment(int dimension, int k, int num_traces);
void runExperiments_track();

void run_large_scale_experiment(int dimension, int k, int num_iterations, double fault_ratio, unsigned int seed);
void runExperiments();

void run_strategy_comparison(int dimension, int k, int num_traces, unsigned int seed);

void test_recursive_strategic_vs_bfs(int n, int k, int trials, double fault_rate);

#endif //RESTART_EXPERIMENTS_H
