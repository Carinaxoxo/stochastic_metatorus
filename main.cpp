#include <cmath>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <random>
#include <stdlib.h>
#include <cstdlib>
#include <chrono>

#include "random.h"
#include "metatorus.h"

#define N 10000
#define DELIVERY_FAIL (-1)

//params = [(n, k, round(0.01 * p, 2)) for n, k, p
//        in itertools.product([2], [64], range(0, 55, 5))]


void test_torus_constructor(int n, int k) {
    random_device rd;
    mt19937 mt(rd());
    try {
        metaTorus torus(n, k);
        std::cout << "Torus initialized." << std::endl;
//        torus.display();
//        torus.setTestedFaultyLinks();

        // Clear any faulty links
//        torus.clearFaultyLinks();

        if (!torus.isFullyConnected()) {
            std::cerr << "Error: The torus is not fully connected." << std::endl;
            return;
        } else{
            std::cout << "The torus is fully connected." << std::endl;
        }

//        int seed=52;
//        int d_route, d_directed_route, d_brute, d_bfs;
//        torus.setRandomFaultyLinks((double)0, &seed);
////        int from, to;
//        int from = mt() % torus.V;
//        int to = mt() % torus.V;
        // Ensure no faulty links
//        torus.printFaultyLinks();

        // Randomly generate two valid nodes
//        do {
//            from = get_random(seed++, torus.V);
//        } while (!torus.nodes[from].validity); // Ensure 'from' is valid
//
//        do {
//            to = get_random(seed++, torus.V);
//        } while (!torus.nodes[to].validity); // Ensure 'to' is valid

//        torus.calcRoutingProbabilities();
//        torus.calcDirectedRoutingProbabilities();

//        cout << "Routing from node " << from << " to node " << to << endl;
//        std::unordered_map<int, bool> visited;
//        d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);
//        std::cout << "BFS distance: " << d_bfs << std::endl;

//        int reachable = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);
//        std::cout << "Reachable nodes from 0: " << visited.size() << " / " << torus.V << std::endl;

//
//        torus.calcRoutingProbabilities();
//        torus.printProbabilities();
//        torus.testDirectedRoutingProbabilities();
//        torus.calcDirectedRoutingProbabilities();
//        torus.printDirectedProbabilities();

//        std::cout << "Routing from node " << from << " to node " << to << std::endl;
//        d_route = torus.route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_directed_route = torus.directed_route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_brute = torus.brute(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>));
//
//        cout <<
//        d_route << ", " <<
//        d_directed_route << " , " <<
//        d_brute << " , " <<
//        d_bfs << " , " << endl;

//        std::cout << "Path: ";
//        for (Node* node : path) {
//            std::cout << node->index << " ";
//        }
//        torus.display();
//        torus.printFaultyLinks();
//        torus.testDirectedRoutingProbabilities();

        std::cout << "Torus initialized successfully with n=" << n << ", k=" << k << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "Error initializing Torus with n=" << n << ", k=" << k << ": " << e.what() << std::endl;
    }
}

void test_strategy(int n, int k, int num_tests) {
    random_device rd;
    mt19937 mt(rd());
    int nonEmptyCount = 0;
    for (int i = 0; i < num_tests; ++i) {
        try {
            metaTorus torus(n, k);
//            std::cout << "Torus initialized." << std::endl;
//        torus.display();
//        torus.setTestedFaultyLinks();

            // Clear any faulty links
//        torus.clearFaultyLinks();

            if (!torus.isFullyConnected()) {
                std::cerr << "Error: The torus is not fully connected." << std::endl;
                return;
            } else{
//                std::cout << "The torus is fully connected." << std::endl;
            }

            int seed=52;
            int d_route, d_directed_route, d_brute, d_bfs;
//            torus.setRandomFaultyLinks((double)0, &seed);
//        int from, to;
            int from = mt() % torus.V;
            int to = mt() % torus.V;
//         Ensure no faulty links
//            torus.printFaultyLinks();

            std::vector<std::vector<std::vector<int>>> result = torus.return_node_which_close_to_goal(&(torus.nodes[from]), &(torus.nodes[to]));
            // Print the input nodes
//            std::cout << "Node c (current): ";
//            for (int val : torus.nodes[from].value) {
//                std::cout << val << " ";
//            }
//            std::cout << std::endl;
//
//            std::cout << "Node t (target): ";
//            for (int val : torus.nodes[to].value) {
//                std::cout << val << " ";
//            }
//            std::cout << std::endl;

            // Print the result
//            if (!result.empty()){
//                std::cout << "Nodes closer to the target in 2d:" << std::endl;
//                for (const auto& node : result) {
//                    for (int val : node) {
//                        std::cout << val << " ";
//                    }
//                    std::cout << std::endl;
//                }
//            } else {
//                std::cout << "No valid neighbor found in 2d." << std::endl;
//            }

            if (result.empty()){
                std::cout << "No valid neighbor found." << std::endl;
                std::cout << "Node c (current): ";
                for (int val : torus.nodes[from].value) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;

                std::cout << "Node t (target): ";
                for (int val : torus.nodes[to].value) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            } else {
                nonEmptyCount++;
            }

            // Print the result
//            if (!result.empty()){
//                nonEmptyCount++;
//                std::cout << "Nodes closer to the target in 2d:" << std::endl;
//                for (const auto& node : result) {
//                    for (int val : node) {
//                        std::cout << val << " ";
//                    }
//                    std::cout << std::endl;
//                }
//            } else {
//                std::cout << "No valid neighbor found in 2d." << std::endl;
//                std::cout << "Node c (current): ";
//                for (int val : torus.nodes[from].value) {
//                    std::cout << val << " ";
//                }
//                std::cout << std::endl;
//
//                std::cout << "Node t (target): ";
//                for (int val : torus.nodes[to].value) {
//                    std::cout << val << " ";
//                }
//                std::cout << std::endl;
//            }


            // Randomly generate two valid nodes
//        do {
//            from = get_random(seed++, torus.V);
//        } while (!torus.nodes[from].validity); // Ensure 'from' is valid
//
//        do {
//            to = get_random(seed++, torus.V);
//        } while (!torus.nodes[to].validity); // Ensure 'to' is valid

//        torus.calcRoutingProbabilities();
//        torus.calcDirectedRoutingProbabilities();

//        cout << "Routing from node " << from << " to node " << to << endl;
//        std::unordered_map<int, bool> visited;
//        d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);
//        std::cout << "BFS distance: " << d_bfs << std::endl;

//        int reachable = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);
//        std::cout << "Reachable nodes from 0: " << visited.size() << " / " << torus.V << std::endl;

//
//        torus.calcRoutingProbabilities();
//        torus.printProbabilities();
//        torus.testDirectedRoutingProbabilities();
//        torus.calcDirectedRoutingProbabilities();
//        torus.printDirectedProbabilities();

//        std::cout << "Routing from node " << from << " to node " << to << std::endl;
//        d_route = torus.route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_directed_route = torus.directed_route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_brute = torus.brute(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
//        d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>));
//
//        cout <<
//        d_route << ", " <<
//        d_directed_route << " , " <<
//        d_brute << " , " <<
//        d_bfs << " , " << endl;

//        std::cout << "Path: ";
//        for (Node* node : path) {
//            std::cout << node->index << " ";
//        }
//        torus.display();
//        torus.printFaultyLinks();
//        torus.testDirectedRoutingProbabilities();

//            std::cout << "Torus initialized successfully with n = " << n << ", k = " << k << std::endl;
        } catch (const std::exception &e) {
            std::cerr << "Error initializing Torus with n=" << n << ", k=" << k << ": " << e.what() << std::endl;
        }
    }
    std::cout << "Non-empty wanted ratio: " << nonEmptyCount << "/" << num_tests<< std::endl;
}

void test_bfs_multiple_times(int n, int k, int num_tests){
    std::cout << "Running BFS test " << num_tests << " times with n=" << n << ", k=" << k << std::endl;
    random_device rd;
    mt19937 mt(rd());
    int seed = 52;
    int successful_tests = 0;
    for (int i = 0; i < num_tests; ++i) {
        try {
            metaTorus torus(n, k);

            if (!torus.isFullyConnected()) {
                std::cerr << "Error: The torus is not fully connected." << std::endl;
                continue;
            }

            int from = mt() % torus.V;
            int to = mt() % torus.V;

            int iterations = 0;
            std::unordered_map<int, bool> visited;
            int d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);

            if (d_bfs != DELIVERY_FAIL) {
                successful_tests++;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error during BFS test: " << e.what() << std::endl;
        }
    }

    std::cout << "BFS test completed. Successful tests: " << successful_tests << " / " << num_tests << std::endl;
}

void test_all_algorithms_multiple_times(int n, int k, int num_tests){
    std::cout << "Running all algorithms test " << num_tests << " times with n=" << n << ", k=" << k << std::endl;
    random_device rd;
    mt19937 mt(rd());
    int seed = 52;
    int successful_tests = 0;

    int n_brute_success = 0, n_route_success = 0, n_directed_route_success = 0;
    double mean_d_route = 0, mean_d_directed_route = 0, mean_d_brute = 0, mean_d_bfs = 0;

    for (int i = 0; i < num_tests; ++i) {
        try{
            metaTorus torus(n, k);

            if (!torus.isFullyConnected()) {
                std::cerr << "Error: The torus is not fully connected." << std::endl;
                continue;
            } else {
                cout << "it is fully connected" << endl;
            }

            torus.setRandomFaultyLinks((double)0, &seed);

            int from = mt() % torus.V;
            int to = mt() % torus.V;
            std::unordered_map<int, bool> visited;
            int d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);

            if (d_bfs != DELIVERY_FAIL) {
                successful_tests++;
            }

            torus.calcRoutingProbabilities();
            torus.calcDirectedRoutingProbabilities();

//            int d_brute = torus.brute(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
            int d_route = torus.route_test(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
            int d_directed_route = torus.directed_route_test(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);

//            if (d_brute != DELIVERY_FAIL) {
//                n_brute_success++;
//                mean_d_brute += d_brute;
//            }

            if (d_route != DELIVERY_FAIL) {
                n_route_success++;
                mean_d_route += d_route;
            }
//
            if (d_directed_route != DELIVERY_FAIL) {
                n_directed_route_success++;
                mean_d_directed_route += d_directed_route;
            }

            mean_d_bfs += d_bfs;
        } catch (const std::exception &e) {
            std::cerr << "Error during test: " << e.what() << std::endl;
        }
    }

    double p_route_success = n_route_success / (double) num_tests;
    double p_directed_route_success = n_directed_route_success / (double) num_tests;
//    double p_brute_success = n_brute_success / (double) num_tests;

    std::cout << "All algorithms test completed." << std::endl;
    std::cout << "Successful tests: " << successful_tests << " / " << num_tests << std::endl;
    std::cout << "Route success rate: " << p_route_success << std::endl;
    std::cout << "Directed route success rate: " << p_directed_route_success << std::endl;
//    std::cout << "Brute force success rate: " << p_brute_success << std::endl;
    std::cout << "Mean BFS distance: " << mean_d_bfs / num_tests << std::endl;
    std::cout << "Mean route distance: " << mean_d_route / n_route_success << std::endl;
//    std::cout << "Mean directed route distance: " << mean_d_directed_route / n_directed_route_success << std::endl;
//    std::cout << "Mean brute force distance: " << mean_d_brute / n_brute_success << std::endl;
}

void test_brute(int n, int k, int num_tests){
    std::cout << "Running brute test " << num_tests << " times with n=" << n << ", k=" << k << std::endl;
    random_device rd;
    mt19937 mt(rd());
    int seed = 52;
    int successful_tests = 0;
    double mean_d_brute = 0;

    for (int i = 0; i < num_tests; ++i) {
        try {
            metaTorus torus(n, k);

            if (!torus.isFullyConnected()) {
                std::cerr << "Error: The torus is not fully connected." << std::endl;
                continue;
            }

            int from = mt() % torus.V;
            int to = mt() % torus.V;

            int iterations = 0;
            std::unordered_map<int, bool> visited;
            int d_brute = torus.brute(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);

            if (d_brute != DELIVERY_FAIL) {
                successful_tests++;
                mean_d_brute += d_brute;
            }

        } catch (const std::exception &e) {
            std::cerr << "Error during brute test: " << e.what() << std::endl;
        }
    }
//    double p_brute_success = successful_tests / (double) num_tests;

    std::cout << "brute test completed. Successful tests: " << successful_tests << " / " << num_tests << std::endl;
    std::cout << "Mean brute force distance: " << mean_d_brute / successful_tests << std::endl;
}

void test_directed(int n, int k, int num_tests){
    std::cout << "Running directed test " << num_tests << " times with n=" << n << ", k=" << k << std::endl;
    random_device rd;
    mt19937 mt(rd());
    int seed = 52;
    int successful_tests = 0;

    int n_directed_route_success = 0;
    double mean_d_directed_route = 0, mean_d_bfs = 0;

    for (int i = 0; i < num_tests; ++i) {
        try{
            metaTorus torus(n, k);

            if (!torus.isFullyConnected()) {
                std::cerr << "Error: The torus is not fully connected." << std::endl;
                continue;
            } else {
                cout << "it is fully connected" << endl;
            }

            int from = mt() % torus.V;
            int to = mt() % torus.V;
            std::unordered_map<int, bool> visited;
            int d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), visited);

            if (d_bfs != DELIVERY_FAIL) {
                successful_tests++;
            }

            torus.testDirectedRoutingProbabilities();
//            torus.printDirectedProbabilities();
            int d_directed_route = torus.directed_route_test(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);

//
            if (d_directed_route != DELIVERY_FAIL) {
                n_directed_route_success++;
                mean_d_directed_route += d_directed_route;
            }

            mean_d_bfs += d_bfs;
        } catch (const std::exception &e) {
            std::cerr << "Error during test: " << e.what() << std::endl;
        }
    }

    double p_directed_route_success = n_directed_route_success / (double) num_tests;

    std::cout << "All algorithms test completed." << std::endl;
    std::cout << "Successful tests: " << successful_tests << " / " << num_tests << std::endl;
    std::cout << "Directed route success rate: " << p_directed_route_success << std::endl;
    std::cout << "Mean BFS distance: " << mean_d_bfs / num_tests << std::endl;
    std::cout << "Mean directed route distance: " << mean_d_directed_route / n_directed_route_success << std::endl;

}


int main(int argc, char* argv[]){
    using namespace std;
    test_strategy(4, 3, 1000);
//     test_strategy(5, 4, 1000);
//    test_all_algorithms_multiple_times(4, 3, 1000);
//    test_all_algorithms_multiple_times(5, 4, 1000);
//    test_bfs_multiple_times(4, 3, 10000); // Run BFS test 10,000 times
//    test_bfs_multiple_times(5, 4, 10000); // Run BFS test 10,000 times
//    test_torus_constructor(4, 3); // Edge case
//    test_torus_constructor(5, 4); // New test case
//     test_brute(4, 3, 1000);
//     test_brute(5, 4, 1000);
//    test_directed(4, 3, 1000);
//     test_directed(5, 4, 1000);
//    test_proposal();
    return 0;
}




/*
int main(int argc, char* argv[]) {
    using namespace std;

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    double p_faulty = atof(argv[3]);
    random_device rd;
    mt19937 mt(rd());
//
    double mean_d_route = 0, mean_d_directed_route = 0, mean_d_brute = 0, mean_d_bfs = 0;
    int n_brute_success = 0, n_route_success = 0, n_directed_route_success = 0;
    int n_loop = N;
    vector<int> brute_p, route_p, bfs_p;
//
    bool *route_success = new bool[n_loop];
    bool *directed_route_success = new bool[n_loop];
    bool *brute_success = new bool[n_loop];

    double *d_routes = new double[n_loop];
    double *d_directed_routes = new double[n_loop];
    double *d_brutes = new double[n_loop];
    double *d_bfss = new double[n_loop];
//
    //uniform_int_distribution<int> dice(0, static_cast<int>(pow(k, n) - 1));

    assert(0.0 <= p_faulty and p_faulty < 1.0);

//     パラレルに実行される 同じデータにアクセスしないように注意
    int seed = atoi(argv[5]);

    for (int i = atoi(argv[4]); i < n_loop; i++) {

        bool has_non_faulty_route = false;
        while (true) {
            cout << "loop: " << i << endl;
//             step 1
            metaTorus *t = new metaTorus(n, k);
            t->setRandomFaultyLinks((double) p_faulty, &seed);

//             step 2
            int from = get_random(seed++, t->V);//dice(mt);
            int to = get_random(seed++, t->V);//dice(mt);
//            cout << "from,to=" << from << "," << to << endl;
            d_bfss[i] = t->bfs(&(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>));
//            get the bfs route
            has_non_faulty_route = (d_bfss[i] != DELIVERY_FAIL);
            if (not has_non_faulty_route) {
                delete t;
                continue;

            }
            //it has to be N times successfully run
//            else {
//                cout <<"has non faulty route, yeah" << endl;
//                cout << "d_bfss[i] = " << d_bfss[i] << endl;
//            }

//             step 3
            d_brutes[i] = t->brute(nullptr, &(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>),
                                   0);
            brute_success[i] = (d_brutes[i] != DELIVERY_FAIL);
//            cout << "brute_success[i] = " << brute_success[i] << endl;
//            cout << "d_brutes[i] = " << d_brutes[i] << endl;
//            get the brute path

            t->calcRoutingProbabilities();
            t->calcDirectedRoutingProbabilities();
//            t->printProbabilities();
//            t->printDirectedProbabilities();
//            return 0;
            d_routes[i] = t->route(nullptr, &(t->nodes[from]), &(t->nodes[to]), *(new std::unordered_map<int, bool>),
                                   0);
            route_success[i] = (d_routes[i] != DELIVERY_FAIL);
//            cout << "route_success[i] = " << route_success[i] << endl;
//            cout << "d_routes[i] = " << d_routes[i] << endl;
            d_directed_routes[i] = t->directed_route(nullptr, &(t->nodes[from]), &(t->nodes[to]),
                                                     *(new std::unordered_map<int, bool>), 0);
            directed_route_success[i] = (d_directed_routes[i] != DELIVERY_FAIL);
//            cout << "directed_route_success[i] = " << directed_route_success[i] << endl;
//            cout << "d_directed_routes[i] = " << d_directed_routes[i] << endl;

            delete t;
            break;
        }

//        cout << i << "," << seed << "," << route_success[i] << "," << brute_success[i];
//        if(route_success[i])
//            cout << "," << d_routes[i];
//        else
//            cout << "," << -1;
//
//        if(directed_route_success[i])
//            cout << "," << d_directed_routes[i];
//        else
//            cout << "," << -1;
//
//        if(brute_success[i])
//            cout << "," << d_brutes[i];
//        else
//            cout << "," << -1;
//        cout << endl;
    }

//     データを集計する
    for (int i = 0; i < n_loop; i++) {
        if (route_success[i]) {
            n_route_success++;
            mean_d_route += d_routes[i];
        }
//        cout << mean_d_route << endl;

        if (directed_route_success[i]) {
            n_directed_route_success++;
            mean_d_directed_route += d_directed_routes[i];
        }
//        cout << mean_d_directed_route << endl;


        if (brute_success[i]) {
            n_brute_success++;
            mean_d_brute += d_brutes[i];
        }

        assert(d_bfss[i] != DELIVERY_FAIL);
        mean_d_bfs += d_bfss[i];
    }
//
    double p_route_success = n_route_success / (double) n_loop;
    double p_directed_route_success = n_directed_route_success / (double) n_loop;
    double p_brute_success = n_brute_success / (double) n_loop;

    cout << p_faulty << ", "
         << n << ", "
         << k << ", "
         << p_route_success << ", "
         << p_directed_route_success << ", "
         << p_brute_success << ", "
         << mean_d_route / n_route_success << ", "
         << mean_d_directed_route / n_directed_route_success << ", "
         << mean_d_brute / n_brute_success << ", "
         << mean_d_bfs / n_loop
         << endl;

}

*/


