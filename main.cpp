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
//        torus.setTestedFaultyLinks();
//        torus.printFaultyLinks();
        int seed=52;
        int d_route, d_directed_route, d_brute, d_bfs;
        int from, to;
        torus.setRandomFaultyLinks((double)0.5, &seed);
//        int from = mt() % torus.V;
//        int to = mt() % torus.V;

        // Randomly generate two valid nodes
        do {
            from = get_random(seed++, torus.V);
        } while (!torus.nodes[from].validity); // Ensure 'from' is valid

        do {
            to = get_random(seed++, torus.V);
        } while (!torus.nodes[to].validity); // Ensure 'to' is valid

        cout << "from,to=" << from << "," << to << endl;

        torus.calcRoutingProbabilities();
//        torus.printProbabilities();
        torus.calcDirectedRoutingProbabilities();
//        torus.printDirectedProbabilities();

        std::cout << "Routing from node " << from << " to node " << to << std::endl;
        d_route = torus.route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
        d_directed_route = torus.directed_route(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
        d_brute = torus.brute(nullptr, &(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>), 0);
        d_bfs = torus.bfs(&(torus.nodes[from]), &(torus.nodes[to]), *(new std::unordered_map<int, bool>));
//
        cout << d_route << ", "
        << d_directed_route << " , "
        << d_brute << " , "
        << d_bfs << " , " << endl;

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

int main(int argc, char* argv[]){
    using namespace std;
//    test_torus_constructor(4, 5); // New test case
//    test_torus_constructor(3, 5); // Working case
//    test_torus_constructor(2, 5); // New test case
    test_torus_constructor(3, 3); // Edge case
//    test_torus_constructor(5, 5); // New test case
//    test_torus_constructor(4, 8); // New test case
//    test_torus_constructor(3, 16); // New test case
//    test_torus_constructor(2, 64); // New test case
//    test_torus_constructor(5, 5); // New test case
    return 0;
}
