#include <stdexcept>
#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_set>
#include <cassert>
#include <unordered_map>
#include <queue>
#include <random>
#include <cstring>
#include <cstdlib>
#include <utility>
#include <iomanip>
#include <set>

#include "bbst.h"
#include "random.h"

//#include "boost/numeric/ublas/matrix_sparse.hpp"
//#include "boost/numeric/ublas/io.hpp"

#include "metatorus.h"

#define DELIVERY_FAIL -1;

metaTorus::metaTorus(int _n, int _k) {
    assert(_k > 2);
    assert(_n >= 0);
    //to make sure it is different from hypercube

    // calculate the parameters for this torus
    n = _n;
    k = _k;
    V = static_cast<int>(pow(k, n));
    E = n*V;
    diameter = static_cast<int>(floor(k/2)*n);

    // 隣接行列と故障リンク行列をV*Vのスパース行列として初期化
    //adjacency = mapped_matrix<int>(V, V);
    //F = mapped_matrix<int>(V, V);
    F.resize(V, vector<int>(V, -1));
//    for (auto& row : F) {
//        row.resize(V);
//    }

    // メモリを確保する
    allocate();
    // nodes->valueをセットする
    setValues();
    // トーラスの接続行列を計算
    createConnection();
}

void metaTorus::display() {
    std::cout << "Torus: n=" << n << ", k=" << k << std::endl;
    std::cout << "Number of vertices (V): " << V << std::endl;
    std::cout << "Number of edges (E): " << E << std::endl;
    std::cout << "Diameter: " << diameter << std::endl;

    std::cout << "Node coordinates and neighbors:" << std::endl;
    for (int i = 0; i < V; ++i) {
        std::cout << "Node " << i << " (";
        for (int j = 0; j < n; ++j) {
            std::cout << nodes[i].value[j];
            if (j < n - 1) std::cout << ",";
        }
        std::cout << "): ";

        int validity = (nodes[i].validity) ? 1 : 0;
        cout << "validity: " << validity << ": ";

        // Print neighbors
        for (auto neighbor : nodes[i].neighbors) {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;
    }

//    std::cout << "Adjacency matrix (F):" << std::endl;
//    for (int i = 0; i < V; ++i) {
//        for (int j = 0; j < V; ++j) {
//            std::cout << F[i][j] << " ";
//        }
//        std::cout << std::endl;
//    }
    std::cout << std::endl;
}

void metaTorus::allocate() {
    // O(n k^n)
//    nodes = (Node*)malloc(V*sizeof(Node));
    nodes.resize(V);
//    nodes.reserve(V);
    for (int i = 0; i < V; ++i) {
        nodes[i].index = i;
        nodes[i].value.resize(n, 0);
        nodes[i].validity = false;
        // Initialize P matrix
        nodes[i].P.resize(n + 1, vector<double>(diameter + 1, -1));

        // Initialize 3D array d_P
        nodes[i].d_P.resize(n + 1, vector<vector<double>>(diameter + 1, vector<double>(4, -1)));

        nodes[i].neighbors.resize(4, -1);
    }
}

// インデックスをk進数に変換するようなアルゴリズム
//stay the same
void metaTorus::setValues() {
    // O(n k^n)
    int m, po;
    for(int i=0; i<V; i++) {
        //go through by id
        m = i;
        // for all dimension
        for(int j=n-1; j>0; j--) {
            //from most significant to least
//            po = int(pow(k, j));
            po = static_cast<int>(pow(k, j));
            (nodes)[i].value[n-j-1] = int(m / po);
            m %= po;
        }

        (nodes)[i].value[n-1] = m;

        nodes[i].validity = (0 <= nodes[i].value[0] && nodes[i].value[0] < (n - 1));

//        std::string v = "";
//        for(int j=0; j<n; j++)
//            v += to_string(nodes[i].value[j]) + ",";
//        cout << v << "=" << i << endl;
//        node_value_index[v] = i;
    }
}

void metaTorus::createConnection() {
    // O(n^2 k^n)
    using namespace std;
    for (int i = 0; i < V; ++i) {
        vector<int> current = nodes[i].value;

//        std::cout << "current " << i << " (";
//        for (int j = 0; j < n; ++j) {
//            std::cout << current[j];
//            if (j < n - 1) std::cout << ",";
//        }
//        std::cout << "): ";
        nodes[i].neighbors.clear();

        // Compute neighbors in the 4 directions
//        std::set<int> unique_neighbors;

        // Change the first bit (+1, -1)
        vector<int> temp = current;
        // Neighbor 1: Increment the first coordinate
        temp[0] = (current[0] + 1) % (n - 1);

        int neighbor1 = getId(temp);

        // Neighbor 2: Decrement the first coordinate
        temp = current; // rest
        temp[0] = (current[0] - 1 + k) % (n - 1);
        int neighbor2 = getId(temp);

        // Neighbor 3: Increment the x-th coordinate
        temp = current;
        int x = (current[0] + 1); // Ensure x is a valid index
        temp[x] = (temp[x] + 1) % k;
        int neighbor3 = getId(temp);

        // Neighbor 4: Decrement the x-th coordinate
        temp = current;
        temp[x] = (temp[x] - 1 + k) % k;
        int neighbor4 = getId(temp);

        nodes[i].neighbors.push_back(neighbor1);
        nodes[i].neighbors.push_back(neighbor2);
        nodes[i].neighbors.push_back(neighbor3);
        nodes[i].neighbors.push_back(neighbor4);

        // Ensure the neighbor list has exactly 4 elements
        if (nodes[i].neighbors.size() != 4) {
            std::cerr << "Error: Node " << i << " has incorrect neighbors. Found "
                      << nodes[i].neighbors.size() << " neighbors." << std::endl;
        }

//        // Update the adjacency matrix
        F[i][neighbor1] = 1;
        F[i][neighbor2] = 1;
        F[i][neighbor3] = 1;
        F[i][neighbor4] = 1;
    }
}

void metaTorus::setTestedFaultyLinks() {

// previous research version
//    setFaultyLink(&nodes[1], &nodes[2]);
//    setFaultyLink(&nodes[2], &nodes[3]);
//    setFaultyLink(&nodes[6], &nodes[7]);
//    setFaultyLink(&nodes[11], &nodes[16]);
//    setFaultyLink(&nodes[12], &nodes[17]);
//    setFaultyLink(&nodes[14], &nodes[10]);
//    setFaultyLink(&nodes[18], &nodes[19]);
//    setFaultyLink(&nodes[19], &nodes[15]);
//    setFaultyLink(&nodes[21], &nodes[1]);
//    setFaultyLink(&nodes[21], &nodes[22]);

// current research version
    setFaultyLink(&nodes[0], &nodes[5]);
    setFaultyLink(&nodes[0], &nodes[4]);
    setFaultyLink(&nodes[10], &nodes[11]);
    setFaultyLink(&nodes[20], &nodes[24]);
    setFaultyLink(&nodes[6], &nodes[11]);
    setFaultyLink(&nodes[16], &nodes[21]);
    setFaultyLink(&nodes[7], &nodes[12]);
    setFaultyLink(&nodes[22], &nodes[23]);
    setFaultyLink(&nodes[8], &nodes[13]);
    setFaultyLink(&nodes[24], &nodes[4]);
}

void metaTorus::setRandomFaultyLinks(double p_faulty, int* seed) {
    std::mt19937 gen(*seed);
    std::vector<std::pair<int, int>> all_links;

    // Generate all unique links in the torus
    for (int i = 0; i < V; ++i) {

        if (!nodes[i].validity) continue; // Skip invalid nodes

        for (int j : nodes[i].neighbors) {
            if (i < j) {
                all_links.emplace_back(i, j);
            }
        }
    }

    // Shuffle the list of links
    std::shuffle(all_links.begin(), all_links.end(), gen);

    // Calculate the number of faulty links
    int num_faulty_links = static_cast<int>(all_links.size() * p_faulty);

    // Select the top num_faulty_links from the shuffled list to be faulty
    for (int i = 0; i < num_faulty_links && i < all_links.size(); ++i) {
        int a_id = all_links[i].first;
        int b_id = all_links[i].second;
        setFaultyLink(&nodes[a_id], &nodes[b_id]);
    }
}

void metaTorus::setFaultyLink(Node *a, Node *b) {
    F[a->index][b->index] = F[b->index][a->index] = 0;
    //two directions
}

bool metaTorus::hasFaultyLink(Node *a, Node *b) {
    // 対称行列で上三角しか使ってないので入れ替えて補う
    bool ab = static_cast<bool>(F[a -> index][b -> index] == 0);
    bool ba = static_cast<bool>(F[b -> index][a -> index] == 0);

    //id is not constent

//    cout << "Checking link between " << a->index << " and " << b->index << ": "
//         << "ab=" << ab << ", ba=" << ba << endl;

//    assert(not (ab and ba));
    //both have to be non-faulty

    return (ab or ba);
}

int metaTorus::getId(vector<int> input){
    int result = 0;
    int power = 1;
    for (int i = n-1; i >= 0; i--){
        result += input[i] * power;
        power *= k;
    }

    return result;
}

bool metaTorus::isNeighbor(Node *a, Node *b) {
    // Number of differing dimensions and index of the differing dimension
    int hamming_dist = 0, i_different = 0, i_check;
    int h_a, h_b;
    //values at the header

    // Values at the differing index
    int a_differ, b_differ;
    h_a = a->value[0];
    h_b = b->value[0];
    i_check = h_a;

    // check all dimensions
    for (int i = 0; i < n; i++) {
        if (a->value[i] != b->value[i]) {
            hamming_dist++;
            i_different = i; // Record the index of the differing dimension
        }
    }

    // If more than one dimension differs, they cannot be neighbors
    if (hamming_dist != 1)
        return false;

    if (i_different == 0 && h_a != h_b){
        if (abs(h_a - h_b) == 1)
            return true;
        // Case 2.2: Wrap-around difference (e.g., 0 and k-1)
        if ((h_a == n - 2 && h_b == 0) || (h_a == 0 && h_b == n - 2))
            return true;
    }

    if (i_different == i_check + 1 && h_a == h_b) {
        // Get the differing values in the differing dimension
        a_differ = a->value[i_different];
        b_differ = b->value[i_different];

        if(abs(a_differ - b_differ) == 1)
            return true;

        // 0とk-1のパターン
        if((a_differ == k-1 and b_differ == 0) or (a_differ == 0 and b_differ == k-1))
            return true;
    }

    // Additional metatorus-specific neighbor rules
    // Here you can add logic to handle any unique metatorus connections.
    // For example, if metatorus has long-range or diagonal connections:
    // if (some_condition_based_on_metatorus_topology)
    //     return true;

    // If no conditions are met, they are not neighbors
    return false;
}

int metaTorus::getSpecificNeighbor(Node *a, int neighborIndex) {
    // Check if the neighborIndex is valid
    if (neighborIndex < 0 || neighborIndex >= nodes.size()) {
        std::cerr << "Invalid neighbor index: " << neighborIndex << std::endl;
        return -1; // Return an invalid index if the neighborIndex is out of bounds
    }

    // Iterate through the neighbors of node `a` to find the index of the specific neighbor
    for (int i = 0; i < a->neighbors.size(); i++) {
        if (a->neighbors[i] == neighborIndex) {
            return i; // Return the index if the neighbor is found
        }
    }

    // If no matching neighbor is found, return -1
    return -1;
}

int metaTorus::distance(Node *a, Node *b) {
    int d = 0;
    for(int i=0; i<n; i++)
        if(abs(a->value[i]-b->value[i]) > floor(k/2))
            d += k-abs(a->value[i]-b->value[i]);
        else
            d += abs(a->value[i]-b->value[i]);

    assert(0 <= d and d <= diameter);
    return d;
}

int metaTorus::hammingDistance(Node *a, Node *b) {
    int d = 0;
    for(int i=0; i<n; i++)
        if(a->value[i] != b->value[i])
            d++;

    assert(0 <= d and d <= n);
    return d;
}

bool metaTorus::inPre(Node *neighbor, Node *a, Node *b) {
    return distance(neighbor, b) < distance(a, b);
}


bool metaTorus::inSpr(Node *neighbor, Node *a, Node *b) {
    return not inPre(neighbor, a, b);
}

void metaTorus::printFaultyLinks() {
    int faultyLinkCount = 0;

    for(int i=0; i < V; i++) {
        for (int j=0; j<V; j++)
            if(F[i][j] == 0){
                faultyLinkCount++;
//                cout << "(" << i << "," << j << ")" << endl;
                std::cout << "((";
                for(int l=0; l<n; l++) {
                    std::cout << nodes[i].value[l];
                    if(l < n-1)
                        std::cout << ", ";
                }
                //std::cout << std::endl;
                std::cout << "), (";
                for(int l=0; l<n; l++){
                    std::cout << nodes[j].value[l];
                    if(l < n-1)
                        std::cout << ", ";
                }
                std::cout << "))" << std::endl;
            };
    }
    cout << "Total number of faulty links: " << faultyLinkCount << endl;
}

void metaTorus::calcRoutingProbabilities() {
    Node *a;
    Node *neighbor;
    double sum, p11;

//    cout << "Starting calculation of routing probabilities..." << endl;
    /*
     * まず全てのノードについてP(a)_{1,1}を計算
     * First, calculate P(a)_{1,1} for all nodes
     */
    for (int i = 0; i < V; i++) {
        a = &nodes[i];
        p11 = 0.0;

        if (!a->validity) continue;

        for (int j = 0; j < 4; j++) {
            if (j >= sizeof(a->neighbors)) {
                cerr << "Error: Neighbor index " << j << " out of bounds for node " << i << " in initial loop" << endl;
                continue;
            }

            neighbor = &nodes[a->neighbors[j]];

            if (!neighbor) {
                cerr << "Error: Neighbor node is null for index " << j << " in initial loop" << endl;
                continue;
            }

            if (not hasFaultyLink(a, neighbor)) {
                p11 += 1.0;
//                cout << "not faulty link" << endl;
            }
        }
        p11 /= 4;
//        setProbability(a, 0, 0, 1);
        setProbability(a, 1, 1, p11);

//        cout << "p_1,1 = " << p11 << endl;
    }

    /*
     * P(a)_{h,d}を計算していく
     * P(a)_{h,d}の計算には, P(a)_{h,d-1}やP(a)_{h-1,d-1}が必要になる.
     * P(a)_{1,1}だけで計算できるP(a)_{1,2}から順番に計算していく(動的計画法)
     *
     * Calculation of P(a)_{h,d}
     * P(a)_{h,d} needs P(a)_{h,d-1} and P(a)_{h-1,d-1}.
     * Calculation starts with P(a)_{1,2}.
     */
    for (int d = 2; d <= diameter; d++) {
//        cout << "coming to the d = " << d << " loop" << endl;
        for (int h = 1; h <= min(n, d); h++) {
//            cout << "coming to the h = " << h << " loop" << endl;
            for (int i = 0; i < V; i++) {
//                cout << "coming to the node = " << i << " loop" << endl;
                a = &(nodes[i]);
//                cout << "current node id is " << i << endl;
                sum = 0.0;

                if (!a->validity) continue;

                if (h == 1) {
                    for (int in = 0; in < 4; in++) {
                        if (a->neighbors[in] >= V) {  // Check array bounds
                            std::cerr << "Neighbor index out of bounds: " << a->neighbors[in] << std::endl;
                            continue;
                        }
                        neighbor = &nodes[a->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                        double f = !hasFaultyLink(a, neighbor) ? 1.0 : 0.0;
                        double temp = getProbability(neighbor, h, d - 1);
                        sum += f * temp;
                    }
                } else if (h == d) {
                    for (int in = 0; in < 4; in++) {
                        if (a->neighbors[in] >= V) {  // Check array bounds
                            std::cerr << "Neighbor index out of bounds: " << a->neighbors[in] << std::endl;
                            continue;
                        }
                        neighbor = &nodes[a->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                        double f = !hasFaultyLink(a, neighbor) ? 1.0 : 0.0;
//                                cout << "faultiness = " << f << endl;
                        double temp = getProbability(neighbor, h - 1, d - 1);
//                                cout << "temp = " << temp << endl;
                        sum += f * temp;
                    }
                } else {
                    for (int in = 0; in < 4; in++) {
                        if (a->neighbors[in] >= V) {
                            std::cerr << "Neighbor index out of bounds: " << a->neighbors[in] << std::endl;
                            continue;
                        }
                        neighbor = &nodes[a->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                        double f = !hasFaultyLink(a, neighbor) ? 1.0 : 0.0;
                        double p1 = getProbability(neighbor, h - 1, d - 1);
                        double p2 = getProbability(neighbor, h, d - 1);
                        sum += (f * (double(h - 1) * p1 + double(d - h) * p2)) / double(d - 1);
                    }
                }
                sum /= 4;
//                p /= pow(2, h) * combination(n, h);
                setProbability(a, h, d, sum);
//                cout << "p_" << h << "," << d << " = " << p << endl;
            }
        }

//    cout << "Finished calculation of routing probabilities." << endl;
    }
}

void metaTorus::calcDirectedRoutingProbabilities() {
    Node *a;
    Node *s;
    Node *neighbor;
    double sum, p, p11;

//    cout << "Starting calculation of directed routing probabilities..." << endl;
    /*
     * まず全てのノードについてP(a)_{1,1}を計算
     * First, calculate P(a)_{1,1} for all nodes
     */
    for(int i=0; i<V; i++) {
//        std::cout << "Processing node " << i << std::endl;  // Debugging output
        s = &nodes[i];

        if (!s->validity) continue;

        for (int j=0; j<4; j++) {
            if (s->neighbors[j] < 0 || s->neighbors[j] >= V) {  // Check array bounds
                std::cerr << "Neighbor index out of bounds: " << s->neighbors[j] << std::endl;
                continue;
            }
            double total = 0.0;
            a = &nodes[s->neighbors[j]];
//            cout << "getting direction node " << a->index << " for h = d = 1" << endl; // Debugging output
            sum = 0.0;

            for (int in=0; in < 4; in++) {
                if (s->neighbors[in] < 0 || s->neighbors[j] >= V) {  // Check array bounds
                    std::cerr << "Neighbor index out of bounds: " << s->neighbors[in] << std::endl;
                    continue;
                }

                neighbor = &nodes[s->neighbors[in]];

//                cout << "getting neighbor node " << neighbor->index << " for h = d = 1" << endl;  // Debugging output
                if (a != neighbor){
                    total += 1.0;
                    sum += !hasFaultyLink(s, a) ? 1.0 : 0.0;
//                    cout << "faultiness = " << hasFaultyLink(s, a) << endl;
                }
            }

            p11 = sum / total;
//            cout << "p1,1 = " << p11 << " for node " << s->index << " -> neighbor " << a->index << endl;
            setDirectedProbability(s, a, 1, 1, p11);
        }

    }

    /*
     * P(a)_{h,d}を計算していく
     * P(a)_{h,d}の計算には, P(a)_{h,d-1}やP(a)_{h-1,d-1}が必要になる.
     * P(a)_{1,1}だけで計算できるP(a)_{1,2}から順番に計算していく(動的計画法)
     *
     * Calculation of P(a)_{h,d}
     * P(a)_{h,d} needs P(a)_{h,d-1} and P(a)_{h-1,d-1}.
     * Calculation starts with P(a)_{1,2}.
     */

//    cout << "come to the second round" << endl;
    for (int d = 2; d <= diameter; d++) {
//        cout << "coming to the d = " << d << " loop" << endl;
        for (int h = 1; h <= min(n, d); h++) {
//            cout << "coming to the h = " << h << " loop" << endl;
            for (int i = 0; i < V; i++) {
//                cout << "coming to the node = " << i << " loop" << endl;
                s = &(nodes[i]);

                if (!s->validity) continue;

//                std::cout << "Processing node " << i << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
//                cout << "current node id is " << i << endl;
                for (int j = 0; j < 4; j++) {
//                    if (j >= sizeof(a->neighbors)) {
//                        cerr << "Error: Neighbor index " << j << " out of bounds for node " << i << endl;
//                        continue;
//                    }
                    if (s->neighbors[j] >= V) {  // Check array bounds
                        std::cerr << "Neighbor index out of bounds: " << s->neighbors[j] << std::endl;
                        continue;
                    }
                    double total = 0.0;
                    //make the neighbor a node such that we an check the faultiness
                    a = &nodes[s->neighbors[j]];
//                    cout << "getting direction node " << a->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                    sum = 0.0;
                    if (h == 1){
                        for (int in = 0; in < 4; in++) {
                            if (s->neighbors[in] >= V) {  // Check array bounds
                                std::cerr << "Neighbor index out of bounds: " << s->neighbors[in] << std::endl;
                                continue;
                            }
                            neighbor = &nodes[s->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                            if (a != neighbor){
                                total += 1.0;
                                double f = !hasFaultyLink(s, neighbor) ? 1.0 : 0.0;
                                double temp = getDirectedProbability(neighbor, s, h, d - 1);
                                sum += f * temp;
                            }
                        }
                    } else if (h == d){
                        for (int in = 0; in < 4; in++) {
                            if (s->neighbors[in] >= V) {  // Check array bounds
                                std::cerr << "Neighbor index out of bounds: " << s->neighbors[in] << std::endl;
                                continue;
                            }
                            neighbor = &nodes[s->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                            if (a != neighbor){
                                total += 1.0;
//                                cout << "not direction node! for h = " << h << ", d = " << d << std::endl;
                                double f = !hasFaultyLink(s, neighbor) ? 1.0 : 0.0;
//                                cout << "faultiness = " << f << endl;
                                double temp = getDirectedProbability(neighbor, s, h - 1, d - 1);
//                                cout << "temp = " << temp << endl;
                                sum += f * temp;
                            }
                        }
                    } else {
                        for (int in = 0; in < 4; in++) {
//                            if (s->neighbors[in] >= V) {
//                                std::cerr << "Neighbor index out of bounds: " << s->neighbors[in] << std::endl;
//                                continue;
//                            }
                            neighbor = &nodes[s->neighbors[in]];
//                            cout << "getting neighbor node " << neighbor->index << " for h = " << h << ", d = " << d << std::endl;  // Debugging output
                            if (a != neighbor){
                                total += 1.0;
                                double f = !hasFaultyLink(s, neighbor) ? 1.0 : 0.0;
                                double p1 = getDirectedProbability(neighbor, s, h - 1, d - 1);
                                double p2 = getDirectedProbability(neighbor, s, h, d - 1);
                                sum += (f * (double(h - 1) * p1 + double(d - h) * p2)) / double(d - 1);
                                // how to
                            }
                        }

                    }
                    p = sum / total;
//                    cout << "here we are at the final calculation = " << p << endl;
//                    cout << " s = " << s->index
//                            << ", a = " << a->index
//                            << ", h = " << h
//                            << ", d = " << d << endl;
                    setDirectedProbability(s, a, h, d, p);
//                    cout << " got p = " << getDirectedProbability(s,a,h,d) << endl;
                }
            }
        }
    }

//    cout << "Finished calculation of routing probabilities." << endl;
}

void metaTorus::testDirectedRoutingProbabilities() {
    Node *a;
    Node *s;
    Node *neighbor;
    double sum, p, p11;

//    cout << "Starting calculation of directed routing probabilities..." << endl;
    /*
     * まず全てのノードについてP(a)_{1,1}を計算
     * First, calculate P(a)_{1,1} for all nodes
     */
    for(int i=0; i<V; i++) {
//        std::cout << "Processing node " << i << std::endl;  // Debugging output
        s = &nodes[i];
        for (int j=0; j<4; j++) {
            if (s->neighbors[j] < 0 || s->neighbors[j] >= V) {  // Check array bounds
                std::cerr << "Neighbor index out of bounds: " << s->neighbors[j] << std::endl;
                continue;
            }
            double total = 0.0;
            a = &nodes[s->neighbors[j]];
//            cout << "getting direction node " << a->index << " for h = d = 1" << endl; // Debugging output
            sum = 0.0;

            for (int in=0; in < 4; in++) {
                if (s->neighbors[in] < 0 || s->neighbors[j] >= V) {  // Check array bounds
                    std::cerr << "Neighbor index out of bounds: " << s->neighbors[in] << std::endl;
                    continue;
                }

                neighbor = &nodes[s->neighbors[in]];

//                cout << "getting neighbor node " << neighbor->index << " for h = d = 1" << endl;  // Debugging output
                if (a != neighbor){
                    total += 1.0;
                    sum += !hasFaultyLink(s, a) ? 1.0 : 0.0;
//                    cout << "faultiness = " << hasFaultyLink(s, a) << endl;
                }
            }

            cout << "sum = " << sum;

            p11 = sum / total;
            cout << "p1,1 = " << p11 << " for node " << s->index << " -> neighbor " << a->index << endl;
            setDirectedProbability(s, a, 1, 1, p11);
        }

    }

    /*
     * P(a)_{h,d}を計算していく
     * P(a)_{h,d}の計算には, P(a)_{h,d-1}やP(a)_{h-1,d-1}が必要になる.
     * P(a)_{1,1}だけで計算できるP(a)_{1,2}から順番に計算していく(動的計画法)
     *
     * Calculation of P(a)_{h,d}
     * P(a)_{h,d} needs P(a)_{h,d-1} and P(a)_{h-1,d-1}.
     * Calculation starts with P(a)_{1,2}.
     */

//    cout << "come to the second round" << endl;

//    cout << "Finished calculation of routing probabilities." << endl;
}

void metaTorus::setProbability(Node *a, int h, int d, double p) {
    assert(0 < h && h <= n);
    assert(0 < d && d <= diameter);
    assert(0.0 <= p && p <= 1.0);
    assert(a->validity);

    // Ensure this is the first time the probability is being set
    assert(nodes[a->index].P[h][d] == -1);
    nodes[a->index].P[h][d] = p;
}

void metaTorus::setDirectedProbability(Node *a, Node *b, int h, int d, double d_p) {
    assert(0 <= h && h <= n);
    assert(0 <= d && d <= diameter);
    assert(0.0 <= d_p && d_p <= 1.0);
    assert(a->validity);
    assert(b->validity);

    int index = getSpecificNeighbor(a, b->index);
    // Ensure that the index is within bounds
    assert(index >= 0 && index < 4);

    // Ensure this is the first time the directed probability is being set
//    assert(nodes[a->index].d_P[h][d][index] == -1);
    nodes[a->index].d_P[h][d][index] = d_p;
}

double metaTorus::getProbability(Node *a, int h, int d) {
    assert(0 <= h && h <= n);
    assert(0 <= d && d <= diameter);
    assert(a->validity);

    if (h == 0) {
        assert(d == 0);
        return 1.0; // Base case: self-loop probability
    }

    // Ensure the probability was previously set
    if (nodes[a->index].P[h][d] == -1) {
        std::cerr << "Invalid P value detected! Node: " << a->index
                  << ", Hamming distance: " << h
                  << ", Distance: " << d << std::endl;
    }
    assert(nodes[a->index].P[h][d] != -1);

    return nodes[a->index].P[h][d];
}

double metaTorus::getDirectedProbability(Node *a, Node *b, int h, int d) {
    assert(0 <= h && h <= n);
    assert(0 <= d && d <= diameter);
    assert(a->validity);
    assert(b->validity);

    if (h == 0) {
        assert(d == 0);
        return 1.0; // Base case: self-loop probability
    }

    int index = getSpecificNeighbor(a, b->index);
    // Check if the index is within bounds
    assert(index >= 0 && 4);

    // Ensure the directed probability was previously set
    assert(nodes[a->index].d_P[h][d][index] != -1);

    return nodes[a->index].d_P[h][d][index];
}

void metaTorus::printProbabilities() {
    Node* a;
    std::cout << "Node ";
    for (int d = 1; d <= floor(k / 2) * n; d++)
        for (int h = 1; h <= std::min(n, d); h++) {
            std::cout << "P_{" << h << "," << d << "} ";
        }
    std::cout << std::endl;

    for (int i = 0; i < V; i++) {
        a = &nodes[i];

        if (!a->validity) continue;

        std::cout << "(";
        for (int j = 0; j < n; j++) {
            std::cout << a->value[j];
            if (j < n - 1)
                std::cout << ",";
        }
        std::cout << ") "; // Print the node coordinate

        for (int d = 1; d <= floor(k / 2) * n; d++)
            for (int h = 1; h <= std::min(n, d); h++) {
                std::cout << int(getProbability(a, h, d) * 100) / 100.0 << " ";
            }
        std::cout << std::endl;
    }
}

void metaTorus::printDirectedProbabilities() {
    Node* s;
    Node* a;
    std::cout << "Node ";
    for (int d = 1; d <= floor(k / 2) * n; d++)
        for (int h = 1; h <= std::min(n, d); h++) {
            std::cout << "P_{" << h << "," << d << "} ";
        }
    std::cout << std::endl;

    for (int i = 0; i < V; i++) {
        s = &nodes[i];

        if (!s->validity) continue;

        for (int j = 0; j < s->neighbors.size(); j++) {
            a = &nodes[s->neighbors[j]];

            std::cout << "to (";
            for (int q = 0; q < n; q++) {
                std::cout << s->value[q];
                if (q < n - 1)
                    std::cout << ",";
            }
            std::cout << ") = " << s->index;
            std::cout << " from (";
            for (int p = 0; p < n; p++) {
                std::cout << a->value[p];
                if (p < n - 1)
                    std::cout << ",";
            }
            std::cout << ") = " << a->index << " ";

            for (int d = 1; d <= floor(k / 2) * n; d++)
                for (int h = 1; h <= std::min(n, d); h++) {
                    std::cout << std::fixed << std::setprecision(2)
                              << getDirectedProbability(s, a, h, d) << " ";
                }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

int metaTorus::route(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    int dist_ct, h;
    double p, p_max;
    Node *neighbor, *max_neighbor = nullptr;

//    cout << "Current node: " << c->index << ", Distance: " << d << endl;

    // Base case
    if(c->index == t->index){
//        cout << "Reached dest node in route: " << t->index << endl;
        return d;
    }
    //path length?

    dist_ct = distance(c, t);
    p_max = -1;
    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(isNeighbor(c, neighbor));

        if(prev and prev->index == neighbor->index) continue;

        if (inPre(neighbor, c, t) and !hasFaultyLink(c, neighbor) ) {
            h = hammingDistance(neighbor, t);
            p = getProbability(neighbor, h, std::max(0, dist_ct - 1));

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(max_neighbor != nullptr && p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
//        cout << "Preferred route chosen: " << max_neighbor->index << endl;
        return route(c, max_neighbor, t, visited, d + 1);
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(isNeighbor(c, neighbor));
        //haslink
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (inSpr(neighbor, c, t) and !hasFaultyLink(c, neighbor)) {
            h = hammingDistance(neighbor, t);
            p = getProbability(neighbor, h, min(dist_ct + 1, diameter));

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(max_neighbor != nullptr && p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
//        cout << "Spare route chosen: " << max_neighbor->index << endl;
        return route(c, max_neighbor, t, visited, d + 1);
    }

//    cout << "No valid route found from node " << c->index << " to node " << t->index << endl;
    return DELIVERY_FAIL;
}

int metaTorus::directed_route(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    int dist_ct, h;
    double p, p_max;
    Node *neighbor, *max_neighbor = nullptr;

//    cout << "in directed Current node: " << c->index << ", Distance: " << d << endl;

    // Base case
    if(c->index == t->index){
//        cout << "Reached dest node in route: " << t->index << endl;
        return d;
    }

    dist_ct = distance(c, t);
    p_max = -1;
    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(isNeighbor(c, neighbor));
        //has link
        if(prev and prev->index == neighbor->index){
            continue;
            //skip this neighbor
        }
        if (inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor) ) {
            h = hammingDistance(neighbor, t);
            p = getDirectedProbability(neighbor, c, h, dist_ct - 1);

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(max_neighbor != nullptr && p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
//        route_p.push_back(max_neighbor->index);
//        cout << " directed route choose" << max_neighbor->index << ", " ;
        return directed_route(c, max_neighbor, t, visited, d + 1);
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        assert(isNeighbor(c, neighbor));
        //haslink
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (inSpr(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            h = hammingDistance(neighbor, t);
            p = getDirectedProbability(neighbor, c, h, min(dist_ct + 1, diameter));

            if (p > p_max) {
                p_max = p;
                max_neighbor = neighbor;
            }
        }
    }
    if(max_neighbor != nullptr && p_max > 0) {
        if (visited.count(max_neighbor->index) > 0)
            return DELIVERY_FAIL;
//        cout << " directed route choose" << max_neighbor->index << ", ";
        return directed_route(c, max_neighbor, t, visited, d + 1);
    }

    return DELIVERY_FAIL;
}

int metaTorus::brute(Node *prev, Node *c, Node *t, std::unordered_map<int, bool> visited, int d) {
    Node *neighbor;

    //std::cout << "c:";
    //for(int i=0; i<n; i++)
    //    std::cout << c->value[i] << " ";
    //std::cout << std::endl;

    if(c->index == t->index){
        return d;
    }

    visited[c->index] = true;

    // Try to deliver message to preferred node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        if(prev and prev->index == neighbor->index){
            continue;
        }

        if(inPre(neighbor, c, t) and not hasFaultyLink(c, neighbor)) {
            if(visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
//            cout << "brute choose " << neighbor->index << ",";
            return brute(c, neighbor, t, visited, d + 1);
        }
    }

    // Try to deliver message to spare node
    for(int i=0; i<2*n; i++) {
        neighbor = &(nodes[c->neighbors[i]]);
        if(prev and prev->index == neighbor->index){
            continue;
        }
        if (not hasFaultyLink(c, neighbor)){
            assert(inSpr(neighbor, c, t));
            if (visited.count(neighbor->index) > 0)
                return DELIVERY_FAIL;
//            cout << "brute choose " << neighbor->index << ",";
            return brute(c, neighbor, t, visited, d + 1);
        }
    }

    return DELIVERY_FAIL;
}

int metaTorus::bfs(Node *c, Node *t, std::unordered_map<int, bool> visited) {
    std::queue<int> q_i, q_d;
    int a_d;
    Node *a, *neighbor;

    if(c->index == t->index){
        return 0;
    }

    visited[c->index] = true;
    q_i.push(c->index);
    q_d.push(0);

    while(not q_i.empty()){

        assert(q_d.size() == q_i.size());

        a = &(nodes[q_i.front()]);
        q_i.pop();
        a_d = q_d.front();
        q_d.pop();

        if(a->index == t->index)
            return a_d;

        for(int i=0; i<2*n; i++) {
            neighbor = &(nodes[a->neighbors[i]]);
            if(not hasFaultyLink(a, neighbor) and visited.count(neighbor->index) == 0){
                visited[neighbor->index] = true;
                q_i.push(neighbor->index);
                q_d.push(a_d+1);
            }
        }
    }

    assert(q_d.empty());
    return DELIVERY_FAIL;
}

void metaTorus::printPath(vector<int> path) {
    std::cout << "Path: ";
    for(const int& n : path) {
        std::cout << n << " ";
    }
    cout << endl;
}

metaTorus::~metaTorus() {
//    adjacency.clear();
    F.clear();
    nodes.clear();

//    for(int i=0; i<V; i++){
//        delete[] nodes[i].value;
//        for(int j=0; j<n; j++){
//            delete[] nodes[i].P[j];
//            for (int k = 0; k < 2*n; ++k) {
//                delete[] nodes[i].d_P[j][k]; // Deallocate memory for neighbor
//            }
//            delete[] nodes[i].d_P[j];
//        }
//        free(nodes[i].P);
//        free(nodes[i].d_P);
//        delete[] nodes[i].neighbors;
//    }
//    free(nodes);
}

//
// Created by Carina Z on 2024/11/17.
//

