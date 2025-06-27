#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
using namespace std;

#ifndef STOCHASTICMETATORUS_METATORUS_H
#define STOCHASTICMETATORUS_METATORUS_H

typedef struct NODE
{
    int index;
//    int *value; // n-dimension vector
    vector<int> value;
    bool validity;
//    double **P;
    vector<vector<double>> P;
//    unordered_map<double, unordered_map<int, int>> **d_P;
//    double ***d_P;
    vector<vector<vector<double>>> d_P;
    //a matrix, with respect to each neighbor
    /*
     * P: matrix of probabilities
     * d_P: matrix of directed probabilities
     * P contain probability of node depends on h and d
     * where 1 <= h <= n, 1 <= d <= diameter
     * d_p contain the probabilities according
     * to specific neighbor
    */

//    int *neighbors;  // indices
    vector<int> neighbors;
} Node;

//struct Link {
//    int a;
//    int b;
//};


// using namespace boost::numeric::ublas;
class metaTorus {

public:

    // ↓ ====== called in constructor
    void allocate();
    void setValues();
    void createConnection();
    void saveNeighbor();
    // ↑ ====== called in constructor

    /*
     * isNeighborとhasLinkは全く同じ挙動になるが、
     * hasLinkはisNeighborの結果を利用する(=adjacencyを参照する)ので速い
     *
     * Function isNeighbor and hasLink has same behavior
     * hasLink is faster than isNeighbor (hasLink uses result of isNeighbor)
     */
    int getId(vector<int> input);
    bool isNeighbor(Node *a, Node *b);
    int getSpecificNeighbor(Node *a, int neighborIndex);

    // ↓ ====== just to read/write data
//    bool hasLink(Node *a, Node *b);
    void setFaultyLink(Node *a, Node *b);
    bool hasFaultyLink(Node *a, Node *b);

    void setProbability(Node *a, int h, int d, double p);
    double getProbability(Node *a, int h, int d);
    void setDirectedProbability(Node *a, Node *b, int h, int d, double p);
    double getDirectedProbability(Node *a, Node *b, int h, int d);
    // ↑ ====== just to read/write data

    // Functions to calculate distance
    int distance(Node* a, Node* b);
    int hammingDistance(Node* a, Node* b);
    int getDistance(Node *a, Node *b);
    int getHops(Node *a, Node *b);

    // Functions to classify neighbor nodes
    bool inPre(Node* neighbor, Node* a, Node* b);
    bool inSpr(Node* neighbor, Node* a, Node* b);

    int k;
    int n;

//    Node *nodes; // nodes
    int V; // number of nodes
    int E; // number of edges
    int diameter;

    // スパース行列
//    mapped_matrix<int> adjacency;
    //mapped_matrix<int> F;
//    vector<pair<int, int>> F;
    vector<Node> nodes;
    vector< vector<int> > F;
    //fina a way to replace mapped_matrix

    // n dimensional vector (string) -> index of a node (in nodes)
    unordered_map<string, int> node_value_index;

    metaTorus(int n, int k); // constructor
    ~metaTorus(); // destructor

    // for debugging
    void display(); // Declaration of display method
    void setTestedFaultyLinks();
    void printAdjMatrix();
    void printFaultyLinks();
    void printProbabilities();
    void printDirectedProbabilities();
    vector<vector<int>> getNeighbors(Node *a);

    // ↓割りと重要な関数たち Very important algorithms
    void setRandomFaultyLinks(double p_faulty, int* seed);
    //we need to sort it out where is it
    //bool isFaulty(Node *a, Node *b);

    // algorithm P
    void calcRoutingProbabilities();
    void calcDirectedRoutingProbabilities();
    void testDirectedRoutingProbabilities();

    //new strategy
    std::vector<std::vector<std::vector<int>>> return_node_which_close_to_goal(Node *c, Node *t);

    // testing routing algorithms
    int route_test(Node *prev, Node* c, Node* t, std::unordered_map<int, bool>& visited, int depth);
    int directed_route_test(Node *prev, Node *c, Node *t, std::unordered_map<int, bool>& visited, int depth);

    // routing algorithms
    int route(Node *prev, Node* c, Node* t, unordered_map<int, bool> visited, int d);
    int directed_route(Node *prev, Node *c, Node *t, unordered_map<int, bool> visited, int d);
    int strategy_route(Node *prev, Node *c, Node *t, unordered_map<int, bool> visited, int d);
    int brute(Node *prev, Node *c, Node *t, unordered_map<int, bool> visited, int d);

    // breath-first search for optimal routing
    int bfs(Node *c, Node *t, unordered_map<int, bool> visited);

    // test if the graph is fully connected
    bool isFullyConnected();
    void clearFaultyLinks();

    void printPath(vector<int> p);
};

void sortCouple(double *a, double *b, int len);

#endif //STOCHASTICMETATORUS_METATORUS_H

//
// Created by Carina Z on 2024/11/17.
//
