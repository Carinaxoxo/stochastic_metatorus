#include <random>

#ifndef STOCHASTICMETATORUS_RANDOM_H
#define STOCHASTICMETATORUS_RANDOM_H

// nがシード、mがmod数
// 0 <= r <= m-1の整数値を返す
static int get_random(int n, int m) {
    // シードが同じなら同じ乱数が出る
    std::mt19937 mt(n);

    // moduloの挙動がc++とpythonで違うので、念の為ここで吸収しておく
    return (int(mt())%m+m)%m;
}

#endif //STOCHASTICMETATORUS_RANDOM_H

//
// Created by Carina Z on 2024/11/17.
//
