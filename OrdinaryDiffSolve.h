#pragma once

#include <vector>
#include <string>
#include <fstream>

#include "VectorOperations.h"

using namespace std;

template <typename T>
using rightSideFunc = vector<T>(*)(T, const vector<T>&);

template<typename T>
struct DiffPrecisionParams {
    T tau0;
    T eps;
};

template<typename P>
struct InitialValueProblemParams {
    int n;
    rightSideFunc<P> f;
    vector<P> u0;
    P t0;
    P T;
};

template<typename T>
struct RungeKuttaParams {
    int stages;
    vector<T> a;
    vector<vector<T>> b;
    vector<T> sigma;

    RungeKuttaParams(int stages) {
        this->stages = stages;
        // (ToDo)
    }
};


template<typename T>
void rungeKuttaMethod(const InitialValueProblemParams<T>& problemParams,
                      const RungeKuttaParams<T>& rkParams,
                      const DiffPrecisionParams<T>& precisionParams,
                      const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto prevY = problemParams.u0;
    auto

}

template<typename T>
vector<T> rungeKuttaStep(const RungeKuttaParams<T>& rkParams, T tau, int n, rightSideFunc<T> f, T t, T y) {
    int m = rkParams.stages;
    vector<vector<T>> k(m);

    auto& a = rkParams.a;
    auto& b = rkParams.b;
    auto& sigma = rkParams.sigma;

    k[0] = f(t, y);

    vector<T> kSum = 0;
    for (int i = 1; i < m; i++) {
        kSum = vector<T>(n, 0);
        for (int j = 0; j < i; ++j) {
            kSum += b[i][j] * k[j];
        }

        k[i] = f(t + a[i] * tau, y + tau * kSum);
    }

    kSum = vector<T>(n, 0);
    for (int j = 0; j < m; ++j) {
        kSum += sigma[j] * k[j];
    }

    auto nextY = y + tau * kSum;

    return nextY;
}



