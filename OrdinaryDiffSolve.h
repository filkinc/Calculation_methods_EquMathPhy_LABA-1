#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "VectorOperations.h"
#include "Norma.h"

using namespace std;

template <typename T>
using rightSideFunc = vector<T>(*)(T, const vector<T>&);

template<typename T>
struct DiffPrecisionParams {
    T tau0;
    T eps;
    T tol;
    vectorNorm<T> vNorm;
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
    int p;

    vector<T> a;
    vector<vector<T>> b;
    vector<T> sigma;

    RungeKuttaParams(int stages, int p) {
        this->stages = stages;
        this->p = p;

        a = vector<T>(stages);
        b = vector<vector<T>>(stages);
        sigma = vector<T>(stages);

        for (int i = 1; i < stages; ++i) {
            b[i] = vector<T>(i);
        }
    }
};

template<typename T>
vector<T> rungeKuttaStep(const RungeKuttaParams<T>& rkParams, T tau, int n, rightSideFunc<T> f, T t, const vector<T>& y) {
    int m = rkParams.stages;
    vector<vector<T>> k(m);

    auto& a = rkParams.a;
    auto& b = rkParams.b;
    auto& sigma = rkParams.sigma;

    k[0] = f(t, y);

    vector<T> kSum = vector<T>(n, 0);
    for (int i = 1; i < m; i++) {
        kSum.assign(n, 0);
        for (int j = 0; j < i; ++j) {
            kSum += b[i][j] * k[j];
        }

        k[i] = f(t + a[i] * tau, y + tau * kSum);
    }

    kSum.assign(n, 0);
    for (int j = 0; j < m; ++j) {
        kSum += sigma[j] * k[j];
    }

    auto nextY = y + tau * kSum;

    return nextY;
}


template<typename T>
void rungeKuttaMethod(const InitialValueProblemParams<T>& problemParams,
                      const RungeKuttaParams<T>& rkParams,
                      const DiffPrecisionParams<T>& precisionParams,
                      const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto prevY = problemParams.u0;

    T t = problemParams.t0;

    output << t;
    for (auto& yi : prevY) {
        output << ' ' << yi;
    }
    output << '\n';

    T errorCoef = pow(2, rkParams.p) / (pow(2, rkParams.p) - 1);

    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> medY;
        vector<T> altNextY;
        vector<T> nextY;

        nextY = rungeKuttaStep<T>(rkParams, tau, problemParams.n, problemParams.f, t, prevY);

        medY = rungeKuttaStep<T>(rkParams, tau / 2, problemParams.n, problemParams.f, t, prevY);
        altNextY = rungeKuttaStep<T>(rkParams, tau / 2, problemParams.n, problemParams.f, t, medY);

        T error = precisionParams.vNorm(errorCoef * (altNextY - nextY));

        if (error <= precisionParams.eps) {
            t += tau;

            output << t;
            for (auto& yi : nextY) {
                output << ' ' << yi;
            }
            output << '\n';

            prevY = nextY;

            if (error <= precisionParams.eps / precisionParams.tol) {
                tau *= 2;
            }

        } else {
            tau /= 2;
        }
    }

    output.close();
}
