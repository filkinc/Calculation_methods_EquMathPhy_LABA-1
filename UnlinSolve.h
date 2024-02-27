#pragma once

#include <cmath>
#include <vector>
#include "VectorOperations.h"
#include "QuadMatrix.h"
#include "LinSolveAlgs.h"

using namespace std;

template<typename T>
using UnlinSolveSystemFunc = vector<T>(*)(const vector<T>&);
template<typename T>
using UnlinSolveFunc = T(*)(T);

template<typename T>
struct UnlinSolveOptions {
    vectorNorm<T> vNorm;
    matrixNorm<T> mNorm;
    T eps;
    int maxIterationCount;
};

// Метод Ньютона!

template<typename T>
T diff12(UnlinSolveFunc<T> f, T x, T tau){
    return (3 * f(x + tau) - 4 * f(x) + f(x - tau)) / (2 * tau);
}

template<typename T>
vector<T> diff12(UnlinSolveSystemFunc<T> f, const vector<T>& x, int coordInd, T tau){
    vector<T> fX = x;
    vector<T> bX = x;

    fX[coordInd] += tau;
    bX[coordInd] -= tau;

    return (3 * f(fX) - 4 * f(x) + f(bX)) / (2 * tau);
}

template<typename T>
QuadMatrix<T> jacobian12(UnlinSolveSystemFunc<T> f, const vector<T>& x, T tau){
    int n = x.size();

    QuadMatrix<T> J(n);

    for (int j = 0; j < n; j++) {
        auto df = diff12(f, x, j, tau);

        for (int i = 0; i < n; ++i) {
            J(i, j) = df[i];
        }
    }

    return J;
}

template<typename T, typename FunctorType>
T templateDiff12(FunctorType f, T x, T tau){
    return (T(3) * f(x + tau) - T(4) * f(x) + f(x - tau)) / (2 * tau);
}

template<typename T, typename FunctorType>
vector<T> templateDiff12(FunctorType f, const vector<T>& x, int coordInd, T tau){
    vector<T> fX = x;
    vector<T> bX = x;

    fX[coordInd] += tau;
    bX[coordInd] -= tau;

    return (T(3) * f(fX) - T(4) * f(x) + f(bX)) / (2 * tau);
}

template<typename T, typename FunctorType>
QuadMatrix<T> templateJacobian12(FunctorType f, const vector<T>& x, T tau){
    int n = x.size();

    QuadMatrix<T> J(n);

    for (int j = 0; j < n; j++) {
        auto df = templateDiff12(f, x, j, tau);

        for (int i = 0; i < n; ++i) {
            J(i, j) = df[i];
        }
    }

    return J;
}

template<typename T>
vector<T> newtonMethod(const vector<T>& x0, UnlinSolveSystemFunc<T> f, const UnlinSolveOptions<T>& opts) {
    vector<T> x = x0;

    for (int iterCount = 1; iterCount <= opts.maxIterationCount && opts.vNorm(x) < opts.eps; ++iterCount) {
        auto J = jacobian12(f, x, opts.eps);
        auto y = gaussLinSolve(J, -f(x));

        x = x + y;
    }

    return x;
}

template<typename T, typename FunctorType>
vector<T> templateNewtonMethod(const vector<T>& x0, FunctorType f, const UnlinSolveOptions<T>& opts) {
    vector<T> x = x0;

    // auto err = opts.vNorm(x);

    auto d = 1 + 1;

    for (int iterCount = 1; iterCount <= opts.maxIterationCount && opts.vNorm(f(x)) >= opts.eps; ++iterCount) {
        auto J = templateJacobian12(f, x, opts.eps);
        auto y = gaussLinSolve(J, -f(x)).first;

        x = x + y;
    }

    return x;
}
