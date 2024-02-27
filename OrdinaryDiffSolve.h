#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>

#include "VectorOperations.h"
#include "UnlinSolve.h"
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
struct multiStepDiffParams {
    int steps;

    vector<T> a;
    vector<vector<T>> b;
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
                      bool useStepRegulation,
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

        if (!useStepRegulation || error <= precisionParams.eps) {
            t += tau;

            output << t;
            for (auto& yi : nextY) {
                output << ' ' << yi;
            }
            output << '\n';

            prevY = nextY;

            if (useStepRegulation && error <= precisionParams.eps / precisionParams.tol) {
                tau *= 2;
            }

        } else {
            tau /= 2;
        }
    }

    output.close();
}

template<typename T>
void explicitEulerMethod(const InitialValueProblemParams<T>& problemParams,
                         const DiffPrecisionParams<T>& precisionParams,
                         const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto prevY = problemParams.u0;
    auto f = problemParams.f;

    T t = problemParams.t0;

    output << t;
    for (auto& yi : prevY) {
        output << ' ' << yi;
    }
    output << '\n';

    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> nextY;

        nextY = prevY + tau * f(t, prevY);

        t += tau;

        output << t;
        for (auto& yi : nextY) {
            output << ' ' << yi;
        }
        output << '\n';

        prevY = nextY;
    }

    output.close();
}

template<typename T>
void implicitEulerMethod(const InitialValueProblemParams<T>& problemParams,
                         const DiffPrecisionParams<T>& precisionParams,
                         const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto prevY = problemParams.u0;
    auto f = problemParams.f;

    T t = problemParams.t0;

    output << t;
    for (auto& yi : prevY) {
        output << ' ' << yi;
    }
    output << '\n';

    UnlinSolveOptions<T> unlinOpts;
    unlinOpts.mNorm = norm_inf;
    unlinOpts.vNorm = norm_inf;
    unlinOpts.maxIterationCount = 10;
    unlinOpts.eps = precisionParams.eps;

    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> nextY;

        auto F = [prevY, t, tau, f](vector<T> y) {
            return y - tau * f(t + tau, y) - prevY;
        };

        nextY = templateNewtonMethod(prevY, F, unlinOpts);

        t += tau;

        output << t;
        for (auto& yi : nextY) {
            output << ' ' << yi;
        }
        output << '\n';

        prevY = nextY;
    }

    output.close();
}

template<typename T>
void symmetricalSchemeMethod(const InitialValueProblemParams<T>& problemParams,
                         const DiffPrecisionParams<T>& precisionParams,
                         const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto prevY = problemParams.u0;
    auto f = problemParams.f;

    T t = problemParams.t0;

    output << t;
    for (auto& yi : prevY) {
        output << ' ' << yi;
    }
    output << '\n';

    UnlinSolveOptions<T> unlinOpts;
    unlinOpts.mNorm = norm_inf;
    unlinOpts.vNorm = norm_inf;
    unlinOpts.maxIterationCount = 10;
    unlinOpts.eps = precisionParams.eps;

    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> nextY;

        auto F = [prevY, t, tau, f](vector<T> y) {
            return y - tau * (f(t + tau, y) + f(t, prevY)) / T(2) - prevY;
        };

        nextY = templateNewtonMethod(prevY, F, unlinOpts);

        t += tau;

        output << t;
        for (auto& yi : nextY) {
            output << ' ' << yi;
        }
        output << '\n';

        prevY = nextY;
    }

    output.close();
}

template<typename T>
void adamsBashfort4Method(const InitialValueProblemParams<T>& problemParams,
                         const DiffPrecisionParams<T>& precisionParams,
                         const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto f = problemParams.f;
    auto n = problemParams.n;

    T t = problemParams.t0;

    // Параметры метода Рунге-Кутты
    RungeKuttaParams<T> rkParams(4, 4);
    rkParams.a = {0, 0.5, 0.5, 1};
    rkParams.b = {
        {},
        {0.5},
        {0, 0.5},
        {0, 0, 1}
    };
    rkParams.sigma = {1./6, 2./6, 2./6, 1./6};

    vector<T> prevYs[4]; // хранит последние 4 y-ка в след. порядке:
                         // {y_n, y_{n-1}, y_{n-2}, y_{n-3}}
    prevYs[3] = problemParams.u0;

    output << t;
    for (auto& yi : prevYs[3]) {
        output << ' ' << yi;
    }
    output << '\n';

    // Поиск y1, y2, y3 методом Рунге-Кутты (параметры метода выше)
    for (int k = 2; k >= 0; --k) {
        t += tau;

        prevYs[k] = rungeKuttaStep(rkParams, tau, n, f, t, prevYs[k + 1]);
        output << t;
        for (auto& yi : prevYs[k]) {
            output << ' ' << yi;
        }
        output << '\n';
    }

    // Основная схема
    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> nextY;

        nextY = prevYs[0] + (tau / T(24)) * (T(55) * f(t, prevYs[0]) -
                T(59) * f(t - tau, prevYs[1]) + T(37) * f(t - T(2) * tau, prevYs[2]) -
                T(9) * f(t - 3 * tau, prevYs[3]));

        t += tau;

        prevYs[3] = nextY; // удаление y_{n-3} из массива / добавление y_{n+1} в массив
        // Циклический сдвиг массива вправо на 1 --- перемещает y_{n+1} в начало
        for (int i = 2; i >= 0; --i) {
            swap(prevYs[i], prevYs[i + 1]);
        }

        output << t;
        for (auto& yi : nextY) {
            output << ' ' << yi;
        }
        output << '\n';
    }

    output.close();
}

template<typename T>
void ab4PredictorCorrectorMethod(const InitialValueProblemParams<T>& problemParams,
                         const DiffPrecisionParams<T>& precisionParams,
                         const string& outputFileName) {

    ofstream output(outputFileName);

    auto tau = precisionParams.tau0;
    auto f = problemParams.f;
    auto n = problemParams.n;

    T t = problemParams.t0;

    // Параметры метода Рунге-Кутты
    RungeKuttaParams<T> rkParams(4, 4);
    rkParams.a = {0, 0.5, 0.5, 1};
    rkParams.b = {
        {},
        {0.5},
        {0, 0.5},
        {0, 0, 1}
    };
    rkParams.sigma = {1./6, 2./6, 2./6, 1./6};

    vector<T> prevYs[4]; // хранит последние 4 y-ка в след. порядке:
                         // {y_n, y_{n-1}, y_{n-2}, y_{n-3}}
    prevYs[3] = problemParams.u0;

    output << t;
    for (auto& yi : prevYs[3]) {
        output << ' ' << yi;
    }
    output << '\n';

    // Поиск y1, y2, y3 методом Рунге-Кутты (параметры метода выше)
    for (int k = 2; k >= 0; --k) {
        t += tau;

        prevYs[k] = rungeKuttaStep(rkParams, tau, n, f, t, prevYs[k + 1]);
        output << t;
        for (auto& yi : prevYs[k]) {
            output << ' ' << yi;
        }
        output << '\n';
    }

    // Основная схема
    int iterCount = 0;
    for (; t <= problemParams.T && iterCount <= 10000; ++iterCount) {
        vector<T> predictedY;
        vector<T> nextY;

        predictedY = prevYs[0] + (tau / T(24)) * (T(55) * f(t, prevYs[0]) -
                T(59) * f(t - tau, prevYs[1]) + T(37) * f(t - T(2) * tau, prevYs[2]) -
                T(9) * f(t - 3 * tau, prevYs[3]));

        nextY = prevYs[0] + (tau / T(24)) * (T(9) * f(t + tau, predictedY) +
                T(19) * f(t, prevYs[0]) - T(5) * f(t - tau, prevYs[1]) +
                f(t - T(2) * tau, prevYs[2]));

        t += tau;

        prevYs[3] = nextY; // удаление y_{n-3} из массива / добавление y_{n+1} в массив
        // Циклический сдвиг массива вправо на 1 --- перемещает y_{n+1} в начало
        for (int i = 2; i >= 0; --i) {
            swap(prevYs[i], prevYs[i + 1]);
        }

        output << t;
        for (auto& yi : predictedY) {
            output << ' ' << yi;
        }
        output << '\n';
    }

    output.close();
}

