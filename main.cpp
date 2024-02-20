#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "OrdinaryDiffSolve.h"
#include "VectorOperations.h"

using namespace std;

const string INPUT_FILE_NAME = "C:\\Users\\Alex\\Documents\\Qt\\BMSTU\\NumMethods Semester 6\\NumMethods_Sem6_Lab_1\\input.txt";
const string METHOD_PARAMS_FILE_NAME = "C:\\Users\\Alex\\Documents\\Qt\\BMSTU\\NumMethods Semester 6\\NumMethods_Sem6_Lab_1\\parameters.txt";
const string OUTPUT_FILE_NAME = "C:\\Users\\Alex\\Documents\\Qt\\BMSTU\\NumMethods Semester 6\\NumMethods_Sem6_Lab_1\\output.txt";

const double PI = acos(-1);
const double PI_SQR = PI * PI;

template<typename T>
inline vector<T> func1 (T, const vector<T>& xs) {
    auto x = xs[0];
    auto y = xs[1];

    return {2 * x + y * y - 1, 6 * x - y * y + 1};
}

template<typename T>
inline vector<T> func2 (T, const vector<T>& xs) {
    auto x = xs[0];
    auto y = xs[1];

    return {1 - x*x - y*y , 2 * x};
}

template<typename T>
inline vector<T> func3 (T, const vector<T>& xs) {
    auto x = xs[0];
    auto y = xs[1];
    auto z = xs[2];

    return {10 * (y - x) , x * (28 - z) - y, x * y - (8/3) * z};
}

template<typename T>
inline vector<T> pendulumFunc (T, const vector<T>& xs) {
    auto x1 = xs[0];
    auto x2 = xs[1];

    return {x2 , -x1};
}

/*
 * t0 = 0.1
 * T = 1
 * u0 = {0, 15.708}
 */

template<typename T>
inline vector<T> xSqrSinFunc (T t, const vector<T>& xs) {
    auto x1 = xs[0];
    auto x2 = xs[1];

    auto t_sqr = t*t;

    return {x2 , -(2 + PI_SQR / (4 * t_sqr)) / t_sqr * x1 + 2 / t * x2};
}

// Хранение функций через метапрограммирование
template<typename T>
struct FuncData {
    static vector<rightSideFunc<T>> rightSides;
};

template<typename T>
vector<rightSideFunc<T>> FuncData<T>::rightSides = {
    func1<T>,
    func2<T>,
    func3<T>,
    pendulumFunc<T>,
    xSqrSinFunc<T>
};

template<typename Type>
InitialValueProblemParams<Type> readInputFile(){
    ifstream input(INPUT_FILE_NAME);

    InitialValueProblemParams<Type> params;
    int f_id;

    input >> f_id >> params.n;
    params.u0 = readVector<double>(input, params.n);
    input >> params.t0 >> params.T;

    input.close();

    params.f = FuncData<Type>::rightSides[f_id];

    return params;
}

template<typename Type>
void runCalculations(const InitialValueProblemParams<Type>& problemParams) {
    ifstream parametersInput(METHOD_PARAMS_FILE_NAME);

    DiffPrecisionParams<Type> precisionParams;
    int method_id;

    parametersInput >> method_id >> precisionParams.tau0 >> precisionParams.tol >> precisionParams.eps;
    precisionParams.vNorm = norm_inf;

    parametersInput.close();

    if (method_id == 0){
        RungeKuttaParams<Type> rkParams(2, 2);
        rkParams.a = {0, 0.5};
        rkParams.b = {
            {},
            {0.5}
        };
        rkParams.sigma = {0, 1};

        rungeKuttaMethod(problemParams, rkParams, precisionParams, OUTPUT_FILE_NAME);

    } else if (method_id == 1) {
        RungeKuttaParams<double> rkParams(4, 4);
        rkParams.a = {0, 0.5, 0.5, 1};
        rkParams.b = {
            {},
            {0.5},
            {0, 0.5},
            {0, 0, 1}
        };
        rkParams.sigma = {1./6, 2./6, 2./6, 1./6};

        rungeKuttaMethod(problemParams, rkParams, precisionParams, OUTPUT_FILE_NAME);

    } else if (method_id == 2) {
        explicitEulerMethod(problemParams, precisionParams, OUTPUT_FILE_NAME);

    } else if (method_id == 3) {
        implicitEulerMethod(problemParams, precisionParams, OUTPUT_FILE_NAME);

    }
}

int main() {
    /*
    InitialValueProblemParams<double> problemParams;
    problemParams.f = func1;
    problemParams.T = 2;
    problemParams.t0 = 0;
    problemParams.u0 = {0, 0};
    problemParams.n = 2;

    DiffPrecisionParams<double> precisionParams;
    precisionParams.eps = 1e-2;
    precisionParams.tau0 = 1e-1;
    precisionParams.tol = 16;
    precisionParams.vNorm = norm_inf;

    RungeKuttaParams<double> rkParams(4, 4);
    rkParams.a = {0, 0.5, 0.5, 1};
    rkParams.b = {
        {},
        {0.5},
        {0, 0.5},
        {0, 0, 1}
    };
    rkParams.sigma = {1./6, 2./6, 2./6, 1./6};

    vector<double> a = {1, 2};
    vector<double> b = {3, 4};

    a += b;
    auto c = 2.0 * b;

    rungeKuttaMethod(problemParams, rkParams, precisionParams, OUTPUT_FILE_NAME);

    */

    auto problemParams = readInputFile<double>();

    runCalculations(problemParams);

    cout << "END";

    #ifdef QT_DEBUG
        system("pause");
    #endif

    return 0;
}
