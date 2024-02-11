#pragma once

#include <vector>

using namespace std;

template<class T>
vector<T> sum(const vector<T>& a, const vector<T>& b) {
    int n = min(a.size(), b.size());
    vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = a[i] + b[i];
    }

    return res;
}

template<class T>
vector<T> diff(const vector<T>& a, const vector<T>& b) {
    int n = min(a.size(), b.size());
    vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = a[i] - b[i];
    }

    return res;
}

template<class T>
vector<T> div(const vector<T>& a, T coef) {
    int n = a.size();
    vector<T> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] / coef;
    }
    return res;
}

template<class T>
vector<T> mul(const vector<T>& a, T coef) {
    int n = a.size();
    vector<T> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = a[i] * coef;
    }
    return res;
}

template<class T>
vector<T> operator+ (const vector<T>& a, const vector<T>& b) {
    sum(a, b);
}

template<class T>
vector<T> operator- (const vector<T>& a, const vector<T>& b) {
    diff(a, b);
}

template<class T>
vector<T> operator* (const vector<T>& a, const vector<T>& b) {
    mul(a, b);
}

template<class T>
vector<T> operator/ (const vector<T>& a, const vector<T>& b) {
    div(a, b);
}
