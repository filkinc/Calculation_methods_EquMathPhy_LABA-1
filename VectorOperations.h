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
    return sum(a, b);
}

template<class T>
vector<T> operator- (const vector<T>& a, const vector<T>& b) {
    return diff(a, b);
}

template<class T>
vector<T> operator- (const vector<T>& a) {
    int n = a.size();
    vector<T> res(n);

    for (int i = 0; i < n; ++i) {
        res[i] = -a[i];
    }

    return res;
}

template<class T>
vector<T> operator* (const vector<T>& v, T coef) {
    return mul(v, coef);
}

template<class T>
vector<T> operator* (T coef, const vector<T>& v) {
    return mul(v, coef);
}

template<class T>
vector<T> operator/ (const vector<T>& v, T coef) {
    return div(v, coef);
}

template<class T>
vector<T> operator+= (vector<T>& a, const vector<T>& b) {
    return a = sum(a, b);
}

template<class T>
vector<T> operator-= (vector<T>& a, const vector<T>& b) {
    return a = diff(a, b);
}

template<class T>
vector<T> operator*= (vector<T>& v, T coef) {
    return v = mul(v, coef);
}

template<class T>
vector<T> operator/= (vector<T>& v, T coef) {
    return v = div(v, coef);
}

template<typename T, typename Stream>
vector<T> readVector(Stream& stream, int size) {
    vector<T> readed(size);
    for (int i = 0; i < size; ++i) {
        stream >> readed[i];
    }

    return readed;
}
