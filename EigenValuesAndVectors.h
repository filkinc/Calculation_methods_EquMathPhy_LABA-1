#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "QuadMatrix.h"
#include "Matrix.h"
#include "Norma.h"
#include "LinSolveAlgs.h"

#define MAX_ITER 10

template<class T>
struct EigenValueSearchResults {
	vector<T> eigenValues;
	vector<vector<T>> eigenVector;
	int iterationCount = 0;
};

template <class T>
EigenValueSearchResults<T> qrSimpleSearchEigenvalue(const QuadMatrix<T>& A, size_t n, T eps) {
	auto [curQ, curR] = qrDecomposition(A);
	size_t iter = 0;
	QuadMatrix<T> nextA = A;
	if (n == 1) {
		EigenValueSearchResults<T> res;
		res.eigenValues = { nextA(0, 0) };

		res.iterationCount = iter;

		return res;
	}
		while (true) {
			nextA = curR * curQ;

			//Временно коментим
			if (iter < 3) {
				cout << endl /*<< "Вращение №: " << iter*/;
				for (int i = 0; i < n; i++) {
					cout << endl;
					for (int j = 0; j < n; j++) {
						cout << nextA(i, j) << " ";
					}
				}
				cout << endl;
			}


			auto [curQ1, curR1] = qrDecomposition(nextA);
			curQ = curQ1;
			curR = curR1;

			++iter;

			if (fabs(nextA(n - 1,n - 2)) <= eps) {
				break;
			}
		}
		
		T curEigen = nextA(n - 1, n - 1);

		nextA.resize(n - 1);

		auto res = qrSimpleSearchEigenvalue(nextA, n - 1, eps);
		res.eigenValues.push_back(curEigen);
		
		res.iterationCount += iter;
		return res;
}

template <class T>
EigenValueSearchResults<T> qrPlusShiftSearchEigenvalue(const QuadMatrix<T>& A, size_t n, T eps) {
	size_t iter = 0;
	//A.print();

	if (n == 1) {
		EigenValueSearchResults<T> res;
		res.eigenValues = { A(0, 0) };

		res.iterationCount = iter;

		return res;
	}
	
	T sigma = A(n - 1, n - 1);

	//Создание матрицы Е
	QuadMatrix<T> E(n);
	for (int i = 0; i < n; ++i) {
		E(i, i) = 1;
	}

	QuadMatrix<T> nextA = A - sigma * E;

	auto [curQ, curR] = qrDecomposition(nextA);
	

	while (true) {

		nextA = (curR * curQ) + (sigma * E);

		if (fabs(nextA(n - 1, n - 2)) <= eps) {
			break;
		}

		//Временно коментим
		if (iter < 3) {
			cout << endl /*<< "Вращение №: " << iter*/;
			for (int i = 0; i < n; i++) {
				cout << endl;
				for (int j = 0; j < n; j++) {
					cout << nextA(i, j) << " ";
				}
			}
			cout << endl;
		}

		sigma = nextA(n - 1, n - 1);
		nextA = nextA - (sigma * E);
		auto [curQ1, curR1] = qrDecomposition(nextA);
		curQ = curQ1;
		curR = curR1;
		
		++iter;
	}

	T curEigen = nextA(n - 1, n - 1);

	nextA.resize(n - 1);

	auto res = qrPlusShiftSearchEigenvalue(nextA, n - 1, eps);
	res.eigenValues.push_back(curEigen);
	res.iterationCount += iter;
	return res;
}

template <class T>
QuadMatrix<T> HesenbergMatrix(const QuadMatrix<T>& A, int k, int l) {
	size_t n = A.order();

	QuadMatrix<T> H(n);

	T alpha = (A(k - 1, k - 2)) / sqrt(A(k - 1, k - 2) * A(k - 1, k - 2) + A(l - 1, k - 2) * A(l - 1, k - 2));
	T beta = (A(l - 1, k - 2)) / sqrt(A(k - 1, k - 2) * A(k - 1, k - 2) + A(l - 1, k - 2) * A(l - 1, k - 2));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {	
			if (i == k - 1) {
				H(i, j) = alpha * A(k - 1, j) + beta * A(l - 1, j);
			}
			else if (i == l - 1) {
				H(i, j) = -beta * A(k - 1, j) + alpha * A(l - 1, j);
			}
			else {
				H(i, j) = A(i, j);
			}
		}
	}

	QuadMatrix<T> rH(n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (j == k - 1) {
				rH(i, j) = alpha * H(i, k - 1) + beta * H(i, l - 1);
			}
			else if (j == l - 1){
				rH(i, j) = -beta * H(i, k - 1) + alpha * H(i, l - 1);
			} 
			else {
				rH(i, j) = H(i, j);
			}
		}
	}

	/*cout << endl << "Коэффициенты:" << endl;
	cout << "alpha: " << alpha << endl;
	cout << "beta: " << beta << endl; 

	for (int i = 0; i < n; ++i) {
		cout << endl;
		for (int j = 0; j < n; ++j) {
			cout << rH(i, j) << " ";
		}
	}
	cout << endl;*/
	return rH;
}

template <class T>
EigenValueSearchResults<T> qrHesenbergSearchEigenvalue(const QuadMatrix<T>& A, int n, T eps) {
	QuadMatrix<T> H = A;
	for (int l = 3; l <= n; ++l) {
		for (int k = 2; k <= l - 1; ++k) {
			H = HesenbergMatrix(H, k, l);
		}
	}
	return qrSimpleSearchEigenvalue(H, n, eps);
}

template <class T>
EigenValueSearchResults<T> qrPlusShiftHesenbergSearchEigenvalue(const QuadMatrix<T>& A, int n, T eps) {
	QuadMatrix<T> H = A;
	for (int l = 3; l <= n; ++l) {
		for (int k = 2; k <= l - 1; ++k) {
			H = HesenbergMatrix(H, k, l);
		}
	}
	return qrPlusShiftSearchEigenvalue(H, n, eps);
}

template <class T>
EigenValueSearchResults<T> eigenVectorReverseIteration(const QuadMatrix<T>& A, vector<T> lambda, T eps) {
	int n = A.order();
	int k = 0;
	vector<T> curVecEigen(n);
	vector<T> nextVecEigen(n);
	curVecEigen[0] = 1;

	EigenValueSearchResults<T> res;

	QuadMatrix<T> E(n);
	for (int i = 0; i < n; ++i) {
		E(i, i) = 1;
	}

	for (int i = 0; i < n; ++i) {
		//int iter = 0;
		int iter = 0;
		QuadMatrix<T> B = (A - lambda[i] * E);
		auto Binv = B.inv();
		while (true) {
			nextVecEigen = Binv * curVecEigen;
			nextVecEigen = div(nextVecEigen, norm_2(nextVecEigen));
			if ( fabs(fabs(compose(nextVecEigen, curVecEigen)) - 1) <= eps ) {
				break;
			}
			/*printVector(nextVecEigen);
			cout << endl;*/
			curVecEigen = nextVecEigen;
			++iter;
		}
		res.iterationCount += iter;

		res.eigenVector.push_back(curVecEigen);
		//res.eigenValues += iter;
	}
	return res;
}

template<class T>
EigenValueSearchResults<T> eigenVectorReley(const QuadMatrix<T>& A, const vector<T>& vecEigen, T eps) {
	int n = A.order();
	int k = 0;
	vector<T> curVecEigen(n);
	vector<T> nextVecEigen(n);
	curVecEigen[0] = 1;

	vector<T> lambda = vecEigen;
	printVector(lambda);
	EigenValueSearchResults<T> res;

	QuadMatrix<T> E(n);
	for (int i = 0; i < n; ++i) {
		E(i, i) = 1;
	}

	for (int i = 0; i < n; ++i) {
		//int iter = 0;
		int iter = 0;
		QuadMatrix<T> B = (A - lambda[i] * E);
		auto Binv = B.inv();
		cout << "labda change " << i << endl;
		while (true) {
			cout << lambda[i] << endl;
			nextVecEigen = Binv * curVecEigen;
			nextVecEigen = div(nextVecEigen, norm_2(nextVecEigen));
			if (fabs(fabs(compose(nextVecEigen, curVecEigen)) - 1) <= eps) {
				break;
			}
			lambda[i] = compose(nextVecEigen, A.transposed() * nextVecEigen) / compose(nextVecEigen, nextVecEigen);
			
			curVecEigen = nextVecEigen;
			++iter;
		}
		res.iterationCount += iter;

		res.eigenVector.push_back(curVecEigen);
		//res.eigenValues += iter;
	}
	return res;
}

template<class T> 
vector<T> checkEigenVector(const QuadMatrix<T>& A, const vector<vector<T>>& vecVectorEigen, const vector<T>& lambda) {
	int n = vecVectorEigen.size();
	vector<T> chekEigenValue;
	for (int i = 0; i < n; ++i) {
		chekEigenValue.push_back( norm_2( diff( A * vecVectorEigen[i], coefcompose(lambda[i], vecVectorEigen[i]) ) ) );
	}
	return chekEigenValue;
}