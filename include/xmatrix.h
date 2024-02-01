#ifndef XMATRIX_H
#define XMATRIX_H

#include <vector>
#include <iostream>
#include "Number.h"


class FMatrix
{
    public:
        int m, n, rk;
        int a, b, c;
        int r1, r2, c1, c2;
        Field* F;
        std::vector<int> free;
        std::vector<int> col_basis;
        /*** free[i]=1 indicates that the i-th component can be chosen freely ***/
        /*** col_basis transfers the information from free to the nullspace ***/
        FMatrix(Field*, int, int);
        FMatrix(Field*, int, int, int);
        FMatrix(Field*, int, int, int, int);
        FMatrix(Field*, int, FMatrix, int);

        int& operator()(int i, int j) { return A[i*n+j]; };
        int& operator()(int i, int j, int k) { return operator()(i*b+j,k); };
        int& operator()(int i, int j, int k, int l) { return operator()(i*r2+j,k*c2+l); };
        FMatrix GaussJordan();
        FMatrix Gauss();
        FMatrix Nullspace();
        FMatrix Transpose();
        FMatrix Submatrix(int, int, int, int);
        void prout();
    private:
        std::vector<int> A;
};


#endif // XMATRIX_H
