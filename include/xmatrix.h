#ifndef XMATRIX_H
#define XMATRIX_H

#include <vector>
#include <iostream>
#include "Number.h"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>



using INT = unsigned int;


using namespace std;

template<typename T>
class Sync_queue {
public:
    void put(const T& val);
    void put(T&& val);
    void get(T& val);
private:
    mutex mtx;
    condition_variable cond;
    list<T> q;
};

struct Message
{
    int a, b, c;
};


class FMatrix21
{
    public:
        int a, b, c;
        Field* F;
        int& operator()(int i, int j, int k) { return A[(i*b+j)*c+k]; };
        FMatrix21(Field* F0, int a0, int b0, int c0) : a(a0), b(b0), c(c0), F(F0) { A.resize(a*b*c);  };
    private:
        std::vector<int> A;
};


class FMatrix22
{
    public:
        int m, n, rk;
        int r1, r2, c1, c2;
        Field* F;
        FMatrix22(Field*, int, int, int, int);
        int& operator()(int i, int j) { return A[i*n+j]; };
        int& opLarge(INT i, INT j) { return A[(INT)i*n+(INT)j]; };
        int& operator()(int i, int j, int k, int l) { return A[(((INT)(i*r2+j)*c1+k)*(INT)c2+(INT)l)]; };
        int gauss();
        int gauss1();
        int gauss2();
    private:
        std::vector<int> A;
        Sync_queue<Message> mq;
        Sync_queue<Message> mq1;
        void loop();
};


class FMatrix
{
    public:
        int m, n, rk;
        Field* F;
        std::vector<int> free;
        std::vector<int> col_basis;
        /*** free[i]=1 indicates that the i-th component can be chosen freely ***/
        /*** col_basis transfers the information from free to the nullspace ***/
        FMatrix(Field*, int, int);
        int& operator()(int i, int j) { return A[i*n+j]; };
        FMatrix gauss_jordan();
        FMatrix nullspace();
        FMatrix transpose();
        FMatrix submatrix(int, int, int, int);
        void print();
    private:
        std::vector<int> A;
};


#endif // XMATRIX_H
