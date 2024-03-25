#pragma once

#include <vector>
#include <iostream>
#include <cstdint>
#include "Number.h"

#include <memory>

using INT = int64_t;

class Matrix21
{
    public:
        inline Matrix21(const Field* F, int a, int b, int c)
            : F(F)
            , a(a)
            , b(b)
            , c(c)
        {
            m_values.resize(a*b*c);
        }

        inline int& operator()(int i, int j, int k)
        {
            return m_values[(i*b+j)*c+k];
        }

        inline int operator()(int i, int j, int k) const
        {
            return m_values[(i*b+j)*c+k];
        };

        const Field* F;
        const int a, b, c;

    private:
        std::vector<int> m_values;
};


class Matrix22
{
    public:
        Matrix22(const Field*, int rows1, int rows2, int cols1, int cols2);
        ~Matrix22(); 

        inline SHORT& opLarge(INT i, INT j)
        {
            return m_values[(INT)i*m_cols+(INT)j];
        }

        inline SHORT& operator()(int i, int j, int k, int l)
        {
            return m_values[(((INT)(i*m_r2+j)*m_c1+k)*(INT)m_c2+(INT)l)];
        }

        int gauss2();

        const Field* F;

    private:
        Matrix22( const Matrix22& ) = delete;
        Matrix22( Matrix22 && ) = delete;
        Matrix22& operator=( const Matrix22 & ) = delete;
        Matrix22& operator=( Matrix22 && ) = delete;

        void loop();

        const int m_rows, m_cols;
        const int m_r1, m_r2, m_c1, m_c2;
        int m_rk;

        std::vector<SHORT> m_values;

        class MessageQueue;
        std::unique_ptr<MessageQueue> m_mq;
        std::unique_ptr<MessageQueue> m_mq1;
};


class FMatrix
{
    public:
        /*** free[i]=1 indicates that the i-th component can be chosen freely ***/
        /*** col_basis transfers the information from free to the nullspace ***/
//        FMatrix(Field*, int, int);
        FMatrix(const Field* F, int rows, int cols)
            : F(F)
            , m_rows(rows)
            , m_cols(cols)
        {
            A.resize(rows*cols);
            m_free.resize(cols);
        }

        int& operator()(int i, int j) { return A[i*m_cols+j]; }
        int operator()(int i, int j) const { return A[i*m_cols+j]; }

        FMatrix& gauss_jordan();
        FMatrix nullspace() const;
        FMatrix transpose() const;
        FMatrix submatrix(int, int, int, int) const;
        void print() const ;

        const Field* const F;

        inline int rowCount() const { return m_rows; }
        inline int columnCount() const { return m_cols; }
        inline int columnBasis( int index ) const { return m_col_basis[ index ]; }

    private:
        const int m_rows, m_cols;
        int m_rk = 0;
        std::vector<bool> m_free;
        std::vector<int> m_col_basis;
        std::vector<int> A;
};
