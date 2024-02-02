#include "xmatrix.h"
#include "Number.h"


FMatrix::FMatrix(Field* F0, int rows, int cols)
{
    m = rows;
    n = cols;
    F = F0;
    A.resize(rows*cols);
    free.resize(n);
};

FMatrix::FMatrix(Field* F0, int a1, int b1, int c1)
{
    a = a1; b = b1; c = c1;
    m = a1*b1; n = c1;
    F = F0;
    A.resize(m*n);
    free.resize(n);
};

FMatrix::FMatrix(Field* F0, int rows1, int rows2, int cols1, int cols2)
{
    r1 = rows1; r2 = rows2; c1 = cols1; c2 = cols2;
    m = rows1*rows2;
    n = cols1*cols2;
    F = F0;
    A.resize(m*n);
    free.resize(n);
};

FMatrix::FMatrix(Field* F0, int rows, FMatrix M, int start)
{
    m = rows;
    n = M.n;
    F = F0;
    A.resize(m*n);
    for (int i=0;i<m*n;i++)
        A[i] = M.A[start+i];
    free.resize(n);
};


FMatrix FMatrix::Submatrix(int row0, int row1, int col0, int col1)
{
    FMatrix B(F, row1-row0, col1-col0);

    for (int i=0; i<B.m; i++)
        for (int j=0; j<B.n; j++)
            B(i,j) = operator()(row0+i, col0+j);
    return B;
};


FMatrix FMatrix::GaussJordan()
{
    int i,j,pi,pj,x;

    pi=pj=0;
    for (i=0;i<n;i++) free[i]=1;

    while (pi<m && pj<n)
    {
        if ( operator()(pi,pj) )
        {
            free[pj]=0;  // x_pj is not a free variable
            x = F->inverse(operator()(pi,pj));

            for (j=pj;j<n;j++)
                operator()(pi,j) = F->product(x,operator()(pi,j));  // multiply row to make pivot =1 (change basis in image)
            for (i=0;i<m;i++) if (pi != i) for (j=n-1;j>=pj;j--)
                operator()(i,j) = F->sum(operator()(i,j), F->neg(F->product(operator()(i,pj), operator()(pi,j))));
                // subtract multiple of pivot row to set the full column to 0  (change basis in image)
            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i != m && operator()(i,pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i != m)   // swap rows pi and i
                for (j=pj;j<n;j++)
                {
                    x = operator()(i,j);
                    operator()(i,j) = operator()(pi,j);
                    operator()(pi,j) = x;
                }
            else pj++;
        }
    }
    rk = pi;

    return *this;
}


FMatrix FMatrix::Nullspace()
{
    int i,j;

    FMatrix N(F, n, n-rk);

    int col = 0;

    N.col_basis.resize(n-rk);

    for (i=0;i<n;i++) if (free[i])
    {
        int equ = 0;
        for (j=0;j<n;j++)
            N(j,col) = (free[j]) ? 0 : F->neg(operator()(equ++,i));
        N(i,col) = 1;
        N.col_basis[col] = i;
        col++;
    }

    return N;
}


FMatrix FMatrix::Gauss()
{
    int i,j,pi,pj,x;

    pi=pj=0;
    for (i=0;i<n;i++) free[i]=1;

    while (pi<m && pj<n)
    {
        if ( operator()(pi,pj) )
        {
            free[pj]=0;  // x_pj is not a free variable
            x = F->inverse(operator()(pi,pj));

            for (j=pj;j<n;j++)
                operator()(pi,j) = F->product(x,operator()(pi,j));  // multiply row to make pivot =1 (change basis in image)
            for (i=pi+1;i<m;i++) if (operator()(i,pj) != 0) {
                int* pt = &(F->exp[F->log[F->neg(operator()(i,pj))]]);
                for (j=n-1;j>=pj;j--) {
       //             operator()(i,j) = F->sum(operator()(i,j), F->neg(F->product(operator()(i,pj), operator()(pi,j))));
                    if (operator()(pi,j) != 0) {
                        int C = F->rebase[pt[F->log[operator()(pi,j)]]];
                        operator()(i,j) = F->unbase[F->rebase[operator()(i,j)]+C];
                    };
                };
            };
            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i != m && operator()(i,pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i != m)   // swap rows pi and i
                for (j=pj;j<n;j++)
                {
                    x = operator()(i,j);
                    operator()(i,j) = operator()(pi,j);
                    operator()(pi,j) = x;
                }
            else pj++;
        }
    }
    rk = pi;

    return *this;
}


FMatrix FMatrix::Transpose()
{
    int i, j;
    FMatrix B(F,n,m);

    for (i=0;i<n;i++)
        for (j=0;j<m;j++)
            B(i,j) = operator()(j,i);

    return B;
}


void FMatrix::prout()
{
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++) std::cout << operator()(i,j) << " ";
        std::cout << "\n";
    };
    std::cout << "\n";
}



