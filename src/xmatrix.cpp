#include "xmatrix.h"
#include "Number.h"


template<typename T>
void Sync_queue<T>::put(const T& val)
{
    lock_guard<mutex> lck(mtx);
    q.push_back(val);
    cond.notify_one();
}

template<typename T>
void Sync_queue<T>::get(T& val)
{
    unique_lock<mutex> lck(mtx);
    cond.wait(lck,[this]{ return !q.empty(); });
    val = q.front();
    q.pop_front();
}

/*** FMatrix   ************************************/

FMatrix FMatrix::submatrix(int row0, int row1, int col0, int col1) const
{
    FMatrix B(F, row1-row0, col1-col0);

    for (int i=0; i<B.m; i++)
        for (int j=0; j<B.n; j++)
            B(i,j) = operator()(row0+i, col0+j);
    return B;
};


FMatrix& FMatrix::gauss_jordan()
{
    int i,j,pi,pj;

    pi=pj=0;
    for (i=0;i<n;i++) free[i]=1;

    while (pi<m && pj<n)
    {
        if ( const auto v = operator()(pi,pj) )
        {
            free[pj]=0;  // x_pj is not a free variable
            const auto x = F->inverse( v );

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
                    const auto x = operator()(i,j);
                    operator()(i,j) = operator()(pi,j);
                    operator()(pi,j) = x;
                }
            else pj++;
        }
    }
    rk = pi;

    return *this;
}


FMatrix FMatrix::nullspace() const
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


FMatrix FMatrix::transpose() const
{
    FMatrix B(F,n,m);

    for (int i=0;i<n;i++)
        for (int j=0;j<m;j++)
            B(i,j) = operator()(j,i);

    return B;
}


void FMatrix::print() const
{
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++) std::cout << F->decode[operator()(i,j)] << " ";
        std::cout << "\n";
    };
    std::cout << "\n";
}

/*** FMatrix22   ********************************/

int Matrix22::gauss()
{
    int i,j,pi,pj;
    SHORT x;

    pi=pj=0;
//    for (i=0;i<n;i++) free[i]=1;

    while (pi<m && pj<n)
    {
        if ( opLarge(pi,pj) )
        {
//            free[pj]=0;  // x_pj is not a free variable
            x = F->inverse(opLarge(pi,pj));

            for (j=pj;j<n;j++)
                opLarge(pi,j) = F->product(x,opLarge(pi,j));  // multiply row to make pivot =1 (change basis in image)
            for (i=pi+1;i<m;i++) if (opLarge(i,pj) != 0) {
                const SHORT* pt = &(F->exp[F->log[F->neg(opLarge(i,pj))]]);
                for (j=n-1;j>=pj;j--) {
                    if (opLarge(pi,j) != 0) {
                        SHORT C = pt[F->log[opLarge(pi,j)]];
                        opLarge(i,j) = F->normalize[opLarge(i,j)+C];
                    };
                };
            };

            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i != m && opLarge(i,pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i != m)   // swap rows pi and i
                for (j=pj;j<n;j++)
                {
                    x = opLarge(i,j);
                    opLarge(i,j) = opLarge(pi,j);
                    opLarge(pi,j) = x;
                }
            else pj++;
        }
    }
    rk = pi;

    return rk;
}


int Matrix22::gauss1()
{
    int i,j,pi,pj;
    int x;

    pi=pj=0;
    std::vector<int> row(m);
    for (int i=0; i<m; i++) row[i] = i;

    while (pi<m && pj<n)
    {
        if ( opLarge(row[pi],pj) )
        {
            for (i=pi+1;i<m;i++) if (opLarge(row[i],pj) != 0) {
                const SHORT* pt = &(F->exp[0]) + F->log[F->neg(opLarge(row[i],pj))] + F->log[F->inverse(opLarge(row[pi],pj))];
                for (j=n-1;j>=pj;j--) {
                    if (opLarge(row[pi],j) != 0) {
                        SHORT C = pt[F->log[opLarge(row[pi],j)]];
                        opLarge(row[i],j) = F->normalize[opLarge(row[i],j)+C];
                    };
                };
            };

            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i < m && opLarge(row[i],pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i < m)   // swap rows pi and i
                {
                    x = row[i];
                    row[i] = row[pi];
                    row[pi] = x;
                }
            else pj++;
        }
    }
    rk = pi;

    return rk;
}


int Matrix22::gauss2()
{
    int i,pi,pj;
    int x;

    pi=pj=0;
    std::vector<int> row(m);
    for (int i=0; i<m; i++) row[i] = i;
    int tasks = thread::hardware_concurrency() - 1;
/**************************************************************
//  spawn threads
***************************************************************/
    thread t[tasks];
    for (int i=0; i<tasks; i++) t[i] = thread{&Matrix22::loop, this};

    while (pi<m && pj<n)
    {
        if ( opLarge(row[pi],pj) != 0)
        {
            int jobs = 0;
            Message ms;
            ms.a = row[pi]; ms.b = pj;
            for (i=pi+1;i<m;i++) if (opLarge(row[i],pj) != 0) {
/**************************************************************
//  add (pi,pj,i) to task queue
***************************************************************/
                ms.c = row[i];
                mq.put(ms);
                jobs++;
            };
/**************************************************************
//  wait for completion (count down jobs)
***************************************************************/
            while (jobs > 0){
                Message ms;
                mq1.get(ms);
                jobs--;
            };

            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i < m && opLarge(row[i],pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i < m) {   // swap rows pi and i
                x = row[i];
                row[i] = row[pi];
                row[pi] = x;
            }
            else pj++;
        }
    }
/**************************************************************
//  send stop messages to queue
***************************************************************/
    Message ms;
    ms.a = ms.b = ms.c = -1;
    for (int i=0;i<tasks;i++) mq.put(ms);  // stop messages
/**************************************************************
//  join threads
***************************************************************/
    for (int i=0; i<tasks; i++) t[i].join();

    rk = pi;
    return rk;
}


void Matrix22::loop()
{
    while( true ) {
        Message ms;
        mq.get(ms);
        if (ms.a == -1) return;  // received termination message
        int pi = ms.a;
        int pj = ms.b;
        int i = ms.c;

        const INT r0 = (INT) pi * (INT) n;
        const INT r1 = (INT )i * (INT) n;

        const SHORT* pt = &(F->exp[0]) + F->log[F->neg(opLarge(i,pj))] + F->log[F->inverse(opLarge(pi,pj))];
        for (int j=n-1;j>=pj;j--) {
            if ( const auto v1 = A[ r0 + j ] )
            {
                const SHORT C = pt[ F->log[v1] ];
                auto& v2 = A[ r1 + j ];
                v2 = F->normalize[ v2 + C ];
            };
        };

        mq1.put(ms);  // report completion of task
    }
}
