#include "xmatrix.h"
#include "Number.h"

#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>

/*** FMatrix   ************************************/

FMatrix FMatrix::submatrix(int row0, int row1, int col0, int col1) const
{
    FMatrix B(F, row1-row0, col1-col0);

    for (int i=0; i<B.m_rows; i++)
        for (int j=0; j<B.m_cols; j++)
            B(i,j) = operator()(row0+i, col0+j);
    return B;
};


FMatrix& FMatrix::gauss_jordan()
{
    int i,j,pi,pj;

    pi=pj=0;

    for (i=0;i<m_cols;i++)
        m_free[i]=true;

    while (pi<m_rows && pj<m_cols)
    {
        if ( const auto v = operator()(pi,pj) )
        {
            m_free[pj]=false;  // x_pj is not a free variable
            const auto x = F->inverse( v );

            for (j=pj;j<m_cols;j++)
                operator()(pi,j) = F->product(x,operator()(pi,j));  // multiply row to make pivot =1 (change basis in image)
            for (i=0;i<m_rows;i++) if (pi != i) for (j=m_cols-1;j>=pj;j--)
                operator()(i,j) = F->sum(operator()(i,j), F->neg(F->product(operator()(i,pj), operator()(pi,j))));
                // subtract multiple of pivot row to set the full column to 0  (change basis in image)
            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i != m_rows && operator()(i,pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i != m_rows)   // swap rows pi and i
                for (j=pj;j<m_cols;j++)
                {
                    const auto x = operator()(i,j);
                    operator()(i,j) = operator()(pi,j);
                    operator()(pi,j) = x;
                }
            else pj++;
        }
    }
    m_rk = pi;

    return *this;
}


FMatrix FMatrix::nullspace() const
{
    int i,j;

    FMatrix N(F, m_cols, m_cols-m_rk);

    int col = 0;

    N.m_col_basis.resize(m_cols-m_rk);

    for (i=0;i<m_cols;i++)
    {
        if (m_free[i])
        {
            int equ = 0;
            for (j=0;j<m_cols;j++)
                N(j,col) = (m_free[j]) ? 0 : F->neg(operator()(equ++,i));
            N(i,col) = 1;
            N.m_col_basis[col] = i;
            col++;
        }
    }

    return N;
}


FMatrix FMatrix::transpose() const
{
    FMatrix B(F,m_cols,m_rows);

    for (int i=0;i<m_cols;i++)
        for (int j=0;j<m_rows;j++)
            B(i,j) = operator()(j,i);

    return B;
}


void FMatrix::print() const
{
    for (int i=0; i<m_rows; i++)
    {
        for (int j=0; j<m_cols; j++) std::cout << F->decoded( operator()(i,j) ) << " ";
        std::cout << "\n";
    };
    std::cout << "\n";
}

/*** FMatrix22   ********************************/

class Matrix22::MessageQueue
{
    public:
        struct Message { int a, b, c; };

        inline void put(const Message& val)
        {
            using namespace std;

            lock_guard<mutex> lck(mtx);
            q.push_back(val);
            cond.notify_one();
        }

        inline void get(Message& val)
        {   
            using namespace std;

            unique_lock<mutex> lck(mtx);
            cond.wait(lck,[this]{ return !q.empty(); });
            val = q.front();
            q.pop_front();
        }

    private:
        std::mutex mtx;
        std::condition_variable cond;
        std::list<Message> q;
};
    
Matrix22::Matrix22(const Field* F, int rows1, int rows2, int cols1, int cols2)
    : F(F)
    , m_rows(rows1*rows2)
    , m_cols(cols1*cols2)
    , m_r1(rows1)
    , m_r2(rows2)
    , m_c1(cols1)
    , m_c2(cols2)
    , m_mq( new MessageQueue() )
    , m_mq1( new MessageQueue() )
{
    m_values.resize((INT)m_rows * (INT)m_cols);
}

Matrix22::~Matrix22() 
{
}

int Matrix22::gauss2()
{
    int i,pi,pj;
    int x;

    pi=pj=0;
    std::vector<int> row(m_rows);
    for (int i=0; i<m_rows; i++) row[i] = i;
    int tasks = std::thread::hardware_concurrency() - 1;
/**************************************************************
//  spawn threads
***************************************************************/
    std::thread t[tasks];
    for (int i=0; i<tasks; i++) t[i] = std::thread{&Matrix22::loop, this};

    while (pi<m_rows && pj<m_cols)
    {
        if ( opLarge(row[pi],pj) != 0)
        {
            int jobs = 0;
            MessageQueue::Message ms;
            ms.a = row[pi]; ms.b = pj;
            for (i=pi+1;i<m_rows;i++) if (opLarge(row[i],pj) != 0) {
/**************************************************************
//  add (pi,pj,i) to task queue
***************************************************************/
                ms.c = row[i];
                m_mq->put(ms);
                jobs++;
            };
/**************************************************************
//  wait for completion (count down jobs)
***************************************************************/
            while (jobs > 0){
                MessageQueue::Message ms;
                m_mq1->get(ms);
                jobs--;
            };

            pi++;
            pj++;
        }
        else
        {
            i = pi+1;
            while (i < m_rows && opLarge(row[i],pj) == 0) i++; // if pivot = 0, search for nonzero entry in same column

            if (i < m_rows) {   // swap rows pi and i
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
    MessageQueue::Message ms;
    ms.a = ms.b = ms.c = -1;
    for (int i=0;i<tasks;i++) m_mq->put(ms);  // stop messages
/**************************************************************
//  join threads
***************************************************************/
    for (int i=0; i<tasks; i++) t[i].join();

    m_rk = pi;
    return m_rk;
}


void Matrix22::loop()
{
    while( true ) {
        MessageQueue::Message ms;
        m_mq->get(ms);
        if (ms.a == -1) return;  // received termination message

        const int pi = ms.a;
        const int pj = ms.b;
        const int i = ms.c;

        const INT r0 = (INT) pi * (INT) m_cols;
        const INT r1 = (INT )i * (INT) m_cols;

        const auto v = opLarge(i,pj);

        const auto offset = F->log( F->neg( v ) ) + F->log( F->inverse( v ) );

        for (int j = m_cols-1; j >= pj; j--)
        {
            if ( const auto v1 = m_values[ r0 + j ] )
            {
                auto& v2 = m_values[ r1 + j ];
                v2 = F->normalized( v2 + F->exp( offset + F->log( v1 ) ) );
            };
        };

        m_mq1->put(ms);  // report completion of task
    }
}
