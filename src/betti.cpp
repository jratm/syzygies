#include "betti.h"
#include "xmatrix.h"
#include <iostream>


static int Fgcd(const Field*, const std::vector <int>&, const std::vector <int>&, int); // return gcd of p1 and p2
static Matrix21 reduced(const Matrix21&);
static void print_line(int,int,int,int,int);


static long choose(int a, const int b)
{
    if (a<b) return 0;
    long ret = 1;
    for (int i=1; i<=b; i++)
    {
        ret = ret * a / i;
        a--;
    };
    return ret;
}


static Matrix21 multTable( const LineBundle& L, const LineBundle& B)
{
    LineBundle LB = LBmult(L,B);
    FMatrix ML = L.C->sections(L);
    FMatrix MB = L.C->sections(B);
    FMatrix MBL = L.C->sections(LB);

// precalculate the multiplication table  H^0(L) x H^0(B) --> H^0(L tensor B)
    Matrix21 M(L.C->F, ML.columnCount(), MB.columnCount(), MBL.columnCount());    //  if L=B=K, then M has size  ( g x g x (3g-3) )

    std::vector<int> Lproduct(LB.degree+1);

    for (auto i=0; i<ML.columnCount(); i++)
        for (auto j=0; j<MB.columnCount(); j++)
        {
            for (auto i=0; i<LB.degree+1; i++) Lproduct[i] = 0;
            for (auto i1=0; i1<=L.degree; i1++)
                for (auto i2=0; i2<=B.degree; i2++)
                    Lproduct[i1+i2] = L.C->F->sum(Lproduct[i1+i2], L.C->F->product(ML(i1,i), MB(i2,j)));
            for (auto i3=0; i3<MBL.columnCount(); i3++)
                M(i, j, i3) = Lproduct[MBL.columnBasis(i3) ];
        };
// we want to reduce the matrix by two hyperplane sections
// by construction, the lower part of ML is the unit matrix
// This necessitates to reduce by the leftmost column (otherwise, the point at infinity is ...
// we also select the column next to it and check (via Fgcd) whether there is a point of C in their intersection

    const int col2 = ML.columnCount();
    const int col1 = col2 - 1;

    std::vector<int> v1(L.degree + 1), v2(L.degree + 1);

    for (int i=0; i<ML.rowCount(); i++)
    {
        v2[i] = ML(i, col2);
        v1[i] = ML(i, col1);
    }

    if ( Fgcd(L.C->F, v1, v2, L.degree) != 0)
        std::cout << "success " << col1 << " " << col2 << "\n";
    else
        std::cout << "failure\n";

    const Matrix21 M0 = reduced(M);
    const Matrix21 M1 = reduced(M0);

    return M1;
};


static Matrix21 reduced(const Matrix21& M)
{
    FMatrix M1(M.F, M.c, M.a * M.b);
    const int g = M.a;
    std::vector<int> exch(M.a);

    for (int i=0; i<M.a; i++) exch[i] = M.a - 1 - i;

    for (int i=0; i<M.a; i++)
        for (int j=0; j<M.b; j++)
            for (int k=0; k<M.c; k++)
//                M1(k, i * M.b + j) = M(i, j, k);
                M1(k, i * M.b + j) = M(exch[i], exch[j], k);

    M1.gauss_jordan();
// Gaussian elimination changes the basis of the image; here
// the left-most columns are tested for inclusion in a basis;
// since we know that the g left-most columns are independent
// and mapped to zero in the quotient
//  ==> the bottom-right submatrix of M1 represents a basis for the quotient

    Matrix21 M2(M.F, M.a-1, M.b-1, M.c-g);

    for (int i=0; i<M2.a; i++)
        for (int j=0; j<M2.b; j++)
            for (int k=0; k<M2.c; k++)
//                M2(i, j, k) = M1(g + k, (i+1)*M.b + j + 1);
                M2(i, j, k) = M1(g + k, exch[i] * M.b + exch[j]);

    return M2;
}


BettiTable::BettiTable(const Curve* C0)
    : C( C0 )
{
// calculate Betti table for the canonical bundle
    const int g = C->genus();
    const LineBundle L = C->canonical();
    initialize(L);

    const Matrix21 M1 = multTable(L,L);

    m_betti[0] = m_betti[4*g-5] = 1;

    for (auto p=1; 2*p<g; p++){
        int r01 = choose(M1.a, p+1);
        int m1 = choose(M1.a, p) * M1.a;
        int m2 = choose(M1.a, p-1) * M1.c;
        int r = m1 - r01 - run(p, M1);
        m_betti[(g-1)+p] = m_betti[(g-1)*2+(g-1-p-1)] = r;
        m_betti[(g-1)*2+(p-1)] = m_betti[(g-1)*2-p] = r + m_chi[p+1];
        print_line(m1, m2, p, 1, r);
    };
};


BettiTable::BettiTable(const Curve* C0, const LineBundle& L)
    : C( C0 )
{
    initialize(L);

    const Matrix21 M1 = multTable(L,L);
    m_betti[0] = 1;

    for (auto p=1; p<m_h0-2; p++){
        const int r01 = choose(M1.a, p+1);
        const int m1 = choose(M1.a, p) * M1.a;
        const int m2 = choose(M1.a, p-1) * M1.c;
        const int r = m1 - r01 - run(p, M1);
        m_betti[(m_h0-1)+p] = r;
        m_betti[(m_h0-1)*2+(p-1)] = r + m_chi[p+1];
        print_line(m1, m2, p, 1, r);
    };
};


void BettiTable::initialize(const LineBundle& L)
{
    const FMatrix ML = L.C->sections(L);
    m_h0 = ML.columnCount(); // h^0 L (= g)

    const int g = C->genus();
    const int deg = L.degree;

    m_betti.resize(4 * (m_h0-1));
    m_dim.resize(m_h0 * m_h0);
    m_chi.resize(m_h0);
    m_coimage.resize(4*(m_h0-1));

    // calculate Euler characteristic of the Koszul complex
    int sign = 1;
    for (auto i=0; i<m_h0; i++){
        int N = (i == 1) ? m_h0 : i*deg +1-g;
        if (i == 0) N = 1;
        for (auto j=0; j<m_h0; j++){
            m_dim[i*m_h0 + j] = N * choose(m_h0, j);  // only i=0..3 relevant
            if (i+j < m_h0) m_chi[i+j] += sign * N * choose(m_h0, j);
        };
        sign = -sign;
    };

    return;
}


void BettiTable::print() const
{
//    int g = C->genus;

    std::cout << "\n     ";
    for (auto p=0; p<m_h0-1; p++){
        std::cout.width(5);
        std::cout << p << " ";
    };
    std::cout << "\n";
    std::cout << "     ";
    for (auto p=0; p<m_h0-1; p++){
        std::cout << "------";
    };
    std::cout << "\n";

    for (auto q=0; q<4; q++){
        std::cout.width(3);
        std::cout << q << " |";
        for (auto p=0; p<m_h0-1; p++){
            std::cout.width(5);
            std::cout << m_betti[q*(m_h0-1)+p] << " ";
        };
        std::cout << "\n";
    };
};


int run(const int t, const Matrix21& M)
{
/*
generate all combinations, using Algorithm T from Knuth, TAOCP Vol. 4B, p. 359
*/
// n = L.dim = M.a, t as function parameter
    std::vector<int> c(t+2);
// step_T1:
    for (int l=0; l<t; l++) c[l] = l;
    c[t] = M.a;
    c[t+1] = 0;

    int l = t-1;
    int x;

    Matrix22 MM(M.F, choose(M.a,t), M.b, choose(M.a,t-1), M.c);
    int srow = 0;

    while (true)
    {
// run work
        for (int k=0; k<t; k++)
        {
            // int multiplier = c[k];
            bool sign = (k % 2) == 0;
            int trow = 0;
            for (int k1=0; k1<t; k1++){
                if (k1 != k)
                    trow += (k1<k) ? choose(c[k1],k1+1) : choose(c[k1],k1);
            };
            // add to result matrix
            for (int k2=0; k2<M.b; k2++)
                for (int k3=0; k3<M.c; k3++)
                {
                    if (sign)
                        MM(srow, k2, trow, k3) = M(c[k], k2, k3);
                    else
                        MM(srow, k2, trow, k3) = M.F->neg(M(c[k], k2, k3));
                };
        };
        srow++;
// step_T2:
        if (l>=0) x = l+1;
        else {
// step_T3:
            if (c[0]+1 < c[1]){
                c[0]++;
                continue;
            };
            l = 0;
// step_T4:
            do {
                c[l] = l;
                l++;
                x = c[l]+1;
            }
            while (x == c[l+1]);
// step_T5:
            if (l>=t) break;
        };
// step_T6:
        c[l] = x;
        l--;
    };

    return MM.gauss2();
}


void print_line(const int m1, const int m2, const int p, const int q, const int r)
{
    using namespace std;

    cout << "Morphism = ";
    cout.width(5);
    cout << m1 << " x ";
    cout.width(5);
    cout << m2 << ":  ";
    cout << "K_(" << p << "," << q << ") = ";
    cout.width(5);
    cout << r << "\n";
};


int Fgcd(const Field* F, const std::vector <int>& p1, const std::vector <int>& p2, int deg) // return gcd of p1 and p2
{
    std::vector <int> a = p2;
    std::vector <int> b = p1;

    int lead_a = deg;
    int lead_b = deg;

    while (true) {

    // determine leading term of b
        for (; lead_b>=0; lead_b--)
            if (b[lead_b] != 0) break;
        if (lead_b == -1) { return 0; }; // failure: b is a common factor of poly and h
        if (lead_b == 0)  { return lead_a;  }; // success: poly and h have no common factor

    //  make lead coefficient 1
        for (int i=0; i<=lead_b; i++) b[i] = F->product(b[i],F->inverse(b[lead_b]));

    //  get remainder of division
        for (int i=lead_a; i>=lead_b; i--){
            const int factor = a[i];
            if (factor != 0)
                for (int j=0; j<=lead_b; j++)
                    a[i+j-lead_b] = F->sum(a[i+j-lead_b], F->neg( F->product(factor, b[j])));
        };
    // exchange a and b
        for (int i=0; i<=lead_b; i++)
            std::swap( a[i], b[i] );
        lead_a = lead_b;
        lead_b --;
    };
    return -1;
};

//******* unused so far ********************

int Koszul(const int p, const int q, const LineBundle& L, const LineBundle& B)
{
    LineBundle B1 = B;

    int qq = q;
    while (qq>0){
        B1 = LBmult(B1, L);
        qq--;
    };
    while (qq<0){
        B1 = LBmult(B1, LBinverse(L));
        qq++;
    };
    const LineBundle B0 = LBmult(B1, LBinverse(L));

    const Matrix21 M2 = multTable(L,B1);
    int r12 = run(p, M2);
    const Matrix21 M1 = multTable(L,B0);
    int r01 = run(p+1, M1);
//    int r01 = choose(M1.a, p+1);

    const int m1 = choose(M1.a, p) * M1.c;
    const int m2 = choose(M2.a, p-1) * M2.c;

    print_line(m1, m2, p, q, m1-r12-r01);

    return m1-r12-r01;
}


int Koszul(const int p, const int q, const LineBundle& L)
{
    return Koszul(p, q, L, LBmult(L,LBinverse(L)));
}

