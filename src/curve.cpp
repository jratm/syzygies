#include "curve.h"
#include "xmatrix.h"
#include <iostream>


using namespace std;


long choose(int a, int b)
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


LineBundle LBmult(LineBundle L_0, LineBundle L_1)
{
    LineBundle res;
    res.degree = L_0.degree + L_1.degree;
    res.ratios.resize(L_0.C->genus);
    res.C = L_0.C;
    for (int i=0; i<L_0.C->genus; i++) res.ratios[i] = L_0.C->F->product(L_0.ratios[i], L_1.ratios[i]);
    return res;
};


LineBundle LBinverse(LineBundle L_0)
{
    LineBundle res;
    res.degree = -L_0.degree;
    res.ratios.resize(L_0.C->genus);
    res.C = L_0.C;
    for (int i=0; i<L_0.C->genus; i++) res.ratios[i] = L_0.C->F->inverse(L_0.ratios[i]);
    return res;
}


LineBundle Curve::canonical()
{
    int i,j,a,b;
    int deg = 2 * genus - 2;
    FMatrix A(F, 2 * genus - 1, 2 * genus);

    for (i=0; i<genus; i++)
    {
        a = b = 1;
        for (j=0; j<deg+1; j++)
        {
            A(j,2*i) = a;
            A(j,2*i+1) = b;
            a = F->product(a, nodes[i].p);
            b = F->product(b, nodes[i].q);
        }
    };

    FMatrix C = A.gauss_jordan().nullspace();

    LineBundle kan;
    kan.degree = deg;
    kan.ratios.resize(genus);
    kan.C = this;

    for (i=0; i<genus; i++)
    {
        int res = F->product(F->neg(C(2*i+1,0)), F->inverse(C(2*i,0)));
        kan.ratios[i] = res;
    }

    return kan;
}


LineBundle Curve::point(int p0)
{
    LineBundle Lp;
    int p1 = F->encode[p0];

    Lp.degree = 1;
    Lp.ratios.resize(genus);
    Lp.C = this;

    for (int i=0; i<genus; i++)
    {
        int res = F->product(F->sum(nodes[i].p, F->neg(p1)), F->inverse(F->sum(nodes[i].q, F->neg(p1))));
        Lp.ratios[i] = res;
    }

    return Lp;
}


FMatrix Curve::sections(LineBundle L)
{
    int i,j,k,l,a,b,c;

    int deg = L.degree;
    FMatrix A(F, genus, deg+1);

    for (i=0;i<genus;i++)
    {
        k = nodes[i].p;
        l = nodes[i].q;
        c = L.ratios[i];
        a = b = 1;
        for (j=0;j<deg+1;j++)
        {
            A(i,j) = F->sum(a, F->neg(F->product(b, c)));
            a = F->product(a, k);
            b = F->product(b, l);
        }
    }

    FMatrix C = A.gauss_jordan().nullspace();

    return C;
};


std::vector<node> Curve::sample_nodes(int g, int q)
{
    int n = 2*g;
    int N = q-2;
    std::vector<int> sample(n);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto rand = bind(uniform_real_distribution<double>{0,1},default_random_engine(seed));

    int m = 0;

    for (int t=0; t<q; t++){
        double u = rand();
        if ((N-t)*u >= n-m) continue;
        sample[m] = t+1;
        m++;
        if (m == n) break;
    }

    int j = n-1;
    while (j>0){
        double u = rand(); int k = floor(j*u) ;
        int x = sample[k]; sample[k] = sample[j]; sample[j] = x;
        j--;
    };

    std::vector<node> nodes(g);
    for (int i=0; i<g; i++) {
        nodes[i].p = sample[2*i];
        nodes[i].q = sample[2*i+1];
    };

    return nodes;
}


void Curve::print()
{
    std::cout << "Curve properties:\n";
    std::cout << "     genus = " << genus << "\n";
    std::cout << "     Nodes as follows: \n";
    for (int i=0; i<genus; i++){
        std::cout << "(" << F->decode[nodes[i].p] << "," << F->decode[nodes[i].q] << ")";
        if (i != genus-1) std::cout << ", ";
    };
    std::cout << "\n\n";

    return;
};


Matrix21 multTable(LineBundle L, LineBundle B)
{
    LineBundle LB = LBmult(L,B);
    FMatrix ML = L.C->sections(L);
    FMatrix MB = L.C->sections(B);
    FMatrix MBL = L.C->sections(LB);

// precalculate the multiplication table  H^0(L) x H^0(B) --> H^0(L tensor B)
    Matrix21 M(L.C->F, ML.n, MB.n, MBL.n);    //  if L=B=K, then M has size  ( g x (3g-3) x (5g-5) )

    std::vector<int> Lproduct(LB.degree+1);

    for (int i=0; i<ML.n; i++)
        for (int j=0; j<MB.n; j++)
        {
            for (auto i=0; i<LB.degree+1; i++) Lproduct[i] = 0;
            for (int i1=0; i1<=L.degree; i1++)
                for (int i2=0; i2<=B.degree; i2++)
                    Lproduct[i1+i2] = L.C->F->sum(Lproduct[i1+i2], L.C->F->product(ML(i1,i), MB(i2,j)));
            for (int i3=0; i3<MBL.n; i3++)
                M(i, j, i3) =  Lproduct[MBL.col_basis[i3]];
        };

    return M;
};


int Koszul(int p, int q, LineBundle L, LineBundle B)
{
    LineBundle B0, B1;

    B1 = B;

    int qq = q;
    while (qq>0){
        B1 = LBmult(B1, L);
        qq--;
    };
    while (qq<0){
        B1 = LBmult(B1, LBinverse(L));
        qq++;
    };
    B0 = LBmult(B1, LBinverse(L));

    Matrix21 M2 = multTable(L,B1);
    int r12 = run(p, M2);
    Matrix21 M1 = multTable(L,B0);
//    int r01 = run(p+1, M1);
    int r01 = choose(M1.a, p+1);

    int m1 = choose(M1.a, p) * M1.c;
    int m2 = choose(M2.a, p-1) * M2.c;

    cout << "Morphism = ";
    cout.width(5);
    cout << m1 << " x ";
    cout.width(5);
    cout << m2 << ":  ";
    cout << "K_(" << p << "," << q << ") = ";
    cout.width(5);
    cout << m1-r12-r01 << "\n";
//    cout << r01 << "," << r12 << "\n";

    return m1-r12-r01;
}


int Koszul(int p, int q, LineBundle L)
{
    return Koszul(p, q, L, LBmult(L,LBinverse(L)));
}


BettiTable::BettiTable(Curve* C0)
{
    C = C0;
    int g = C->genus;
    LineBundle L = C->canonical();

    Matrix21 M = multTable(L,L);
    int n = M.a;

    betti.resize(4 * (n-1));
    dim.resize(g * n);
    chi.resize(g);
    coimage.resize(4*(n-1));

    int sign = 1;
    for (int i=0; i<g; i++){
        int N = (i == 1) ? g : (2*i-1)*(g-1);
        if (i == 0) N = 1;
        for (int j=0; j<n; j++){
            dim[i*g + j] = N * choose(g, j);  // only i=0..3 relevant
            if (i+j < g) chi[i+j] += sign * N * choose(g, j);
        };
        sign = -sign;
    };

    betti[0] = betti[4*g-5] = 1;
    for (int p=1; 2*p<g; p++){
        int r01 = choose(M.a, p+1);
        int m1 = choose(M.a, p) * M.a;
        int r = m1 - r01 - run(p, M);
        betti[(g-1)+p] = betti[(g-1)*2+(g-1-p-1)] = r;
        betti[(g-1)*2+(p-1)] = betti[(g-1)*2-p] = r + chi[p+1];

        int m2 = choose(M.a, p-1) * M.c;
        cout << "Morphism = ";
        cout.width(5);
        cout << m1 << " x ";
        cout.width(5);
        cout << m2 << ":  ";
        cout << "K_(" << p << "," << "1" << ") = ";
        cout.width(5);
        cout << r << "\n";
    };
};


void BettiTable::print()
{
    int g = C->genus;

    std::cout << "\n     ";
    for (int p=0; p<g-1; p++){
        std::cout.width(4);
        std::cout << p << " ";
    };
    std::cout << "\n";
    std::cout << "     ";
    for (int p=0; p<g-1; p++){
        std::cout << "-----";
    };
    std::cout << "\n";

    for (int q=0; q<4; q++){
        std::cout.width(3);
        std::cout << q << " |";
        for (int p=0; p<g-1; p++){
            std::cout.width(4);
            std::cout << betti[q*(g-1)+p] << " ";
        };
        std::cout << "\n";
    };
};


int run(int t, Matrix21& M)
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
                        MM(srow, k2, trow, k3) = M.F->sum(MM(srow, k2, trow, k3), M(c[k], k2, k3));
                    else
                        MM(srow, k2, trow, k3) = M.F->sum(MM(srow, k2, trow, k3), M.F->neg(M(c[k], k2, k3)));
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

//    MM.gauss();
    MM.gauss2();

    return MM.rk;
}
