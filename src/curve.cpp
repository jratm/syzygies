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


FLineBundle LBmult(FLineBundle L_0, FLineBundle L_1)
{
    FLineBundle res;
    res.degree = L_0.degree + L_1.degree;
    res.ratios.resize(L_0.C->genus);
    res.C = L_0.C;
    for (int i=0; i<L_0.C->genus; i++) res.ratios[i] = L_0.C->F->product(L_0.ratios[i], L_1.ratios[i]);
    return res;
};


FLineBundle LBinverse(FLineBundle L_0)
{
    FLineBundle res;
    res.degree = -L_0.degree;
    res.ratios.resize(L_0.C->genus);
    res.C = L_0.C;
    for (int i=0; i<L_0.C->genus; i++) res.ratios[i] = L_0.C->F->inverse(L_0.ratios[i]);
    return res;
}


FLineBundle FCurve::canonical()
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

    FMatrix C = A.GaussJordan().Nullspace();

    FLineBundle kan;
    kan.degree = deg;
    kan.ratios.resize(genus);
    kan.C = this;

    for (i=0; i<genus; i++)
    {
        int res = F->product(F->neg(C(2*i+1,0)), F->inverse(C(2*i,0)));
        kan.ratios[i] = res;
    }
//    morphism(kan);

    return kan;
}


FLineBundle FCurve::point(int p0)
{
    FLineBundle Lp;
    int p1 = F->encode[p0];

    Lp.degree = 1;
    Lp.ratios.resize(genus);
    Lp.C = this;

    for (int i=0; i<genus; i++)
    {
        int res = F->product(F->sum(nodes[i].p, F->neg(p1)), F->inverse(F->sum(nodes[i].q, F->neg(p1))));
        Lp.ratios[i] = res;
    }

    morphism(Lp);

    return Lp;
}


FMatrix FCurve::morphism(FLineBundle& L)
/***
Each of the g nodes imposes one condition on the degree d=deg functions
on P^1. Each row of A contains the conditions on the coefficients of
such a function, using the standard basis for H^0(O(d)).

return value C: The columns of C correspond to degree d functions on P^1
respecting the prescribed nodes.  Assuming the nodes impose
independent conditions, there are r+1 of them, providing a map to P^r.
C has size ( d+1, r+1 ), and the image of the point x (on P^1) is C.x;
a column (a_0 a_1 ...)^T of C encodes the rational
function a_0 + a_1 t + a_2 t^2 + ... on P^1
***/
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

    FMatrix C = A.GaussJordan().Nullspace();
    L.dim = C.n;

    return C;
};


int FCurve::syzygy(int p, int q, FLineBundle L, FLineBundle B)
{
    FLineBundle B1, B2;

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
    B2 = LBmult(B1, L);

    int r12 = run_K(L, B1, B2, p);
    int r01 = choose(L.dim, p+1);

    int m1 = choose(L.dim,p)*B1.dim;
    int m2 = choose(L.dim,p-1)*B2.dim;

    cout << "Morphism = ";
    cout.width(5);
    cout << m1 << " x ";
    cout.width(5);
    cout << m2 << ":  ";
    cout << "K_(" << p << "," << q << ") = " << m1-r12-r01 << "\n";

    return m1-r12-r01;
}


int FCurve::syzygy(int p, int q, FLineBundle L)
{
    return syzygy(p, q, L, LBmult(L,LBinverse(L)));
}


int FCurve::run_K(FLineBundle& L, FLineBundle& B, FLineBundle& LB, int t)
{

    FMatrix ML = morphism(L);
    FMatrix MB = morphism(B);
    FMatrix MBL = morphism(LB);

// precalculate the multiplication table  H^0(L) x H^0(B) --> H^0(L tensor B)
    FMatrix M(F, L.dim, B.dim, LB.dim);    //      M has size  ( g x (3g-3) x (5g-5) )

    for (int i=0; i<L.dim; i++)
        for (int j=0; j<B.dim; j++)
        {
            std::vector<int> Lproduct(LB.degree+1);
            for (int i1=0; i1<=L.degree; i1++)
                for (int i2=0; i2<=B.degree; i2++)
                    Lproduct[i1+i2] = F->sum(Lproduct[i1+i2], F->product(ML(i1,i), MB(i2,j)));
            for (int i3=0; i3<LB.dim; i3++)
                M(i, j, i3) =  Lproduct[MBL.col_basis[i3]];
        };

/*
generate all combinations, using Algorithm T from Knuth, TAOCP Vol. 4B, p. 359
*/
// n = L.dim, t as function parameter
    std::vector<int> c(t+2);
// step_T1:
    for (int l=0; l<t; l++) c[l] = l;
    c[t] = L.dim;
    c[t+1] = 0;

    int l = t-1;
    int x;

    FMatrix MM(F, choose(L.dim,t), B.dim, choose(L.dim,t-1), LB.dim);
    int srow = 0;

    while (true)
    {
// run work
//        int srow = 0;
//        for (int k1=0; k1<t; k1++)
//            srow += choose(c[k1],k1+1);

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
            for (int k2=0; k2<B.dim; k2++)
                for (int k3=0; k3<LB.dim; k3++)
                {
                    if (sign)
                        MM(srow, k2, trow, k3) = F->sum(MM(srow, k2, trow, k3), M(c[k], k2, k3));
                    else
                        MM(srow, k2, trow, k3) = F->sum(MM(srow, k2, trow, k3), F->neg(M(c[k], k2, k3)));
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

    MM.Gauss();
//    MM.GaussJordan();
    return MM.rk;
}
