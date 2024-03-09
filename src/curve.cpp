#include "curve.h"
#include "xmatrix.h"
#include <iostream>


using namespace std;


LineBundle LBmult(LineBundle L0, LineBundle L1)
{
    LineBundle res(L0.C, L0.degree + L1.degree);
    for (int i=0; i<L0.C->genus; i++) res.ratios[i] = L0.C->F->product(L0.ratios[i], L1.ratios[i]);
    return res;
};


LineBundle LBinverse(LineBundle L0)
{
    LineBundle res(L0.C, -L0.degree);
    for (int i=0; i<L0.C->genus; i++) res.ratios[i] = L0.C->F->inverse(L0.ratios[i]);
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

    LineBundle kan(this, 2 * genus - 2);

    for (i=0; i<genus; i++)
    {
        int res = F->product(F->neg(C(2*i+1,0)), F->inverse(C(2*i,0)));
        kan.ratios[i] = res;
    }

    return kan;
}


LineBundle Curve::trivial()
{
    LineBundle Lp(this, 0);

    for (int i=0; i<genus; i++) Lp.ratios[i] = 1;

    return Lp;
}

LineBundle Curve::pt()
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto rand = bind(uniform_real_distribution<double>{0,1},default_random_engine(seed));

    int x;
    while (true)
    {
        double u = rand();
        x = floor(F->q*u);
        bool success = true;
        for (int i=0; i<genus; i++)
        {
            if (x == F->decode[nodes[i].p]) success = false;
            if (x == F->decode[nodes[i].q]) success = false;
        };
        if (  success == true) break;
    };

    return point(x);
}


LineBundle Curve::point(int p0)
{
    LineBundle Lp(this, 1);

    int p1 = F->encode[p0];

    for (int i=0; i<genus; i++)
    {
        int res = F->product(F->sum(nodes[i].p, F->neg(p1)), F->inverse(F->sum(nodes[i].q, F->neg(p1))));
        Lp.ratios[i] = res;
    }

    return Lp;
}


LineBundle Curve::modify(int a, int b)
{
    LineBundle L = trivial();

    for (int i=0; i<a; i++) L = LBmult(L, pt());
    for (int i=0; i<b; i++) L = LBmult(L, LBinverse(pt()));

    return L;
}



FMatrix Curve::sections(LineBundle& L)
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
    std::cout << "     Nodes in the following points: \n";
    for (int i=0; i<genus; i++){
        if (i % 6 == 0) std::cout << "          ";
        std::cout << "(" << F->decode[nodes[i].p] << "," << F->decode[nodes[i].q] << ")";
        if (i != genus-1) std::cout << ", ";
        if (i % 6 == 5) std::cout << "\n";
    };
    std::cout << "\n\n";

    return;
};


