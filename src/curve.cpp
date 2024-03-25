#include "curve.h"
#include "xmatrix.h"
#include <iostream>


using namespace std;


LineBundle LBmult( const LineBundle& L0, const LineBundle& L1)
{
    LineBundle res(L0.C, L0.degree + L1.degree);

    for (int i=0; i<L0.C->genus(); i++)
        res.ratios[i] = L0.C->F->product(L0.ratios[i], L1.ratios[i]);

    return res;
};


LineBundle LBinverse(const LineBundle& L0)
{
    LineBundle res(L0.C, -L0.degree);
    for (int i=0; i<L0.C->genus(); i++)
        res.ratios[i] = L0.C->F->inverse(L0.ratios[i]);

    return res;
}

Curve::Curve(const Field* F0, int g)
    : Curve( F0, g, sample_nodes(g, F0->elementCount() ) )
{
}

Curve::Curve(const Field* F0, int g, const vector<Node>& Cnodes)
    : F(F0)
    , m_genus(g)
{ 
    nodes.resize(g); 

    for (int i=0; i<g; i++) {
        nodes[i].p = F->encoded( Cnodes[i].p );
        nodes[i].q = F->encoded( Cnodes[i].q );
    } 
};

LineBundle Curve::canonical() const
{
    const int deg = 2 * m_genus - 2;
    FMatrix A(F, 2 * m_genus - 1, 2 * m_genus);

    for ( int i=0; i<m_genus; i++)
    {
        int a = 1;
        int b = 1;

        for ( int j=0; j<deg+1; j++)
        {
            A(j,2*i) = a;
            A(j,2*i+1) = b;
            a = F->product(a, nodes[i].p);
            b = F->product(b, nodes[i].q);
        }
    };

    const FMatrix C = A.gauss_jordan().nullspace();

    LineBundle kan(this, 2 * m_genus - 2);

    for ( int i=0; i<m_genus; i++)
    {
        const int res = F->product(F->neg(C(2*i+1,0)), F->inverse(C(2*i,0)));
        kan.ratios[i] = res;
    }

    return kan;
}


LineBundle Curve::trivial() const
{
    LineBundle Lp(this, 0);

    for (int i=0; i<m_genus; i++) Lp.ratios[i] = 1;

    return Lp;
}

LineBundle Curve::pt() const
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto rand = bind(uniform_real_distribution<double>{0,1},default_random_engine(seed));

    int x;
    while (true)
    {
        const double u = rand();

        x = floor(F->elementCount() * u);
        bool success = true;
        for (int i=0; i<m_genus; i++)
        {
            const auto& node = nodes[i];

            if ( ( x == F->decoded( node.p ) ) || ( x == F->decoded( node.q ) ) )
                success = false;
        };
        if ( success == true) break;
    };

    return point(x);
}


LineBundle Curve::point(int p0) const
{
    LineBundle Lp(this, 1);

    const auto p1 = F->encoded(p0);

    for (int i=0; i<m_genus; i++)
    {
        const int res = F->product(F->sum(nodes[i].p, F->neg(p1)), F->inverse(F->sum(nodes[i].q, F->neg(p1))));
        Lp.ratios[i] = res;
    }

    return Lp;
}


LineBundle Curve::modify(const int a, const int b) const
{
    LineBundle L = trivial();

    for (int i=0; i<a; i++) L = LBmult(L, pt());
    for (int i=0; i<b; i++) L = LBmult(L, LBinverse(pt()));

    return L;
}



FMatrix Curve::sections(const LineBundle& L) const
{
    const int deg = L.degree;
    FMatrix A(F, m_genus, deg+1);

    for ( int i=0;i<m_genus;i++)
    {
        const int k = nodes[i].p;
        const int l = nodes[i].q;
        const int c = L.ratios[i];

        int a = 1;
        int b = 1;

        for ( int j=0;j<deg+1;j++)
        {
            A(i,j) = F->sum(a, F->neg(F->product(b, c)));
            a = F->product(a, k);
            b = F->product(b, l);
        }
    }

    return A.gauss_jordan().nullspace();
};


std::vector<Curve::Node> Curve::sample_nodes(const int g, const int q)
{
    const int n = 2*g;
    const int N = q-2;

    std::vector<int> sample(n);
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
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

    std::vector<Node> nodes(g);
    for (int i=0; i<g; i++)
        nodes[i] = { sample[2*i], sample[2*i+1] };

    return nodes;
}


void Curve::print() const
{
    std::cout << "Curve properties:\n";
    std::cout << "     genus = " << m_genus << "\n";
    std::cout << "     Nodes in the following points: \n";
    for (int i=0; i<m_genus; i++){

        const auto& node = nodes[i];

        if (i % 6 == 0) std::cout << "          ";
        std::cout << "(" << F->decoded( node.p ) << "," << F->decoded( node.q ) << ")";
        if (i != m_genus-1) std::cout << ", ";
        if (i % 6 == 5) std::cout << "\n";
    };
    std::cout << "\n\n";

    return;
};


