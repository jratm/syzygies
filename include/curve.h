#ifndef CURVE_H
#define CURVE_H

#include "xmatrix.h"
#include <vector>


struct node
{
    int p, q;
};

class FCurve;

struct FLineBundle
{
    int degree;
    int dim;
    std::vector<int> ratios;
    FCurve* C;
};

class FCurve
{
    public:
        FCurve(Field* F0, int g, vector<node> Cnodes) : genus(g), F(F0), nodes(Cnodes) {};
        FMatrix sections(FLineBundle&);
        FLineBundle canonical();
        FLineBundle point(int);
        int syzygy(int,int,FLineBundle,FLineBundle);
        int syzygy(int,int,FLineBundle);
        int run_K(FLineBundle&,FLineBundle&,FLineBundle&,int);
        int genus;
        Field* F;
        void print();
    private:
        std::vector<node> nodes;
};

class BettiTable
{
    public:
        BettiTable(FCurve*);
        FCurve* C;
        void print();
    private:
        std::vector<int> betti;
        std::vector<int> dim;
        std::vector<int> chi;
        std::vector<int> coimage;
        int run(int,FMatrix21&);
        int Koszul(int,int,FLineBundle);
        int Koszul(int,int,FLineBundle,FLineBundle);
};

#endif // CURVE_H
