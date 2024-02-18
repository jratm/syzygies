#ifndef CURVE_H
#define CURVE_H

#include "xmatrix.h"
#include <vector>
#include <random>
#include <functional>
#include <chrono>


int run(int, FMatrix21&);


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


int Koszul(int,int,FLineBundle);
int Koszul(int,int,FLineBundle,FLineBundle);


class FCurve
{
    public:
        FCurve(Field*,int);
        FCurve(Field* F0, int g, vector<node> Cnodes) : genus(g), F(F0) {
            nodes.resize(g);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            }
        };
        FMatrix sections(FLineBundle);
        FLineBundle canonical();
        FLineBundle point(int);
        int genus;
        Field* F;
        void print();
    private:
        std::vector<node> sample_nodes(int, int);
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
};

#endif // CURVE_H
