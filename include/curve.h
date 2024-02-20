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

class Curve;

struct LineBundle
{
    int degree;
    int dim;
    std::vector<int> ratios;
    Curve* C;
};


int Koszul(int,int,LineBundle);
int Koszul(int,int,LineBundle,LineBundle);


class Curve
{
    public:
        Curve(Field*,int);
        Curve(Field* F0, int g, vector<node> Cnodes) : genus(g), F(F0) {
            nodes.resize(g);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            }
        };
        FMatrix sections(LineBundle);
        LineBundle canonical();
        LineBundle point(int);
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
        BettiTable(Curve*);
        Curve* C;
        void print();
    private:
        std::vector<int> betti;
        std::vector<int> dim;
        std::vector<int> chi;
        std::vector<int> coimage;
};

#endif // CURVE_H
