#ifndef CURVE_H
#define CURVE_H

#include "xmatrix.h"
#include <vector>
#include <random>
#include <functional>
#include <chrono>


int run(int, Matrix21&);

struct node
{
    int p, q;
};

class LineBundle;
LineBundle LBinverse(LineBundle);
LineBundle LBmult(LineBundle,LineBundle);

class Curve
{
    public:
        Curve(Field* F0, int g) : genus(g), F(F0) {
            nodes.resize(g);
            std::vector<node> Cnodes = sample_nodes(g, F0->q);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            };
        };
        Curve(Field* F0, int g, vector<node> Cnodes) : genus(g), F(F0) {
            nodes.resize(g);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            }
        };
        FMatrix sections(LineBundle&);
        LineBundle canonical();
        LineBundle trivial();
        LineBundle pt();
        LineBundle point(int);
        const int genus;
        Field* F;
        void print();
    private:
        std::vector<node> sample_nodes(int, int);
        std::vector<node> nodes;
};


class LineBundle
{
    public:
        LineBundle(Curve* C0, int d) : C(C0), degree(d) { ratios.resize(C0->genus); };
//        LineBundle() { sections = new FMatrix();};
        Curve* C;
        int degree;
        int dim;
        std::vector<int> ratios;
//        FMatrix* sections;
};


#endif // CURVE_H
