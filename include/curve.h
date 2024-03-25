#pragma once

#include "xmatrix.h"
#include <vector>
#include <random>
#include <functional>
#include <chrono>

int run(int, const Matrix21&);

struct node
{
    int p, q;
};

class LineBundle;

LineBundle LBinverse(const LineBundle&);
LineBundle LBmult(const LineBundle&, const LineBundle&);

class Curve
{
    public:
        Curve(const Field* F0, int g) : genus(g), F(F0) {
            nodes.resize(g);
            const std::vector<node> Cnodes = sample_nodes(g, F0->q);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            };
        };
        Curve(const Field* F0, int g, const vector<node>& Cnodes) : genus(g), F(F0) {
            nodes.resize(g);
            for (int i=0; i<g; i++) {
                nodes[i].p = F->encode[Cnodes[i].p];
                nodes[i].q = F->encode[Cnodes[i].q];
            }
        };
        FMatrix sections(const LineBundle&) const;
        LineBundle canonical() const ;
        LineBundle trivial() const;
        LineBundle pt() const;
        LineBundle point(int) const;
        LineBundle modify(int,int) const;
        const int genus;
        const Field* const F;
        void print() const;
    private:
        std::vector<node> sample_nodes(int, int);
        std::vector<node> nodes;
};


class LineBundle
{
    public:
        LineBundle(const Curve* C0, int d) : C(C0), degree(d) { ratios.resize(C0->genus); };
//        LineBundle() { sections = new FMatrix();};
        const Curve* C;
        int degree;
        int dim;
        std::vector<int> ratios;
//        FMatrix* sections;
};
