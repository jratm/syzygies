#pragma once

#include "xmatrix.h"
#include <vector>
#include <random>
#include <functional>
#include <chrono>

int run(int, const Matrix21&);

class LineBundle;

LineBundle LBinverse(const LineBundle&);
LineBundle LBmult(const LineBundle&, const LineBundle&);

class Curve
{
        struct Node { int p, q; };

    public:
        Curve(const Field* F0, int g);
        Curve(const Field* F0, int g, const std::vector<Node>& Cnodes);

        FMatrix sections(const LineBundle&) const;
        LineBundle canonical() const ;
        LineBundle trivial() const;
        LineBundle pt() const;
        LineBundle point(int) const;
        LineBundle modify(int,int) const;

        void print() const;

        const Field* const F;

        inline int genus() const { return m_genus; }

    private:
        const int m_genus;

        std::vector<Node> sample_nodes(int, int);
        std::vector<Node> nodes;
};


class LineBundle
{
    public:
        LineBundle(const Curve* C0, int d)
            : C(C0)
            , degree(d)
        {
            ratios.resize(C0->genus());
        }

        const Curve* C;
        int degree;
        std::vector<int> ratios;
};
