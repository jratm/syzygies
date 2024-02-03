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
        FCurve(Field* F0, int g) : genus(g) {
                F = F0; nodes.resize(g);
                if (F->new_way == false)
                    for (int i=0; i<g; i++) {  nodes[i].p = 2*i+1; nodes[i].q = 2*i+2; }
                else
                    for (int i=0; i<g; i++) {  nodes[i].p = F->rebase[2*i+1]; nodes[i].q = F->rebase[2*i+2];
                };
        }
        FMatrix morphism(FLineBundle&);
        FLineBundle canonical();
        FLineBundle point(int);
        int syzygy(int,int,FLineBundle,FLineBundle);
        int syzygy(int,int,FLineBundle);
        int run_K(FLineBundle&,FLineBundle&,FLineBundle&,int);
        int genus;
        Field* F;
    private:
        std::vector<node> nodes;
};

#endif // CURVE_H
