#pragma once

#include "curve.h"

class LineBundle;

int Koszul(int,int,const LineBundle&);
int Koszul(int,int,const LineBundle&,const LineBundle&);

class BettiTable
{
    public:
        BettiTable(const Curve*);
        BettiTable(const Curve*, const LineBundle&);
        const Curve* const C;
        void print() const;
    private:
        void initialize(const LineBundle&);
        int h0;  // h0 of line bundle L
        std::vector<int> betti;
        std::vector<int> dim;
        std::vector<int> chi;
        std::vector<int> coimage;
};
