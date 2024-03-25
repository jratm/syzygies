#pragma once

#include "curve.h"

class LineBundle;


int Koszul(int,int,LineBundle&);
int Koszul(int,int,LineBundle&,LineBundle);


class BettiTable
{
    public:
        BettiTable(Curve*);
        BettiTable(Curve*,LineBundle);
        Curve* C;
        void print();
    private:
        void initialize(LineBundle);
        int h0;  // h0 of line bundle L
        std::vector<int> betti;
        std::vector<int> dim;
        std::vector<int> chi;
        std::vector<int> coimage;
};
