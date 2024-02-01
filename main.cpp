#include <iostream>
#include "curve.h"
#include "Number.h"


//using namespace std;

long choose(int a, int b);


int main()
{
    int g = 13;

    Field F(3,5);
    FCurve C(&F, g);
    FLineBundle L = C.canonical();
    std::cout << "\ngenus = " << g << "\n\n";
//    for (int p=1; p<=g-2; p++) for (int q=1; q<=3; q++)
//        C.syzygy(p,q,L);

    std::vector<int> chi(g);
    int sign = 1;
    for (int i=0; i<g; i++){
        int N = (i == 1) ? g : (2*i-1)*(g-1);
        if (i == 0) N = 1;
        for (int j=0; j<g; j++){
            if (i+j < g) chi[i+j] += sign * N * choose(g, j);
        };
        sign = -sign;
    };
    for (int i=0; i<g; i++ ) std::cout << chi[i] << " ";
    std::cout << "\n\n";

    std::vector<int> syz(4*(g-1));
    syz[0] = syz[4*g-5] = 1;
    for (int p=1; 2*p<g; p++){
        int r = C.syzygy(p,1,L);
        syz[(g-1)+p] = syz[(g-1)*2+(g-1-p-1)] = r;
        syz[(g-1)*2+(p-1)] = syz[(g-1)*2-p] = r + chi[p+1];
    };

    std::cout << "\n     ";
    for (int p=0; p<g-1; p++){
        std::cout.width(4);
        std::cout << p << " ";
    };
    std::cout << "\n";
    std::cout << "     ";
    for (int p=0; p<g-1; p++){
        std::cout << "-----";
    };
    std::cout << "\n";

    for (int q=0; q<4; q++){
        std::cout.width(3);
        std::cout << q << " |";
        for (int p=0; p<g-1; p++){
            std::cout.width(4);
            std::cout << syz[q*(g-1)+p] << " ";
        };
        std::cout << "\n";
    };

    std::cout << "\ncalculation complete\n";
    return 0;
}

