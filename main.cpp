#include <iostream>
#include <random>
#include <functional>
#include <chrono>
#include "curve.h"
#include "Number.h"


void test1();
void test2();



int main()
{

    test1();
//    test2();
    return 0;
}


void test1()
{
    Field F(3,5);
    F.print();

    int g = 14;
    Curve C(&F, g);
    C.print();

    BettiTable K(&C);
    K.print();
    return;
}


void test2()
{
    Field F(5,4);
    F.print();

    int g = 13;

    for (int i=0; i<20; i++){
        Curve C(&F, g);
//        C.print();
    BettiTable K(&C);
    K.print();
//        LineBundle L = C.canonical();
//        Koszul(3,1,L);
    };
//    std::vector<node> nodes(g);
//    for (int i=0; i<g; i++) {
//        nodes[i].p = 2*i+1;
//        nodes[i].q = 2*i+2;
//    };

//    for (int i=2*g; i<100; i++){
//        nodes[g-1].q = i;
//        Curve C(&F, g, nodes);
//        C.print();
//        std::cout << i << ":  ";
//        LineBundle L = C.canonical();
//        Koszul(5,1,L);
//    }

    return;
}
