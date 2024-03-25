#include "curve.h"
#include "betti.h"
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
    const Field F(3,5);
    F.print();

    const Curve C(&F, 12);
    C.print();

//    BettiTable K(&C);
    const LineBundle L0 = C.canonical();
    const LineBundle L1 = C.modify(9, 6);
    const LineBundle L = LBmult(L0, L1);
    const BettiTable K(&C, L);
    K.print();

//    Field F1(5,3);
//    F1.print();
//
//    Curve C1(&F1, 17);
//    C1.print();
//
//    BettiTable K1(&C1);
//    K1.print();

    return;
}


void test2()
{
    const Field F(3,5);
    F.print();

    for (int i=0; i<2; i++){
        const Curve C(&F, 17);
        C.print();
        const BettiTable K(&C);
        K.print();
//        const LineBundle L = C.canonical();
//        Koszul(3,1,L);
    };
//    std::vector<node> nodes(g);
//    for (int i=0; i<g; i++) {
//        nodes[i].p = 2*i+1;
//        nodes[i].q = 2*i+2;
//    };

//    for (int i=2*g; i<100; i++){
//        nodes[g-1].q = i;
//        const Curve C(&F, g, nodes);
//        C.print();
//        std::cout << i << ":  ";
//        const LineBundle L = C.canonical();
//        Koszul(5,1,L);
//    }

    return;
}
