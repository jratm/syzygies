#pragma once

#include <vector>
#include <cstdint>

using SHORT = int16_t;

class Field
{
    public:
        int p, f, q;
        SHORT inverse(SHORT);
        SHORT neg(SHORT);
        SHORT product(SHORT,SHORT);
        SHORT sum(SHORT,SHORT);
        std::vector<SHORT> encode;   //
        std::vector<int> decode;   //
        Field(int,int);
        void print();
    private:
        std::vector<int> poly;  //  F is the splitting field of poly
        std::vector<int> gen;   //  gen is a generator of the multiplicative group of F
        void create();
        void generator();
        void tables();
        std::vector<SHORT> exp;   //
        std::vector<SHORT> log;   //
        std::vector<SHORT> negative;
        std::vector<SHORT> normalize;   //
        std::vector<int> mult(std::vector<int>, std::vector<int>);
        int p_inverse(int);
        int gcd(std::vector<int>);
        std::vector<int> power(std::vector<int>, int);
    friend class Matrix22;
};
