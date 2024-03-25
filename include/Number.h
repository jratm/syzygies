#pragma once

#include <vector>
#include <cstdint>

using SHORT = int16_t;

class Field
{
    public:
        const int p, f, q;
        SHORT inverse(SHORT) const;
        SHORT neg(SHORT) const;
        SHORT product(SHORT,SHORT) const;
        SHORT sum(SHORT,SHORT) const;
        std::vector<SHORT> encode;   //
        std::vector<int> decode;   //
        Field(int,int);
        void print() const;
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
        std::vector<int> mult(const std::vector<int>&, const std::vector<int>& ) const;
        int p_inverse(int) const;
        int gcd(const std::vector<int>&) const;
        std::vector<int> power(const std::vector<int>&, int) const;
    friend class Matrix22;
};
