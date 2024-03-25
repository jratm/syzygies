#pragma once

#include <vector>
#include <cstdint>

using SHORT = int16_t;

class Field
{
    public:
        Field(int,int);

        SHORT log( SHORT index ) const;
        SHORT exp( SHORT index ) const;
        SHORT normalized( SHORT index ) const;
        SHORT neg(SHORT) const;

        SHORT inverse(SHORT) const;
        SHORT product(SHORT,SHORT) const;
        SHORT sum(SHORT,SHORT) const;

        void print() const;

        inline int elementCount() const { return m_q; }

        inline int decoded( int index ) const { return m_decode[index]; }
        inline SHORT encoded( int index ) const { return m_encode[index]; }

    private:
        void create();
        void generator();
        void tables();

        std::vector<int> mult(const std::vector<int>&, const std::vector<int>& ) const;

        int p_inverse(int) const;
        int gcd(const std::vector<int>&) const;

        std::vector<int> power(const std::vector<int>&, int) const;

        std::vector<SHORT> m_encode;
        std::vector<int> m_decode;

        std::vector<SHORT> m_exp;   //
        std::vector<SHORT> m_log;   //
        std::vector<SHORT> m_negative;
        std::vector<SHORT> m_normalize;   //

        std::vector<int> m_poly;  //  F is the splitting field of poly
        std::vector<int> m_gen;   //  gen is a generator of the multiplicative group of F

        const int m_p, m_f, m_q;
};

inline SHORT Field::log( const SHORT index ) const
{
    return m_log[index];
}

inline SHORT Field::normalized( const SHORT index ) const
{
    return m_normalize[ index ];
}

inline SHORT Field::exp( const SHORT index ) const
{
    return m_exp[ index ];
}

inline SHORT Field::sum(const SHORT x, const SHORT y) const
{
    return m_normalize[x+y];
}

inline SHORT Field::product(const SHORT x, const SHORT y) const
{
    if (x == 0 || y == 0) return 0;
    return m_exp[m_log[x]+m_log[y]];
}

inline SHORT Field::inverse(const SHORT x) const
{
    return m_exp[ m_q - 1 - m_log[x] ];
}

inline SHORT Field::neg(const SHORT x) const
{
    return m_negative[x];
}


