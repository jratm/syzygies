#include "Number.h"
#include <iostream>

static inline int multiplied( const int p, const int f )
{
    int q = 1;
    for (int i=0; i<f; i++) q *= p;

    return q;
}

Field::Field(const int p0, const int f0)
    : m_p( p0 )
    , m_f( f0 )
    , m_q( multiplied( p0, f0 ) )
{
    m_poly.resize(m_f+1);
    create();
    generator();
    tables();
};


int Field::gcd(const std::vector <int>& h) const // return gcd of poly and h
{
    std::vector <int> a = m_poly;
    std::vector <int> b = h;

    int lead_a = m_f;
    int lead_b = m_f-1;

    while (true) {

    // determine leading term of b
        for (; lead_b>=0; lead_b--)
            if (b[lead_b] != 0) break;
        if (lead_b == -1) { return 0; }; // failure: b is a common factor of poly and h
        if (lead_b == 0)  { return lead_a;  }; // success: poly and h have no common factor

    //  make lead coefficient 1
        int mm = p_inverse(b[lead_b]);
        for (int i=0; i<=lead_b; i++) b[i] = (b[i] *mm) % m_p;

    //  get remainder of division
        for (int i=lead_a; i>=lead_b; i--){
            int factor = a[i];
            if (factor != 0){
                for (int j=0; j<=lead_b; j++){
                    int term = (a[i+j-lead_b] - factor * b[j]) % m_p;
                    if (term < 0) term += m_p;
                    a[i+j-lead_b] = term;
                };
            };
        };
    // exchange a and b
        for (int i=0; i<=lead_b; i++)
            std::swap( a[i], b[i] );
        lead_a = lead_b;
        lead_b --;
    };
    return -1;
};


std::vector<int> Field::power(const std::vector <int>& h, const int n) const // return h^n mod poly
{
    std::vector<int> result(2*m_f+1); result[0]=1;
    std::vector<int> sq(2*m_f+1);

    for (int i=0; i<=m_f; i++) sq[i] = h[i];

    if (n % 2 != 0) result = mult(result,sq);

// calculate x^p mod poly by repeated squaring
    for (int f0=n/2; f0>0; f0 /= 2){
        sq = mult(sq,sq);
        if (f0 % 2 != 0) result = mult(result,sq);
    };

    return result;
};


int Field::p_inverse(int x) const
{
    int a=1, b=0, g=m_p, u=0, v=1, w=x;
    int q, t;

    if (w<0) w += m_p;
    while (w>0)
    {
        q=g/w;
        t=u; u=a-q*u; a=t;
        t=v; v=b-q*v; b=t;
        t=w; w=g-q*w; g=t;
    }
    if (b<0) b += m_p;
    return b % m_p;
}


std::vector<int> Field::mult(const std::vector<int>& a, const std::vector<int>& b) const
{
    std::vector<int> result(2*m_f+1);

    for (int i=0; i<m_f; i++)
        for (int j=0; j<m_f; j++)
            result[i+j] = (result[i+j] + a[i]*b[j]) % m_p;
// reduce modulo poly
    for (int i=2*m_f; i>=m_f; i--){
        int factor = result[i];
        if (factor != 0){
            for (int j=0; j<=m_f; j++){
                int term = (result[i+j-m_f] - factor * m_poly[j]) % m_p;
                if (term < 0) term += m_p;
                result[i+j-m_f] = term;
            };
        };
    };

    return result;
}


void Field::create()
{
// calculate with polynomials of degree up to f
    m_poly[m_f] = 1;        // set leading term
    bool reducible;

    for (int n=0; n<m_q; n++)     // loop over all possible lower terms
    {
        int n0=n;
        for (int i=0; i<m_f; i++){
            m_poly[i] = n0 % m_p;
            n0 /= m_p;
        };

        std::vector<int> res(2*m_f+1); res[1]=1;  // res = X
        std::vector<int> test(2*m_f+1);

        reducible = false;

        for (int j=1; 2*j<=m_f; j++)  //  test for factor of degree j
        {
            res = power(res,m_p);
            test = res; test[1]--; if (test[1] < 0) test[1] = m_p-1;

            if ( gcd(test) == 0){
                reducible = true;
                break;
            };
        };

        if (reducible == false) break;
    };
};


void Field::generator()
{
    std::vector<int> primes;
    m_gen.resize(m_f);

// determine prime factors of q-1
    int remainder = m_q-1;

    for (int i=2; i*i<=remainder; i++)
        if (remainder % i == 0){
            primes.push_back(i);
            while (remainder % i == 0) remainder /= i;
        };
    if (remainder > 1) primes.push_back(remainder);

//    std::cout << "factors of " << q-1 << " are: ";
//    for (auto x : primes) std::cout << x << ", ";
//    std::cout << "\n";

    for (int n=1; n<m_q; n++)     // loop over all field elements
    {
        int n0=n;
        for (int i=0; i<m_f; i++){
            m_gen[i] = n0 % m_p;
            n0 /= m_p;
        };

        bool potential_generator = true;
    // calculate gen^((q-1)/primes[j]) mod poly
        for (auto x : primes){
            std::vector<int> res(2*m_f+1);
            res = power(m_gen, (m_q-1)/x);
    // if res == 1 report failure
            res[0]--; if (res[0]<0) res[0] = m_p-1;
            bool power_is_one = true;
            for (int k=0; k<m_f; k++) if ( res[k] != 0 )  { power_is_one = false; break; };
            if (power_is_one == false) continue;  // (q-1)/p_j -th power not = 1, proceed to next factor
    //  now (q-1)/p_j -th power is = 1, so n does not generate the multiplicative group of F
            potential_generator = false;
            break;      // no need to check other primes; continue with next higher n
        };
        if (potential_generator == true) return;
    };
    std::cout << "failed to find multiplicative generator\n";
};


void Field::tables()
{
    const int p1 = 2*m_p-1; int q1=1;
    for (int i=0; i<m_f; i++) q1 *= p1;
    m_encode.resize(m_q);  // p-respresentation --> p1-representation
    m_decode.resize(q1);  // dirty p1-representation --> (clean) p-representation
    m_normalize.resize(q1);  // dirty p1-representation --> clean p1-representation


    std::vector<int> y(m_f);
    for (int index=0; index<m_q; index++) {
        int z = index;
        for (int i=0; i<m_f; i++){
            y[i] = z % m_p;
            z /= m_p;
        };
        int n = 0; for (int i=m_f-1; i>=0; i--) n = n * p1 + y[i];
        m_encode[index] = n;
    };
    for (int index=0; index <q1; index++){
        int z = index;
        for (int i=0; i<m_f; i++){
            y[i] = z % p1;
            z /= p1;
        };
        int n = 0; for (int i=m_f-1; i>=0; i--) n = n * m_p + y[i] % m_p;
        m_decode[index] = n;
//        int n = 0; for (int i=f-1; i>=0; i--) n = n * p1 + y[i] % p1;
//        unbase0[index] = n;
    };
    for (int index=0; index<q1; index++){
        m_normalize[index] = m_encode[m_decode[index]];
    };

    m_exp.resize(3*m_q);
    m_log.resize(q1);

    std::vector<int> x = m_gen;
    m_negative.resize(q1);

    m_exp[0] = 1;
    for (int index=1; index<m_q; index++) {
        // find corresponding integer
        int n = 0; for (int i=m_f-1; i>=0; i--) n = n*m_p + x[i];  //  ???
        int m = 0; for (int i=m_f-1; i>=0; i--) m = m*m_p + ( x[i]==0 ? 0 : m_p-x[i]); // ??????
        m_exp[index] = m_exp[index+m_q-1] = m_exp[index + 2*m_q - 2] = m_encode[n];
        m_log[m_encode[n]] = index;
        x = mult(x,m_gen);
        m_negative[m_encode[n]] = m_encode[m];
    };

    m_log[1]=0;

    return;
};

void Field::print() const
{
    std::cout << "Base field properties:\n";
    std::cout << "     Characteristc        p = " << m_p << "\n";
    std::cout << "     Exponent             f = " << m_f << "\n";
    std::cout << "     Number of elements   q = " << m_q << "\n";
    std::cout << "     Splitting field of the polynomial " ;
    std::cout << "f =  X^" << m_f;
    for (int i=m_f-1; i>=0; i-- ) if (m_poly[i] != 0) {
        std::cout << " + ";
        if (m_poly[i]!=1 || i==0) std::cout << m_poly[i];
        if (i>0){
            std::cout << " X";
            if (i>1){
                std::cout << "^" << i;
            };
        };
    };
    std::cout << "\n";

    bool first_term = true;
    std::cout << "     Generator of the multiplicative group:  g = ";
    for (int i=m_f-1; i>=0; i-- ) if (m_gen[i] != 0) {
        if (first_term == false) std::cout << " + ";
        first_term = false;
        if (m_gen[i]!=1 || i==0) std::cout << m_gen[i];
        if (i>0){
            std::cout << " X";
            if (i>1){
                std::cout << "^" << i;
            };
        };
    };
    std::cout << "\n\nNOTE: Field elements correspond to degree " << m_f-1 << " polynomials over the prime\n";
    std::cout << "      field Z/" << m_p <<". Writing these coefficients in sequence, the resulting word\n";
    std::cout << "      represents a unique basis " << m_p << " number in the range between 0 and " << m_q-1 << ".\n\n";

    return;
};
