#include "Number.h"
#include <iostream>


Field::Field(int p0, int f0)
{
    p = p0; f = f0;
    q = 1;
    for (int i=0; i<f; i++) q *= p;
    poly.resize(f+1);
    create();
    generator();
    tables();
};


int Field::gcd(std::vector <int> h) // return gcd of poly and h
{
    std::vector<int> a(f+1);
    std::vector<int> b(f+1);
    std::vector<int> c(f+1);

    for (int i=0; i<=f; i++) a[i]=poly[i];
    int lead_a = f;

    for (int i=0; i<f; i++) b[i] = h[i];
    int lead_b = f-1;

    while (true) {

    // determine leading term of b
        for (; lead_b>=0; lead_b--)
            if (b[lead_b] != 0) break;
        if (lead_b == -1) { return 0; }; // failure: b is a common factor of poly and h
        if (lead_b == 0)  { return lead_a;  }; // success: poly and h have no common factor

    //  make lead coefficient 1
        int mm = p_inverse(b[lead_b]);
        for (int i=0; i<=lead_b; i++) b[i] = (b[i] *mm) % p;

    //  get remainder of division
        for (int i=lead_a; i>=lead_b; i--){
            int factor = a[i];
            if (factor != 0){
                for (int j=0; j<=lead_b; j++){
                    int term = (a[i+j-lead_b] - factor * b[j]) % p;
                    if (term < 0) term += p;
                    a[i+j-lead_b] = term;
                };
            };
        };
    // exchange a and b
        for (int i=0; i<=lead_b; i++) c[i]=a[i];
        for (int i=0; i<=lead_b; i++) a[i]=b[i];
        for (int i=0; i<=lead_b; i++) b[i]=c[i];
        lead_a = lead_b;
        lead_b --;
    };
    return -1;
};


std::vector<int> Field::power(std::vector <int> h, int n) // return h^n mod poly
{
    std::vector<int> result(2*f+1); result[0]=1;
    std::vector<int> sq(2*f+1);

    for (int i=0; i<=f; i++) sq[i] = h[i];

    if (n % 2 != 0) result = mult(result,sq);

// calculate x^p mod poly by repeated squaring
    for (int f0=n/2; f0>0; f0 /= 2){
        sq = mult(sq,sq);
        if (f0 % 2 != 0) result = mult(result,sq);
    };

    return result;
};


int Field::p_inverse(int x)
{
    int a=1, b=0, g=p, u=0, v=1, w=x;
    int q, t;

    if (w<0) w += p;
    while (w>0)
    {
        q=g/w;
        t=u; u=a-q*u; a=t;
        t=v; v=b-q*v; b=t;
        t=w; w=g-q*w; g=t;
    }
    if (b<0) b += p;
    return b % p;
}


std::vector<int> Field::mult(std::vector<int> a, std::vector<int> b)
{
    std::vector<int> result(2*f+1);

    for (int i=0; i<f; i++)
        for (int j=0; j<f; j++)
            result[i+j] = (result[i+j] + a[i]*b[j]) % p;
// reduce modulo poly
    for (int i=2*f; i>=f; i--){
        int factor = result[i];
        if (factor != 0){
            for (int j=0; j<=f; j++){
                int term = (result[i+j-f] - factor*poly[j]) % p;
                if (term < 0) term += p;
                result[i+j-f] = term;
            };
        };
    };

    return result;
}


void Field::create()
{
// calculate with polynomials of degree up to f
    poly[f] = 1;        // set leading term
    bool reducible;

    for (int n=0; n<q; n++)     // loop over all possible lower terms
    {
        int n0=n;
        for (int i=0; i<f; i++){
            poly[i] = n0 % p;
            n0 /= p;
        };

        std::vector<int> res(2*f+1); res[1]=1;  // res = X
        std::vector<int> test(2*f+1);

        reducible = false;

        for (int j=1; 2*j<=f; j++)  //  test for factor of degree j
        {
            res = power(res,p);
            test = res; test[1]--; if (test[1] < 0) test[1] = p-1;

            if ( gcd(test) == 0){
                reducible = true;
                break;
            };
        };

        if (reducible == false) {
            std::cout << "irreducible: ";
            for (int i=f; i>=0; i-- ) std::cout << poly[i] << " ";
            std::cout << "\n";
            break;
        };
    };
//    std::cout << "failed to find irreducible polynomial\n";
};


void Field::generator()
{
    std::vector<int> primes;
    gen.resize(f);

// determine prime factors of q-1
    int remainder = q-1;

    for (int i=2; i*i<=remainder; i++)
        if (remainder % i == 0){
            primes.push_back(i);
            while (remainder % i == 0) remainder /= i;
        };
    if (remainder > 1) primes.push_back(remainder);

    std::cout << "factors of " << q-1 << " are: ";
    for (auto x : primes) std::cout << x << ", ";
    std::cout << "\n";

    for (int n=1; n<q; n++)     // loop over all field elements
    {
        int n0=n;
        for (int i=0; i<f; i++){
            gen[i] = n0 % p;
            n0 /= p;
        };

        bool potential_generator = true;
    // calculate gen^((q-1)/primes[j]) mod poly
        for (auto x : primes){
            std::vector<int> res(2*f+1);
            res = power(gen, (q-1)/x);
    // if res == 1 report failure
            res[0]--; if (res[0]<0) res[0] = p-1;
            bool power_is_one = true;
            for (int k=0; k<f; k++) if ( res[k] != 0 )  { power_is_one = false; break; };
            if (power_is_one == false) continue;  // (q-1)/p_j -th power not = 1, proceed to next factor
    //  now (q-1)/p_j -th power is = 1, so n does not generate the multiplicative group of F
            potential_generator = false;
            break;      // no need to check other primes; continue with next higher n
        };
        if (potential_generator == true) {
            std::cout << "found generator: ";
            for (int i=f-1; i>=0; i-- ) std::cout << gen[i] << " ";
            std::cout << "\n";
            return;
        };
    };
    std::cout << "failed to find multiplicative generator\n";
};


void Field::tables()
{
    exp.resize(2*q);
    log.resize(q);
    std::vector<int> x = gen;
    negative.resize(q);

    // multiplication
    exp[0] = 1;
    for (int index=1; index<q; index++) {
        // find corresponding integer
        int n = 0; for (int i=f-1; i>=0; i--) n = n*p + x[i];
        int m = 0; for (int i=f-1; i>=0; i--) m = m*p + ( x[i]==0 ? 0 : p-x[i]);
        exp[index] = exp[index+q-1] = n;
        log[n] = index;
        x = mult(x,gen);
        negative[n] = m;
    };
    log[1]=0;

    // addition
    int p1 = 2*p-1; int q1=1;
    for (int i=0; i<f; i++) q1 *= p1;
    rebase.resize(q);
    unbase.resize(q1);

    std::vector<int> y(f);
    for (int index=0; index<q; index++) {
        int z = index;
        for (int i=0; i<f; i++){
            y[i] = z % p;
            z /= p;
        };
        int n = 0; for (int i=f-1; i>=0; i--) n = n * p1 + y[i];
        rebase[index] = n;
    };
    for (int index=0; index <q1; index++){
        int z = index;
        for (int i=0; i<f; i++){
            y[i] = z % p1;
            z /= p1;
        };
        int n = 0; for (int i=f-1; i>=0; i--) n = n * p + y[i] % p;
        unbase[index] = n;
    };
    std::cout << "all tables completed\n";

    return;
};


int Field::inverse(int x)
{
    return exp[q-1-log[x]];
};


int Field::neg(int x)
{
    return negative[x];
};


int Field::product(int x, int y)
{
    if (x == 0 || y == 0) return 0;
    return exp[log[x]+log[y]];
};


int Field::sum(int x, int y)
//int Field::sum0(int x, int y)
{
    return unbase[rebase[x]+rebase[y]];
};


//int Field::sum(int x, int y)
int Field::sum2(int x, int y)
{
    return x ^ y;
};


//int Field::sum(int x, int y)
int Field::sump(int x, int y)
{
    int z = x + y;
    return  (z >= p) ? z - p : z;
};
