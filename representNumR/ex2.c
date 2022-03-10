// https://moodle.c3sl.ufpr.br/course/view.php?id=493

#include <stdio.h>
#include <stdint.h>
#include <float.h>
#include <math.h> // compile with -lm

typedef union
{
    int32_t i;
    float f;
    // Bitfields for exploration. Do not use in production code.
    struct 
    { 
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    } parts;
} Float_t;


int bin2dec(long long n) {
    int dec = 0, i = 0, rem;

    while (n!=0) {
        rem = n % 10;
        n /= 10;
        dec += rem * pow(2, i);
        ++i;
    }

    return dec;
}

void printFloat_t( Float_t num )
{
    printf("%f -> f:%1.9e, ix:0x%08X, s:%d, e:%3d, ve:%3d mx:0x%06X\n",
    num.f, num.f, num.i,
    num.parts.sign, num.parts.exponent, num.parts.exponent - 127, num.parts.mantissa);

    
    // printf("\n\ntry: f:%1.9e, id:%ld, s:%ld, e:%ld, mx:%ld\n",
    // num.f, num.i,
    // num.parts.sign, num.parts.exponent, num.parts.mantissa);


}


int main()
{
    // printf("\nEpsilon: %1.15f\n\n", FLT_EPSILON);
    Float_t f;
    f.f = 0.0000001;
    printFloat_t(f);
    f.f = 0.0000002;
    printFloat_t(f);

    f.f = 2048.001;
    printFloat_t(f);
    f.f = 2048.002;
    printFloat_t(f);

    // f.f = 0.0;
    // printFloat_t(f);
    // f.f = 1.40129846e-45;
    // printFloat_t(f);
    // f.f = 1.17549435e-38;
    // printFloat_t(f);
    // f.f = 0.2;
    // printFloat_t(f);
    // f.f = 1.0;
    // printFloat_t(f);
    // f.f = 1.5;
    // printFloat_t(f);
    // f.f = 1.75;
    // printFloat_t(f);

}