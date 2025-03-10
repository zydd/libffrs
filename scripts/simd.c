// #include <x86intrin.h>
#include <immintrin.h>
#include <stdio.h>

int main() {

    __m128i a = {0x0101010101010102, 0x010102};
    __m128i b = {0x3f3f3f3f3f3f3f, 0x3f3f3f3f3f3f3f};

    int expect1 = a[1] * b[1];
    int expect2 = a[0] * b[0];
    int expect3 = a[0] * b[1];

    // __m128i result1 = _mm_clmulepi64_si128(a, b, 0xff);

    // _mm_gf2p8affine_epi64_epi8(0,0,0);
    __m128i result1 = _mm_gf2p8mul_epi8(a, b);

    printf("ex1: %x\n", expect1);

    printf("a: %llx %llx\n", a[0], a[1]);
    printf("b: %llx %llx\n", b[0], b[1]);
    printf("rs1: %llx %llx\n", result1[0], result1[1]);

    return 0;
}
