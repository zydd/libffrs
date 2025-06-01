// #include <x86intrin.h>
#include <immintrin.h>
#include <stdio.h>

int main() {
    // {
    //     __m128i a = {0x0101010101010102, 0x010102};
    //     __m128i b = {0x3f3f3f3f3f3f3f, 0x3f3f3f3f3f3f3f};

    //     int expect1 = a[1] * b[1];
    //     int expect2 = a[0] * b[0];
    //     int expect3 = a[0] * b[1];

    //     // __m128i result1 = _mm_clmulepi64_si128(a, b, 0xff);

    //     // _mm_gf2p8affine_epi64_epi8(0,0,0);
    //     __m128i result1 = _mm_gf2p8mul_epi8(a, b);

    //     printf("ex1: %x\n", expect1);

    //     printf("a: %llx %llx\n", a[0], a[1]);
    //     printf("b: %llx %llx\n", b[0], b[1]);
    //     printf("rs1: %llx %llx\n", result1[0], result1[1]);
    // }
    printf("\n");
    {
        typedef int32_t i4x4 __attribute__ ((vector_size(16)));
        typedef int32_t i4x8 __attribute__ ((vector_size(32)));
        typedef int32_t i4x16 __attribute__ ((vector_size(64)));
        typedef int32_t i4x32 __attribute__ ((vector_size(128)));
        // union {
        //     i4x4 v;
        //     int32_t i4[4];
        // } a;
        i4x16 a = {1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1};
        i4x16 b = {1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,16};
        a <<= b;
        printf("a: %d %d %d %d  %d %d %d %d  %d %d %d %d  %d %d %d %d\n",
               a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7],
               a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15]);

    }

    return 0;
}
