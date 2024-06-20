#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

typedef struct {
    double real;
    double imag;
} Complex;

Complex multiply(Complex* a, Complex* b) {
    Complex res = {a->real*b->real - a->imag*b->imag, a->imag*b->real + a->real*b->imag};
    return res;
}

// function to perform the bit reversal of the given index
size_t _permurate_index(size_t x, int N) {
    int n = 0;
    int log2n = log2(N);
    for (int i = 0; i < log2n; i++) {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n;
}

void permurate(Complex* input, int N) {
    for (int i = 0; i < N; i++) {
        int j = _permurate_index(i, N);
        if (j > i) { Complex temp = input[i]; input[i] = input[j]; input[j] = temp; }
    }
}

/// @brief get W_N^k vector
/// @param W Complex vector.
/// @param num_factors the length of W
/// @param N the N point
void _get_rotate_factors(Complex* W, size_t num_factors, int N) {
    Complex base_rotate_factor = {1.0, 0.0}; // start with W_N^0
    Complex rotate_angle = {cos(PI / num_factors), -sin(PI / num_factors)}; // e^(-2pi/m)

    W[0] = base_rotate_factor;
    for (size_t i=1; i < num_factors; i ++) {
        W[i] = multiply(&W[i-1], &rotate_angle);
    }
}

void print_complex_array(Complex *x, int N) {
    for (int i = 0; i < N; i++) {
        printf("(%.2f, %.2f)\n", x[i].real, x[i].imag);
    }
}

void _get_node_ab_index(int *idx_a, int *idx_b, int i, int k, int num_factors) {
    // step = num_factors = 1<<level
    // num_flies = N/2
    // a = i * (2*num_factors) + k = 2 * i * step + k
    // b = a + step
    int step = num_factors;
    *idx_a = 2 * i * step + k;
    *idx_b = *idx_a + step;
}

inline void add(Complex* output, Complex* x) {
    output->real += x->real;
    output->imag += x->imag;
}
inline void minus(Complex* output, Complex* x, Complex* y) {
    output->real = x->real - y->real;
    output->imag = x->imag - y->imag;
}

void fft(Complex* input, int N) {
    // inplace version

    permurate(input, N);

    int num_flies = N/2; // number of butterflies in each level
    for (size_t level=0; level < log2(N); level ++) {

        // 1. prepare W
        int num_factors = 1<<level;
        Complex* W = (Complex*)malloc(sizeof(Complex)*num_factors);
        _get_rotate_factors(W, num_factors, N); // then we can use W[k] as W_N^k

        // 2. through the rotate factors, and calculate ABC
        for (size_t k=0; k < num_factors; k ++) {
            // current rotate factor is W[k], and W[k] will be used for num_flies/num_factors times.
            for (size_t i=0; i < num_flies/num_factors; i ++) {
                int idx_a, idx_b;
                _get_node_ab_index(&idx_a, &idx_b, i, k, num_factors);
                Complex BC = multiply(&input[idx_b], &W[k]);
                minus(&input[idx_b], &input[idx_a], &BC);
                add(&input[idx_a], &BC);
            }
        }
    }
}


//################//
//##### Test #####//
//################//

//##### Pass #####//
int test_peruerate_index() {
    int N = 8;
    for (size_t i=0; i < N; i ++) {
        printf("%zu -> %zu\n", i, _permurate_index(i, N));
    }
    return 0;
}
int test_rotate_factors() {
    int level = 3; // level = 3+1
    int N = 16;
    int num_factors = 1<<level;
    Complex W[num_factors];
    _get_rotate_factors(W, num_factors, N); // then we can use W[k] as W_N^k
    print_complex_array(W, N/2);
    return 0;
}
int test_get_node_ab_index() {
    // i k l --> a b step
    // 0 0 0     0 1  1
    // 1 0 0     2 3  1
    // 2 0 0     4 5  1
    // 3 0 0     6 7  1
    //               
    // 0 0 1     0 2  2
    // 0 1 1     1 3  2
    // 1 0 1     4 6  2
    // 1 1 1     5 7  2
    //               
    // 0 0 2     0 4  4
    // 0 1 2     1 5  4
    // 0 2 2     2 6  4
    // 0 3 2     3 7  4
    int I[] = {0,1,2,3,0,0,1,1,0,0,0,0};
    int K[] = {0,0,0,0,0,1,0,1,0,1,2,3};
    int L[] = {0,0,0,0,1,1,1,1,2,2,2,2};
    printf("i k l --> a b step\n");
    for (size_t j=0; j < 12; j ++) {
        int a, b;
        _get_node_ab_index(&a,&b,I[j], K[j], L[j]);
        printf("%d %d %d     %d %d  %d\n", I[j], K[j], L[j], a, b, 1<<L[j]);
    }
    return 0;
}
void test_fft() {
    const int test_N = 1000000;
    Complex* test_input = (Complex*)malloc(sizeof(Complex)*test_N); // so... a stack can't hold this...
    for (int i = 0; i < test_N; i++) {
        test_input[i].real = i;
        test_input[i].imag = 0;
    }

    // printf("Input:\n");
    // print_complex_array(x, N);

    fft(test_input, test_N);

    // printf("\nOutput:\n");
    // print_complex_array(x, N);
    free(test_input);
}

//##### Fail #####//


//##### main #####//
// fft usage:
// 1. create a signal sequence. all elements should be Complex
// 2. define the number of points in FFT
// 3. use  print_complex_array to print input and output.

int main() {

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    test_fft();

    end = clock();
    clock_t ticks = end - start;
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("total ticks:%zu\n", ticks);
    printf("total cpu time:%.10f", cpu_time_used);

    return 0;
}