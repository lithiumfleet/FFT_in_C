#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define PI 3.14159265358979323846

// structure to hold complex numbers
typedef struct {
    double real;
    double imag;
} Complex;

Complex multiply(Complex* a, Complex* b) {
    Complex c;
    c.real = a->real*b->real - a->imag*b->imag;
    c.imag = a->real*b->imag + a->imag*b->real;
    return c;
}

inline void add(Complex* a, Complex* b) {
    a->real += b->real;
    a->imag += b->imag;
}

// function to compute the DFT
void DFT(Complex* input, Complex* output, int N) {
    for (int k = 0; k < N; k++) { // for each output element
        output[k].real = 0;
        output[k].imag = 0;
        for (int n = 0; n < N; n++) { // for each input element
            double angle = 2 * PI * k * n / N;
            output[k].real += input[n].real * cos(angle) + input[n].imag * sin(angle);
            output[k].imag += input[n].imag * cos(angle) - input[n].real * sin(angle);
        }
    }
}

// function to print complex numbers
void print_complex_array(Complex* arr, int N) {
    for (int i = 0; i < N; i++) {
        printf("Element %d: %.5f + %.5fi\n", i, arr[i].real, arr[i].imag);
    }
}
void print_complex(Complex a) {
    printf("%.5f + %.5fi\n", a.real, a.imag);
}

// return where the index i should be
inline size_t _get_correct_index(size_t i, int block_size, int num_blocks) {
    int p = num_blocks;
    int r = i%block_size;
    int l = i/block_size;
    return p*r+l;
}

// auxiliary function to reorder the index of input.
// FIXME: need to simplify
void _mod_N_perumerate(Complex* input, Complex* output, int block_size, int num_blocks) {
    size_t index = 0;
    size_t base = 0;
    for (size_t i=0; i < block_size; i ++) {
        for (size_t j=0; j < num_blocks; j ++) {
            output[base+j*block_size] = input[index];
            index ++;
        }
        base ++;
    }
}

// lazy... this function from dft_v1
Complex calculate_W(int k, int N) {
    // pure W, no special process
    double angle = 2 * PI * k / N; 
    Complex W_kN;
    W_kN.real = cos(angle);
    W_kN.imag = -sin(angle);
    return W_kN;
}

// function to merge the last level
void merge(Complex* input, Complex* output, int num_blocks, int block_size, int N) {
    for (int k = 0; k < N; k++) { // for each output element
        output[k].real = 0; output[k].imag = 0;
        // to get output[k] need 
        // {... input[k-block_size], input[k], input[k+block_size], ...}.
        // and an array of W_N^k
        // this is what 'i' go through
        for (int i = 0; i < num_blocks; i ++) {
            size_t index = k%block_size + i*block_size;
            Complex w = calculate_W(i*k,N);
            Complex temp_production = multiply(&input[index], &w);
            add(&output[k], &temp_production);
        }
    }
}

// parallel dft: split original signal to from N to max_concurrency * N/max_concurrency, aka N = block_size * num_blocks
void parallel_dft(Complex* input, Complex* middel, Complex* output, int N, int max_concurrency) {
    // check args
    assert(N % max_concurrency == 0);

    int block_size = N/max_concurrency;
    int num_bocks = max_concurrency;

    // perumerate all index to correct position
    _mod_N_perumerate(input, middel, block_size, num_bocks);

    // concurrently run blocks
    for (int block_idx=0; block_idx < num_bocks; block_idx ++) {
        int offset = block_idx*block_size;
        DFT(middel+offset, input+offset, block_size);
    }

    // merge the last level
    merge(input, output, num_bocks, block_size, N);
}


//################
//##### Test #####
//################
int test_get_correct_index_1() {
// int main() {
    // when 3 blocks in size 2
    // i --> j
    // 0     0
    // 1     3
    // 2     1
    // 3     4
    // 4     2
    // 5     5
    int N = 6;
    int num_blocks = 3;
    int block_size = 2;
    for (size_t i=0; i < N; i ++) {
        size_t j = _get_correct_index(i, block_size, num_blocks);
        printf("%zu -> %zu\n", i, j);
    }
    return 0;
}

int test_get_correct_index_2() {
// int main() {
    // when 2 blocks in size 3
    // i --> j
    // 0     0
    // 1     2
    // 2     4
    // 3     1
    // 4     3
    // 5     5
    int N = 6;
    int num_blocks = 2;
    int block_size = 3;
    for (size_t i=0; i < N; i ++) {
        size_t j = _get_correct_index(i, block_size, num_blocks);
        printf("%zu -> %zu\n", i, j);
    }
    return 0;
}

int test_split_chunck_no_cross() {
// int main() {

    int N = 18;
    int max_concurrency = 3;
    int block_size = N/max_concurrency;
    int num_bocks = N/block_size;

    Complex input[N];
    Complex output[N];
    for (int i = 0; i < N; i++) {
        input[i].real = i;
        input[i].imag = 0;
    }

    printf("Input:\n");
    print_complex_array(input, N);
    // concurrently run blocks
    for (int block_idx=0; block_idx < num_bocks; block_idx ++) {
        int offset = block_idx*block_size;
        DFT(input+offset, output+offset, block_size);
    }

    // Print the output
    printf("chunck DFT Output:\n");
    print_complex_array(output, N);
    return 0;
}

int run_dft() {
// int main() {
    int N = 18;
    int max_concurrency = 3;
    Complex input[N];
    Complex middel[N];
    Complex output[N];
    for (int i = 0; i < N; i++) {
        input[i].real = i;
        input[i].imag = 0;
    }

    printf("Input:\n");
    print_complex_array(input, N);

    parallel_dft(input, middel, output, N, max_concurrency);

    printf("perumerate input(Middel):\n");
    print_complex_array(middel, N);
    printf("before merge(input):\n");
    print_complex_array(input, N);
    printf("DFT Output(output):\n");
    print_complex_array(output, N);
    return 0;
}