#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <windows.h>
#include <stdlib.h>


#define PI 3.14159265358979323846


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

// #################################################


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

void _get_node_ab_index(int *idx_a, int *idx_b, int i, int k, int num_factors) {
    // step = num_factors = 1<<level
    // num_flies = N/2
    // a = i * (2*num_factors) + k = 2 * i * step + k
    // b = a + step
    int step = num_factors;
    *idx_a = 2 * i * step + k;
    *idx_b = *idx_a + step;
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

// ####################################################

void block_fft(Complex* input, Complex* output, int N) {
    for (size_t i=0; i < N; i ++) {
        output[i] = input[i];
    }
    fft(output, N);
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

typedef struct {
    Complex* middel;
    Complex* input;
    int block_size;
    int offset;
    int id;
} ThreadData;

DWORD WINAPI dft_thread(LPVOID arg) {
    ThreadData* data = (ThreadData*)arg;
    printf("start thread%d\n", data->id);
    DFT(data->middel+data->offset, data->input+data->offset, data->block_size);
    printf("finish thread%d\n", data->id);
    return 0;
}

DWORD WINAPI fft_thread(LPVOID arg) {
    ThreadData* data = (ThreadData*)arg;
    printf("start thread%d\n", data->id);
    block_fft(data->middel+data->offset, data->input+data->offset, data->block_size);
    printf("finish thread%d\n", data->id);
    return 0;
}


// parallel version for sequencial_dft, open the loop
void parallel_dft(Complex* input, Complex* middel, Complex* output, int N, int max_concurrency) {
    // check args
    assert(N % max_concurrency == 0);

    int block_size = N/max_concurrency;
    int num_blocks = max_concurrency;

    // perumerate all index to correct position
    _mod_N_perumerate(input, middel, block_size, num_blocks);

    // Create an array to hold thread handles
    HANDLE* threads = malloc(num_blocks * sizeof(HANDLE));
    ThreadData* thread_data = malloc(num_blocks * sizeof(ThreadData));

    // concurrently run blocks
    for (int block_idx = 0; block_idx < num_blocks; block_idx++) {
        thread_data[block_idx].id = block_idx;
        thread_data[block_idx].middel = middel;
        thread_data[block_idx].input = input;
        thread_data[block_idx].block_size = block_size;
        thread_data[block_idx].offset = block_idx * block_size;
        threads[block_idx] = CreateThread(NULL, 0, fft_thread, &thread_data[block_idx], 0, NULL);
    }

    // wait for all threads to finish
    WaitForMultipleObjects(num_blocks, threads, TRUE, INFINITE);
    printf("[info] all threads finish.\n");

    // close thread handles
    for (int block_idx = 0; block_idx < num_blocks; block_idx++) {
        CloseHandle(threads[block_idx]);
    }

    // clean up
    free(threads);
    free(thread_data);    

    printf("[info] all threads end. start merge.\n");
    // merge the last level
    merge(input, output, num_blocks, block_size, N);
    printf("[info] finish DFT\n");
}

// split original signal to from N to max_concurrency * N/max_concurrency, aka N = block_size * num_blocks
void sequencial_dft(Complex* input, Complex* middel, Complex* output, int N, int max_concurrency) {
    // check args
    assert(N % max_concurrency == 0);

    int block_size = N/max_concurrency;
    int num_blocks = max_concurrency;

    // perumerate all index to correct position
    _mod_N_perumerate(input, middel, block_size, num_blocks);

    // concurrently run blocks
    for (int block_idx=0; block_idx < num_blocks; block_idx ++) {
        int offset = block_idx*block_size;
        DFT(middel+offset, input+offset, block_size);
    }

    // merge the last level
    merge(input, output, num_blocks, block_size, N);
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

int test_sequencial_dft() {
// int main() {
    int N = 100;
    int max_concurrency = 25;
    Complex* input = (Complex*)malloc(sizeof(Complex)*N);
    Complex* middel = (Complex*)malloc(sizeof(Complex)*N);
    Complex* output = (Complex*)malloc(sizeof(Complex)*N);
    for (int i = 0; i < N; i++) {
        input[i].real = i;
        input[i].imag = 0;
    }

    printf("Input:\n");
    print_complex_array(input, N);

    sequencial_dft(input, middel, output, N, max_concurrency);

    printf("perumerate input(Middel):\n");
    print_complex_array(middel, N);
    printf("before merge(input):\n");
    print_complex_array(input, N);
    printf("DFT Output(output):\n");
    print_complex_array(output, N);
    return 0;
}

// int test_parallel_dft() {
int main() {
    int N = 4096*12;
    int max_concurrency = 12;
    Complex* input = (Complex*)malloc(sizeof(Complex)*N);
    Complex* middel = (Complex*)malloc(sizeof(Complex)*N);
    Complex* output = (Complex*)malloc(sizeof(Complex)*N);
    for (int i = 0; i < N; i++) {
        input[i].real = i;
        input[i].imag = 0;
    }

    parallel_dft(input, middel, output, N, max_concurrency);

    // printf("DFT Output(output):\n");
    // print_complex_array(output, N);
    return 0;
}