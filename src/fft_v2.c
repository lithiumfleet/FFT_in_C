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
    Complex res = {a->real*b->real - a->imag*b->imag, a->imag*b->real+a->real*b->imag};
    return res;
}

// Function to perform the bit reversal of the given index
unsigned int bit_reverse(unsigned int x, int log2n) {
    int n = 0;
    for (int i = 0; i < log2n; i++) {
        n <<= 1;
        n |= (x & 1);
        x >>= 1;
    }
    return n;
}

// Function to compute the FFT
void fft(Complex *x, int N) {
    int log2n = log2(N);

    // Bit-reversal permutation the input
    for (int i = 0; i < N; i++) {
        int rotate_index = bit_reverse(i, log2n);
        if (rotate_index > i) {
            Complex temp = x[i];
            x[i] = x[rotate_index];
            x[rotate_index] = temp;
        }
    }

    // Cooley-Tukey FFT
    for (int level = 1; level <= log2n; level++) {
        int num_rotate_factor = 1 << level; // num_rotate_factor = 2^level
        Complex w = {1.0, 0.0}; // Start with W_N^0
        Complex rotate_factor = {cos(PI / num_rotate_factor/2), -sin(PI / num_rotate_factor/2)}; // rotate_factor = e^(-2pi/num_rotate_factor)

        for (int rotate_index = 0; rotate_index < num_rotate_factor/2; rotate_index++) {
            
            // calulate one butterfly node
            for (int i = rotate_index; i < N; i += num_rotate_factor) {
                Complex t = multiply(&w, &x[i+num_rotate_factor/2]); // calculate C*B
                Complex u = x[i]; // calculate A
                x[i].real = u.real + t.real; x[i].imag = u.imag + t.imag; // A = A + C*B
                x[i + num_rotate_factor/2].real = u.real - t.real; x[i + num_rotate_factor/2].imag = u.imag - t.imag; // B = A - C*B
            }
            // printf("level:%d  rotate_index:%d/%d\n", level, rotate_index+1, num_rotate_factor/2);

            // Update w to the next value, or "rotate" the w with angle PI/num_rotate_factor/2
            w = multiply(&w, &rotate_factor);
        }
    }
}

// Function to print the complex array
void print_complex_array(Complex *x, int N) {
    for (int i = 0; i < N; i++) {
        printf("(%f, %f)\n", x[i].real, x[i].imag);
    }
}


//################//
//##### Test #####//
//################//


//##### Pass #####//
int main() {
    int N = 8; // Length of input array (must be power of 2)
    Complex x[N];
    for (size_t i=0; i < N; i ++) {
        Complex a = {i%10,(i*2+3)%20};
        x[i] = a;
    }
    printf("Input:\n");
    print_complex_array(x, N);
    fft(x, N);
    printf("\nOutput:\n");
    print_complex_array(x, N);
    return 0;
}
int test_bit_reverse() {
    int N = 8;
    for (size_t i=0; i < N; i ++) {
        printf("%zu -> %u\n", i, bit_reverse(i, log2(N)));
    }
    return 0;
}

//##### Fail #####//
