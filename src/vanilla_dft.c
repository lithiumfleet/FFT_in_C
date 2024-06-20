#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

// Structure to hold complex numbers
typedef struct {
    double real;
    double imag;
} Complex;

// Function to compute the DFT
void DFT(Complex* input, Complex* output, int N) {
    for (int k = 0; k < N; k++) { // For each output element
        output[k].real = 0;
        output[k].imag = 0;
        for (int n = 0; n < N; n++) { // For each input element
            double angle = 2 * PI * k * n / N;
            output[k].real += input[n].real * cos(angle) + input[n].imag * sin(angle);
            output[k].imag += input[n].imag * cos(angle) - input[n].real * sin(angle);
        }
    }
}

// Function to print complex numbers
void printComplexArray(Complex* arr, int N) {
    for (int i = 0; i < N; i++) {
        printf("Element %d: %.5f + %.5fi\n", i, arr[i].real, arr[i].imag);
    }
}

void run_dft() {
    int N = 7; // Number of points in the DFT
    Complex input[N];
    Complex output[N];
    for (int i = 0; i < N; i++) {
        input[i].real = i;
        input[i].imag = 0;
    }

    printf("Input:\n");
    printComplexArray(input, N);

    DFT(input, output, N);

    // Print the output
    printf("DFT Output:\n");
    printComplexArray(output, N);
}

int main() {
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    run_dft();

    end = clock();
    clock_t ticks = end - start;
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("total ticks:%zu\n", ticks);
    printf("total cpu time:%.10f", cpu_time_used);
    return 0;
}

