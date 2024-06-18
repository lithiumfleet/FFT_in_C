#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct {
    double real;
    double imag;
} Complex;

typedef struct {
    size_t n_rows;
    size_t n_cols;
    Complex** data;
} CMatrix;

Complex multiply(Complex a, Complex b) {
    Complex c;
    c.real = a.real*b.real - a.imag*b.imag;
    c.imag = a.real*b.imag + a.imag*b.real;
    return c;
}

Complex add(Complex a, Complex b) {
    Complex c;
    c.real = a.real + b.real;
    c.imag = a.imag + b.imag;
    return c;
}

void print_complex(Complex a) {
    printf("%.4f+%.4fj ", a.real, a.imag);
}

void print_cmatrix(CMatrix m) {
    if (!m.n_rows > 4 && m.n_cols > 4) {
    printf(\
"%.1f+%.1fj %.1f+%.1fj ... %.1f+%.1fj\n\
%.1f+%.1fj %.1f+%.1fj ... %.1f+%.1fj\n\
......\n\
%.1f+%.1fj %.1f+%.1fj ... %.1f+%.1fj\n\
%.1f+%.1fj %.1f+%.1fj ... %.1f+%.1fj\n", 
        m.data[0][0].real,                   m.data[0][0].imag,
        m.data[0][1].real,                   m.data[0][1].imag,
        m.data[0][m.n_cols-1].real,          m.data[0][m.n_cols-1].imag,
        m.data[1][0].real,                   m.data[1][0].imag,
        m.data[1][1].real,                   m.data[1][1].imag,
        m.data[1][m.n_cols-1].real,          m.data[1][m.n_cols-1].imag,
        m.data[m.n_rows-2][0].real,          m.data[m.n_rows-2][0].imag,
        m.data[m.n_rows-2][1].real,          m.data[m.n_rows-2][1].imag,
        m.data[m.n_rows-2][m.n_cols-1].real, m.data[m.n_rows-2][m.n_cols-1].imag,
        m.data[m.n_rows-1][0].real,          m.data[m.n_rows-1][0].imag,
        m.data[m.n_rows-1][1].real,          m.data[m.n_rows-1][1].imag,
        m.data[m.n_rows-1][m.n_cols-1].real, m.data[m.n_rows-1][m.n_cols-1].imag
        );
    }
    else {
        for (size_t i=0; i < m.n_rows; i ++) {
            for (size_t j = 0; j < m.n_cols; j++) {
                print_complex(m.data[i][j]);
            }
            printf("\n");
        }
    }
}

Complex new_complex(double real, double imag) {
    Complex res;
    res.real = real;
    res.imag = imag;
    return res;
}

CMatrix new_cmatrix(size_t n_rows, size_t n_cols) {

    CMatrix matrix;
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;

    matrix.data = (Complex**)malloc(n_rows * sizeof(Complex*));
    for (size_t i = 0; i < n_rows; i++) {
        matrix.data[i] = (Complex*)malloc(n_cols * sizeof(Complex));
    }

    for (size_t i = 0; i < n_rows; i ++) {
        for (size_t j = 0; j < n_cols; j ++) {
            matrix.data[i][j].imag = 0;
            matrix.data[i][j].real = 0;
        }
    }

    return matrix;
}

CMatrix matmul(CMatrix p, CMatrix q) {
    assert(p.n_cols == q.n_rows);
    CMatrix result = new_cmatrix(p.n_rows, q.n_cols);
    for (size_t i=0; i < p.n_rows; i ++) {
        for (size_t j=0; j < q.n_cols; j ++)  {
            Complex temp_sum = new_complex(0,0);
            for (size_t k=0; k < q.n_rows; k ++) {
                Complex product = multiply(p.data[i][k], q.data[k][j]);
                temp_sum = add(temp_sum, product);
            }
            result.data[i][j] = temp_sum;
        }
    }
    return result;
}

CMatrix new_R_2() {
    CMatrix R_2 = new_cmatrix(2,2);
    R_2.data[0][0] = new_complex(1,0);
    R_2.data[0][1] = new_complex(1,0);
    R_2.data[1][0] = new_complex(1,0);
    R_2.data[1][1] = new_complex(-1,0);
    return R_2;
}

Complex calculate_W(int k, int N) {
    // pure W, no special process
    int need_reverse = 0;
    if (k > N/2) {
        k %= N/2;
        need_reverse = 1;
    }

    double angle = 2 * PI * k / N; 
    Complex W_kN = new_complex(cos(angle), -sin(angle));

    if (need_reverse) { W_kN.imag *= -1; W_kN.real *= -1; }
    return W_kN;
}

Complex _get_identity_matrix_element() {
    return new_complex(1,0);
}

Complex _get_W_matrix_element(int i, int N) {
    Complex W_matrix_element;
    if (i < N/2) {
        W_matrix_element = calculate_W(i,N);
    } else {
        W_matrix_element = calculate_W(i-N/2,N);
        W_matrix_element.real *= -1;
        W_matrix_element.imag *= -1;
    }
    return W_matrix_element;
}

size_t _get_permutated_index(size_t j, int N) {
    return j%2==0 ? j/2 : (j-1)/2+N/2;
}

/// @brief return matrix R in shape (N, N)
/// @param R previous matrix R in shape (N/2, N/2)
/// @param N target N
/// @return matrix R in shape (N, N)
CMatrix FFT(CMatrix R, int N) {

    assert(R.n_cols==N/2 && R.n_rows==N/2);

    CMatrix result = new_cmatrix(N, N);

    for (int i=0; i < N; i ++) { // scan the lines
        // one line only contains two elements to calculate
        Complex identity_matrix_element = _get_identity_matrix_element();
        Complex W_matrix_element = _get_W_matrix_element(i, N);

        for (size_t j=0; j < N; j ++) {
            // permutation: 0,1,2,3,4,5,6,7 --> 0,2,4,6,1,3,5,7
            size_t permutation_j = _get_permutated_index(j, N);
            if ((i < N/2 && j < N/2) || (i >= N/2 && j < N/2)) {
                result.data[i][permutation_j] = multiply(identity_matrix_element, R.data[i%(N/2)][j%(N/2)]);
            } else {
                result.data[i][permutation_j] = multiply(W_matrix_element, R.data[i%(N/2)][j%(N/2)]);
            }
        }
    }
    return result;
}


/// @brief recursive version for fft
/// @param N target N
/// @return matrix R in shape (N, N)
CMatrix get_R(int N) {
    if (N == 2) return new_R_2();
    CMatrix R = FFT(get_R(N/2), N);
    return R;
}

CMatrix run_fft(CMatrix input, int N) {
    // FIXME: trunck X to N*1
    CMatrix R = get_R(N);
    CMatrix result = matmul(R, input);
    return result;
}

// Function to print complex numbers
void print_complex_array(Complex* arr, int N) {
    for (int i = 0; i < N; i++) {
        printf("Element %d: %.5f + %.5fi\n", i, arr[i].real, arr[i].imag);
    }
}

// ##########################
// ######### Tests ##########
// ##########################

// ######## Passed ##########
void test_rotate_factor() {
    int N = 8;
    for (size_t i = 0; i < N; i++) print_complex(calculate_W(i, N));
}

void test_new_R_2() {
    print_cmatrix(new_R_2());
}

void test_multiply() {
    Complex a = new_complex(1,2);
    Complex b = new_complex(3,-1);
    printf("\na:\n");
    print_complex(a);
    printf("\nb:\n");
    print_complex(b);
    printf("\nres:\n");
    print_complex(multiply(a,b));
}

void test_matmul() {
    CMatrix a = new_cmatrix(2,3);
    CMatrix b = new_cmatrix(3,1);
    a.data[0][0] = new_complex(2,0);
    a.data[1][1] = new_complex(4,0);
    a.data[1][2] = new_complex(1,0);
    b.data[0][0] = new_complex(1,0);
    b.data[2][0] = new_complex(2,0);
    printf("matrix a\n");
    print_cmatrix(a);
    printf("matrix b\n");
    print_cmatrix(b);
    printf("result \n");
    print_cmatrix(matmul(a,b));
}

void test_permutation() {
    int a[] = {0,1,2,3,4,5,6,7};
    for (size_t i=0; i < 8; i ++) {
        printf("%d -> %zu\n", a[i], _get_permutated_index(i, 8));
    }
}

// ####### Failed ###########

void test_get_R() {
    CMatrix R = get_R(8);
    print_cmatrix(R); // ATTENTION: print_cm is changed for debug, after it please remove the ! in if expr!!!!!!!!!!
}

void test_fft() {
    int N = 8; // Number of points in the DFT
    CMatrix input = new_cmatrix(N, 1);
    for (int i = 0; i < N; i++) {
        input.data[i][0].real = i + 1;
        input.data[i][0].imag = 0;
    }
    printf("input:\n");
    print_cmatrix(input);
    CMatrix output = run_fft(input, N);
    printf("FFT Output:\n");
    print_cmatrix(output);
}

int main() {
    test_get_R();
    return 0;
}