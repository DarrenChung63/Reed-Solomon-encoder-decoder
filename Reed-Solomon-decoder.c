#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#define PRIME 929  // Finite field GF(929)
#define T 4        // Error correction capability
#define MESSAGE_LEN 3
#define CODE_LEN 7  // RS(7,3)

// Modular operations in GF(929)
int mod(int a, int p) {
    return ((a % p) + p) % p;
}

// Multiplication in GF(929)
int gf_mult(int a, int b) {
    return mod(a * b, PRIME);
}

// Power in GF(929)
int gf_power(int a, int b) {
    int result = 1;
    if (b == 0) {
        return 1;
    }
    else if (b > 1) {
        for (int i = 0; i < b; i++) {
            result = gf_mult(result, a);
        }
    }
    else {
        for (int i = 0; i < PRIME + b - 1; i++) {
            result = gf_mult(result, a);
        }
    }

    return result;
}

// Division in GF(929)
int gf_div(int a, int b) {
    return mod(a * gf_power(b, PRIME - 2), PRIME);
}

// Function to compute syndromes
int* compute_syndromes(int alpha, int* received_code, int error_num) {
    int* syndromes = (int*)malloc(2 * error_num * sizeof(int));
    memset(syndromes, 0, 2 * error_num * sizeof(int));
    for (int i = 0; i < 2 * error_num; i++) {
        for (int j = 0; j < CODE_LEN; j++) {
            syndromes[i] = mod(syndromes[i] + gf_mult(received_code[j], gf_power(alpha, (i + 1) * (CODE_LEN - 1 - j))), PRIME);
        }
    }

    return syndromes;
}

// Function to comput error locator polynomial coefficients
int* compute_error_locator(int* syndromes, int error_num) {
    // create syndrome matrix
    int syndrome_mat[error_num][error_num + 1] = {0};
    for (int i = 0; i < error_num; i++) {
        for (int j = 0; j < error_num; j++) {
            syndrome_mat[i][j] = syndromes[i + j];
        }
    }
    for (int i = 0; i < error_num; i++) {
        syndrome_mat[i][error_num] = mod(-syndromes[error_num + i], PRIME);
    }

    // Gaussian elimination
    for (int i = 0; i < error_num; i++) {
        // Pivoting: Make the diagonal element 1
        int pivot = syndrome_mat[i][i];
        for (int j = i; j < error_num + 1; j++) {
            if (pivot != 0) {
                syndrome_mat[i][j] = gf_div(syndrome_mat[i][j], pivot);
            }
        }

        // Elimination: Make the element below the diagonal 0
        for (int k = i + 1; k < error_num; k++) {
            int factor = syndrome_mat[k][i];
            for (int j = i; j < error_num + 1; j++) {
                syndrome_mat[k][j] = mod(syndrome_mat[k][j] + gf_mult(syndrome_mat[i][j], -factor), PRIME);
            }
        }
    }

    // Back substitution
    for (int i = error_num - 1; i > 0; i--) {
        for (int k = i - 1; k >= 0; k--) {
            int factor = syndrome_mat[k][i];
            for (int j = i; j < error_num + 1; j++) {
                syndrome_mat[k][j] = mod(gf_mult(-factor, syndrome_mat[i][j]) + syndrome_mat[k][j], PRIME);
            }
        }
    }

    int* error_locator = (int*)malloc((error_num + 1) * sizeof(int));
    memset(error_locator, 0, (error_num + 1) * sizeof(int));
    for (int i = 0; i < error_num; i++) {
        error_locator[i] = syndrome_mat[i][error_num];
    }
    error_locator[error_num] = 1;

    return error_locator;
}

// Calculate differential of a polynomial
int* gf_diff(int* poly, int poly_order) {
    int* diff_poly = (int*)malloc((poly_order - 1) * sizeof(int));
    memset(diff_poly, 0, (poly_order - 1) * sizeof(int));
    for (int i = 0; i < poly_order - 1; i++) {
        diff_poly[i] = gf_mult(poly_order - 1 - i, poly[i]);
    }

    return diff_poly;
}

// Forney Algorithm
int* forney_alg(int* syndrome, int* error_locator, int* error_positions, int error_num, int alpha) {
    // Compute error evaluation polynomial
    int* omega_temp = (int*)malloc(3 * error_num * sizeof(int));
    memset(omega_temp, 0, 3 * error_num * sizeof(int));
    for (int i = 0; i < 2 * error_num; i++) {
        for (int j = 0; j < error_num + 1; j++) {
            omega_temp[i + j] = mod(omega_temp[i + j] + gf_mult(syndrome[2 * error_num - i - 1], error_locator[j]), PRIME);
        }
    }

    int* omega = (int*)malloc(T * sizeof(int));
    memset(omega, 0, T * sizeof(int));
    for (int i = 3 * error_num - T; i < 3 * error_num; i++) {
        omega[i - 3 * error_num + T] = omega_temp[i];
    }

    // compute differential of error locator polynomial
    int* lambda_diff = (int*)malloc(error_num * sizeof(int));
    memset(lambda_diff, 0, error_num * sizeof(int));

    lambda_diff = gf_diff(error_locator, error_num + 1);

    // compute the error values
    int* error_values = (int*)malloc(CODE_LEN * sizeof(int));
    memset(error_values, 0, CODE_LEN * sizeof(int));
    int alpha_inv = gf_power(alpha, -1);
    for (int i = 0; i < error_num; i++) {
        int omega_val = 0;
        for (int j = 0; j < T; j++) {
            omega_val = mod(omega_val + gf_mult(omega[j], gf_power(gf_power(alpha_inv, error_positions[i]), T - 1 - j)), PRIME);
        }
        int lambda_diff_val = 0;
        for (int j = 0; j < error_num; j++) {
            lambda_diff_val = mod(lambda_diff_val + gf_mult(lambda_diff[j], gf_power(gf_power(alpha_inv, error_positions[i]), error_num - 1 - j)), PRIME);
        }

        if (lambda_diff_val != 0) {
            error_values[CODE_LEN - error_positions[i] - 1] = gf_div(-omega_val, lambda_diff_val);
        }
    }

    return error_values;
}

int main() {
    int alpha = 3;
    
    int received_code[CODE_LEN] = {3, 2, 123, 456, 191, 487, 474};
    int error_num = 2;
    int error_positions[error_num] = {3, 4}; // or some other values
    
    // compute syndromes
    int* syndromes = compute_syndromes(alpha, received_code, error_num);

    // compute error locator polynomial coefficients
    int* error_locator = compute_error_locator(syndromes, error_num);
    
    // compute error values using Forney Algorithm
    int* error_values = forney_alg(syndromes, error_locator, error_positions, error_num, alpha);
    // for (int i = 0; i < CODE_LEN; i++) {
    //     printf("%d ", error_values[i]);
    // }

    // compute the correct code word
    int* corrected_code = (int*)malloc(CODE_LEN * sizeof(int));
    memset(corrected_code, 0, CODE_LEN * sizeof(int));
    for (int i = 0; i < CODE_LEN; i++) {
        corrected_code[i] = received_code[i] - error_values[i];
        printf("%d ", corrected_code[i]);
    }

    return 0;
}
