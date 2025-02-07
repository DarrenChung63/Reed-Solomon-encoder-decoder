#include <stdio.h>
#include <stdlib.h>

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

// Compute generator polynomial g(x)
void compute_generator_polynomial(int g[]) {
    int roots[T] = {3, 9, 27, 81}; // α = 3, α^2 = 9, α^3 = 27, α^4 = 81 in GF(929)
    g[T] = -3; 
    g[T - 1] = 1;

    int tempPoly1[T + 1] = {0};
    int tempPoly2[T + 1] = {0};
    for (int i = 1; i < T; i++) {
        for (int j = T; j > 0; j--) {
            tempPoly1[j - 1] = g[j];
        }
        for (int j = 0; j < T + 1; j++) {
            tempPoly2[j] = gf_mult(-roots[i], g[j]);
        }
        for (int j = 0; j < T + 1; j++) {
            g[j] = mod(tempPoly1[j] + tempPoly2[j], PRIME);
        }
    }
}

// Encode message polynomial p(x)
void encode_rs(int message[], int codeword[]) {
    int g[T + 1] = {0};
    compute_generator_polynomial(g);
    int temp[CODE_LEN] = {0};
    int s_r[CODE_LEN] = {0};
    
    // Copy message to temp
    for (int i = 0; i < MESSAGE_LEN; i++) {
        temp[i] = message[i];
        s_r[i] = message[i];
    }

    // Compute remainder
    for (int i = 0; i < CODE_LEN - T; i++) {
        int coef = s_r[i] / g[0];
        for (int j = 0; j < T + 1; j++) {
            s_r[i + j] = mod(s_r[i + j] - coef * g[j], PRIME);
        }
    }

    // Form the final codeword
    for (int i = 0; i < CODE_LEN; i++) {
        codeword[i] = mod(temp[i] - s_r[i], PRIME);
    }
}

// Main function to test the encoder
int main() {
    int message[MESSAGE_LEN] = {3, 2, 1};  // p(x) = 3x^2 + 2x + 1
    int codeword[CODE_LEN] = {0};

    encode_rs(message, codeword);

    printf("Encoded Codeword: ");
    for (int i = 0; i < CODE_LEN; i++) {
        printf("%d ", codeword[i]);
    }
    printf("\n");


    return 0;
}
