# Reed-Solomon Encoder and Decoder

This repository contains a simple implementation of a Reed-Solomon encoder and decoder in C. The implementation is based on a finite field GF(929) and can correct up to 4 errors.

## Features
- **Reed-Solomon Encoder**: Encodes a given message polynomial into a codeword.
- **Reed-Solomon Decoder**: Detects and corrects errors using syndrome computation, error locator polynomial, and the Forney algorithm.

## Files
- `Reed-Solomon-encoder.c` - Implements the encoding of a message into a Reed-Solomon codeword.
- `Reed-Solomon-decoder.c` - Implements error detection and correction for a received codeword.

## Reed-Solomon Parameters
- **Finite Field**: GF(929)
- **Code Length (n)**: 7
- **Message Length (k)**: 3
- **Error Correction Capability (t)**: 4

## Compilation & Execution
To compile the encoder and decoder, use the following commands:

```sh
# Compile the encoder
gcc -o encoder Reed-Solomon-encoder.c

# Compile the decoder
gcc -o decoder Reed-Solomon-decoder.c
```

### Running the Encoder
```sh
./encoder
```
This will generate an encoded Reed-Solomon codeword.

### Running the Decoder
```sh
./decoder
```
This will correct errors in the received codeword and print the corrected output.

## How It Works
### Encoding Process
1. Defines a generator polynomial.
2. Computes the parity symbols.
3. Constructs the final codeword.

### Decoding Process
1. Computes syndromes.
2. Uses the Berlekamp-Massey algorithm to find the error locator polynomial.
3. Identifies error locations.
4. Corrects the errors using the Forney algorithm.

## Example
### Encoding Output
```
Encoded Codeword: 3 2 1 382 191 487 474
```

### Decoding Output
```
Corrected Codeword: 3 2 1 382 191 487 474
```

## License
This project is open-source and available for educational purposes.
