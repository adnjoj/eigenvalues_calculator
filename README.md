# QR Eigenvalue Solver with Double Shift

This program computes all eigenvalues of a real square matrix using the QR algorithm with an implicit double shift QR algorithm. It first reduces the matrix to upper Hessenberg form and then applies the shifted QR iteration, efficiently handling both real eigenvalues and complex conjugate pairs.

## Input Format

Path to the input file is passed as the first argument. The input must follow this structure:

- First line: a single integer `n` – the size of the square matrix.
- Next `n` lines: each contains `n` space‑separated real numbers representing a row of the matrix.

Example for a 3×3 matrix:

```
3
4 1 0
1 3 2
0 2 5
```


## Output Format

The path to the output file is passed as the second argument. The computed eigenvalues are printed one per line. Each eigenvalue is given as a complex number in the form:

- For real eigenvalues: `real` (e.g., `2.71828`)
- For complex eigenvalues: `real +imag i` (e.g., `1.0 -2.5i`, `0 +1i`)

The order of eigenvalues is arbitrary (non‑sorted).

Example output for the above matrix:

```
6
3
3
```
