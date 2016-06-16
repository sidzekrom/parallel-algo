#include <bits/stdc++.h>
using namespace std;

/* Notes on the Discrete Fourier Transform */

/* the following is an implementation of the discrete fourier transform,
   that is, the Fourier Transform of functions from {-1, 1} --> R where R
   are the real numbers. */

/* Let A be a matrix that is 2^n by 2^n.
   index the rows and columns by vectors from F_2^n.
   And let column x contain the Fourier characteristic function
   for vector x, denoted chi_X. and row y of the column is the value
   chi_X takes on y. A[x,y] = (-1)^{x.y} */

/* The Fourier Transform of a function f from F_2^n -> R is given by
   Af, where f is denoted as a 2^n dimensional vector */

/* Naively computing Af can gives a runtime of O(2^{2n}) or O(m^2), where
   m is the size of the field we are working over. By expoiting properties
   of A, we can speed this up to O(m log m). And do even faster in parallel
   (parallel span is O(log m) */

/* we write a single threaded version of the DFT first */

