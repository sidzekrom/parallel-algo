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
/* the reason we write this is because the multithreaded version is actually
   slower for small input sizes! So once we are down to a case where the matrix
   is 16x16, we apply the sequential version instead of staying multithreaded
   until we are all the way down to the base case: the overhead is too much
   then */
vector<long double>
    dft_sequential(vector<long double> func, int start_ind, int end_ind,
        long double normalization) {

    /* just because of the nature of the Fourier transform, we require that
     * end_ind - start_ind be a power of 2 */
    if (start_ind == end_ind) {
        /* this is our base case where we want the Fourier Transform of a
         * function on one variable (the zero dimensional case) */
        return {func[start_ind]};
    }

    /* it is an invariant in the F_2 case that end_ind - start_ind is a
     * power of 2. We maintain it */
    int range = end_ind - start_ind;
    int mid_ind = start_ind + (range / 2);

    /* it is easy to see the power of 2 invariant is maintained */
    /* the following is the 'divide' step */
    vector<long double> dft_top = dft_sequential(func, start_ind, mid_ind);
    vector<long double> dft_bottom = dft_sequential(func, mid_ind, end_ind);
    
    /* the following is the conquer step */
    vector<long double> discrete_ft(end_ind - start_ind);

    /* make sure the scaling factors are correct */
    for (int i = 0; i < mid_ind - start_ind; i++) {
        discrete_ft[i] = normalization * (dft_top[i] + dft_bottom[i]);
        discrete_ft[mid_ind - start_ind + i] = normalization * (dft_top[i] -
            dft_bottom[i]);
    }

    return discrete_ft;
}

int main() {
    return 0;
}
