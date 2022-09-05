#include "bin_interp.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * @brief Convert an integer into a binary array
 * 
 * @param in 
 * @param count 
 * @param out 
 */
void int_to_bin_digit(int num, int count, int* out) {
    /* assert: count <= sizeof(int)*CHAR_BIT */
    unsigned int mask = 1U << (count-1);
    int i;
    for (i = 0; i < count; i++) {
        out[i] = (num & mask) ? 1 : 0;
        num <<= 1;
    }
}

/**
 * @brief Find minimum element in array. Returns double
 * 
 * @param vec 
 * @param count 
 * @return double 
 */
double vector_min_double(double *vec, int count) {
    double min_val = vec[0];
    for (int i = 0; i < count; i++) {
        if (vec[i] < min_val) min_val = vec[i];
    }
    return min_val;
}

/**
 * @brief Find minimum element in array. Returns int
 * 
 * @param vec 
 * @param count 
 * @return int 
 */
int vector_min_int(int *vec, int count) {
    int min_val = vec[0];
    for (int i = 0; i < count; i++) {
        if (vec[i] < min_val) min_val = vec[i];
    }
    return min_val;
}

/**
 * @brief Add two double vectors, u and v, of same length and type count into out.
 * 
 * @param u 
 * @param v 
 * @param count 
 * @param out 
 */
void vector_add_double(double *u, double *v, int count, double *out) {
    for (int i = 0; i < count; i++) {
        out[i] = u[i] + v[i];
    }
}

/**
 * @brief Add two integer vectors, u and v, of same length and type count into out.
 * 
 * @param u 
 * @param v 
 * @param count 
 * @param out 
 */
void vector_add_int(int *u, int *v, int count, int *out) {
    for (int i = 0; i < count; i++) {
        out[i] = u[i] + v[i];
    }
}

/**
 * @brief Multiply two double vectors elementwise of same length into out.
 * 
 * @param u 
 * @param v 
 * @param count 
 * @param out 
 */
void vector_mult_double(double *u, double *v, int count, double *out) {
    for (int i = 0; i < count; i++) {
        out[i] = u[i] * v[i];
    }
}

/**
 * @brief Multiply two integer vectors elementwise of same length and into out.
 * 
 * @param u 
 * @param v 
 * @param count 
 * @param out 
 */
void vector_mult_int(int *u, int *v, int count, int *out) {
    for (int i = 0; i < count; i++) {
        out[i] = u[i] * v[i];
    }
}

/**
 * @brief Return the dot product of two vectors of the same type and length
 * 
 * @param u 
 * @param v 
 * @param count 
 * @return double 
 */
double dot_product_double(double *u, double *v, int count) {
    double sum = 0;
    for (int i = 0; i < count; i++) {
        sum += u[i] * v[i];
    }
    return sum;
}

int dot_product_int(int *u, int *v, int count) {
    double sum = 0;
    for (int i = 0; i < count; i++) {
        sum += u[i] * v[i];
    }
    return sum;
}

void linear_interpolation(
    double *points,  // 1d array of points laid out in chunks of n_dims
    int *nearest_points,
    double *offset_vectors,
    double *amp,
    double *temp_amp,
    int *dim_lengths,
    int n_dims,
    int n_points
) {
    int i, j, k, dir, n_dirs, lin_idx, dir_idx, 
        to_linear_index_factors[n_dims],
        prev_directions[n_dims],
        temp_int_arr1[n_dims],
        temp_int_arr2[n_dims];
    double delta, maximum_vector[n_dims], temp_double_arr[n_dims];

    int total_bins = 1; // can move this up decleration

    // maximum_vector is furthest corner index of bin, but in negative direction
    // total_bins used to instantiate temporary amp array
    for (i = 0; i < n_dims; i++) {
        maximum_vector[i] = (double) (1 - dim_lengths[i]);
        total_bins *= dim_lengths[i];
    }

    // Taking dot product of to_linear_index_factors and some cartesian index
    // will give the correct index along a flattened array
    to_linear_index_factors[n_dims-1] = 1;
    for (i = n_dims - 2; i >= 0; i--) {
        to_linear_index_factors[i] = to_linear_index_factors[i+1] * dim_lengths[i];
    }

    for (i = 0; i < n_points * n_dims; i += n_dims) {

        printf("Point: (%f, %f) \n", points[i], points[i+1]);

        // Ensure point is within bounds of grid, otherwise skip binning point
        vector_add_double(points+i, maximum_vector, n_dims, temp_double_arr);
        if (vector_min_double(points+i, n_dims) < 0 || vector_min_double(temp_double_arr, n_dims) > 0) {
            printf("min %f\n", vector_min_double(points+i, n_dims));
            printf("max %f\n", vector_min_double(temp_double_arr, n_dims));
            printf("Out of bounds\n");
            continue;
        }

        memset(prev_directions, 0, sizeof(prev_directions));
        for (j = 0; j < n_dims; j++) {
            n_dirs = (int) pow((double) 2, (double) j);  // 2^j different combos
            delta = offset_vectors[i+j];

            printf("%f\n", delta);
            
            // Check if delta is zero or right or left 
            if (delta == 0) continue;
            if (delta < 0) {
                dir = -1;
                delta = -delta;
            }
            else dir = 1;

            printf("next\n");
            
            // Generate all combinations of previous directions without actually storing them
            memset(temp_int_arr1, 0, sizeof(temp_int_arr1));
            memset(temp_amp, 0, sizeof(temp_amp));
            temp_amp[dot_product_int(to_linear_index_factors, nearest_points+i, n_dims)] = 1;
            for (k = 0; k < n_dirs; k++) {

                printf("k: %d\n", k);

                // temp_int_arr1 now holds binary array of 0 and 1
                int_to_bin_digit(k, n_dims, temp_int_arr1);

                // temp_int_arr2 now holds one combination of direction
                vector_mult_int(temp_int_arr1, prev_directions, n_dims, temp_int_arr2);
                
                // temp_int_arr1 now holds cartesian index to nearest bin with offset in one direction
                vector_add_int(temp_int_arr2, nearest_points+i, n_dims, temp_int_arr1);

                // index to center bin and bin pointed to by direction combo
                lin_idx = dot_product_int(to_linear_index_factors, temp_int_arr1, n_dims);
                dir_idx = lin_idx + (to_linear_index_factors[j] * dir);

                temp_amp[dir_idx] = temp_amp[lin_idx] * delta;
                temp_amp[lin_idx] *= 1 - delta;
            }
            prev_directions[j] = dir;
            vector_add_double(amp, temp_amp, total_bins, amp);
        }
    }
}

// int main() {
//     printf("Hello World!\n");

//     // // double x[4] = {1.0, 2.0, 3.0, 4.0};
//     // // double y[4] = {2.0, 3.0, 4.0, 5.0};
//     // int x[8] = {1, 2, 3, 4, 5, 6, 7, 8};
//     // int y[4] = {2, 3, 4, 5};
//     // int res[4];

//     // int n = sizeof(x) / sizeof(x[0]);
//     // printf("%d\n", sizeof x);

//     // for (int i = 0; i < 8; i++) {
//     //     printf("%d ", x[i]);
//     // }

//     // memset(x, 0, sizeof(x));

//     // for (int i = 0; i < 8; i++) {
//     //     printf("%d ", x[i]);
//     // }


//     // elementwise_vector_multiplication(res, x, y, n);
//     // for (int loop = 0; loop < 4; loop++) {
//     //     printf("%d ", res[loop]);
//     // }
//     // printf("\n");

//     // int prod = dot_product(x, y, n);
//     // printf("%d\n", prod);

//     // int dims[5] = {10, 10, 10, 10, 10};
//     // linear_interpolation(NULL, dims, 5, NULL, NULL, NULL);

//     int x[4] = {1, 2, 3, 4};
//     int y[4] = {2, 4, 6, 8};
//     vector_add_int(x, y, 4, x);
//     for (int loop = 0; loop < 4; loop++) {
//         printf("%d ", x[loop]);
//     }
//     printf("\n");

//     return 0;
// }