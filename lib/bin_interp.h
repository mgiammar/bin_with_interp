#ifndef BIN_INTERP_H_
#define BIN_INTERP_H_

int vector_min_int(int *vec, int count);
int dot_product_int(int *u, int *v, int count);
void vector_add_int(int *u, int *v, int count, int *out);
void vector_mult_int(int *u, int *v, int count, int *out);
void vector_add_double(double *u, double *v, int count, double *out);
void vector_mult_double(double *u, double *v, int count, double *out);
void int_to_bin_digit(int num, int count, int* out);
double vector_min_double(double *vec, int count);
double dot_product_double(double *u, double *v, int count);
void linear_interpolation(
    double *points,
    int *nearest_points,
    double *offset_vectors,
    double *amp,
    double *temp_amp,
    int *dim_lengths,
    int n_dims,
    int n_points
);

#endif