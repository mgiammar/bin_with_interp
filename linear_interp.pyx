from numpy cimport ndarray
import numpy as np
cimport numpy as cnp

cdef extern from "bin_interp.h":
    void linear_interpolation(
        double *points,
        int *nearest_points,
        double *offset_vectors,
        double *amp,
        double *temp_amp,
        int *dim_lengths,
        int n_dims,
        int n_points
    )
    int vector_min_int(int *vec, int count)
    int dot_product_int(int *u, int *v, int count)
    void vector_add_int(int *u, int *v, int count, int *out)
    void vector_mult_int(int *u, int *v, int count, int *out)
    void vector_add_double(double *u, double *v, int count, double *out)
    void vector_mult_double(double *u, double *v, int count, double *out)
    void int_to_bin_digit(int num, int count, int* out)
    double vector_min_double(double *vec, int count)
    double dot_product_double(double *u, double *v, int count)

def py_vector_min_int(list a):
    cdef int count = len(a)
    cdef ndarray[int, ndim=1] vec = np.array(a, dtype=np.int32)

    return vector_min_int(&vec[0], count)

def py_dot_product_int(list a, list b):
    cdef int count = len(a)
    cdef ndarray[int, ndim=1] u = np.array(a, dtype=np.int32)
    cdef ndarray[int, ndim=1] v = np.array(b, dtype=np.int32)

    return dot_product_int(&u[0], &v[0], count)

def py_vector_add_int(list a, list b):
    cdef int count = len(a)
    cdef ndarray[int, ndim=1] u = np.array(a, dtype=np.int32)
    cdef ndarray[int, ndim=1] v = np.array(b, dtype=np.int32)
    cdef ndarray[int, ndim=1] out = np.zeros((count), dtype=np.int32)

    vector_add_int(&u[0], &v[0], count, &out[0])
    return out

def py_vector_mult_int(list a, list b):
    cdef int count = len(a)
    cdef ndarray[int, ndim=1] u = np.array(a, dtype=np.int32)
    cdef ndarray[int, ndim=1] v = np.array(b, dtype=np.int32)
    cdef ndarray[int, ndim=1] out = np.zeros((count), dtype=np.int32)

    vector_mult_int(&u[0], &v[0], count, &out[0])
    return out

def py_vector_min_double(list a):
    cdef int count = len(a)
    cdef ndarray[cnp.float64_t, ndim=1] vec = np.array(a, dtype=np.float64)

    return vector_min_double(&vec[0], count)

def py_vector_add_double(list a, list b):
    cdef int count = len(a)
    cdef ndarray[cnp.float64_t, ndim=1] u = np.array(a, dtype=np.float64)
    cdef ndarray[cnp.float64_t, ndim=1] v = np.array(b, dtype=np.float64)
    cdef ndarray[cnp.float64_t, ndim=1] out = np.zeros((count), dtype=np.float64)

    vector_add_double(&u[0], &v[0], count, &out[0])
    return out

def py_vector_mult_double(list a, list b):
    cdef int count = len(a)
    cdef ndarray[cnp.float64_t, ndim=1] u = np.array(a, dtype=np.float64)
    cdef ndarray[cnp.float64_t, ndim=1] v = np.array(b, dtype=np.float64)
    cdef ndarray[cnp.float64_t, ndim=1] out = np.zeros((count), dtype=np.float64)

    vector_mult_double(&u[0], &v[0], count, &out[0])
    return out

def py_dot_product_double(list a, list b):
    cdef int count = len(a)
    cdef ndarray[cnp.float64_t, ndim=1] u = np.array(a, dtype=np.float64)
    cdef ndarray[cnp.float64_t, ndim=1] v = np.array(b, dtype=np.float64)

    return dot_product_double(&u[0], &v[0], count)

def py_int_to_bin_digit(int _num):
    cdef int num = _num
    cdef int count = 32
    cdef ndarray[int, ndim=1] out = np.zeros((32), dtype=np.int32)

    int_to_bin_digit(num, count, &out[0])
    return out

def py_dot_product_int(list a, list b):
    cdef int count = len(a)
    cdef ndarray[int, ndim=1] u = np.array(a, dtype=np.int32)
    cdef ndarray[int, ndim=1] v = np.array(b, dtype=np.int32)
    
    return dot_product_int(&u[0], &v[0], count)


def py_numpy_pointer_tests(
    ndarray[int, ndim=1] a,
    ndarray[int, ndim=2] b,
    ndarray[int, ndim=3] c,
):
    b = b.T
    print(a)
    print(b)
    print(c)

    cdef int* a_ptr = &a[0]
    cdef int* b_ptr = &b[0, 0]
    cdef int* c_ptr = &c[0, 0, 0]
    print(<size_t>a_ptr)
    print(<size_t>b_ptr)
    print(<size_t>c_ptr)

    print("a")
    print(a_ptr[0])
    print(a_ptr[0])
    print(a_ptr[0])
    print(a_ptr[0])

    print("b")
    print(b_ptr[0])
    print(b_ptr[1])
    print(b_ptr[2])
    print(b_ptr[3])
    print(b_ptr[4])
    print(b_ptr[5])
    print(b_ptr[6])
    print(b_ptr[7])
    print(b_ptr[8])
    print(b_ptr[9])

    print("c")
    print(c_ptr[0])
    print(c_ptr[1])
    print(c_ptr[2])
    print(c_ptr[3])
    print(c_ptr[4])


def py_bin_with_interpolation(
    ndarray[cnp.float64_t, ndim=2] _points,
    ndarray[cnp.float64_t, ndim=2] _grid_dims
):
    """
    Arguments:
        (np.ndarray) points: random points
        (np.ndarray) grid_dims: Numpy array of each dimension on the grid
    """
    # Transform grid into integers and random points onto that grid
    for i in range(len(_grid_dims)):
        dim = _grid_dims[i]
        min = dim.min()
        width = (dim.max() - min) / dim.size
        _points[i] = (_points[i] - min) / width

    _points = np.vstack(_points).T
    nearest = np.rint(_points)
    distance = _points - nearest
    cdef ndarray[cnp.float64_t, ndim=2] points = _points.copy()
    cdef ndarray[cnp.float64_t, ndim=2] nearest_distance = distance
    cdef ndarray[int, ndim=2] nearest_points = nearest.astype(np.int32)

    cdef ndarray[int, ndim=1] dim_sizes = np.asarray(
        [dim.size for dim in _grid_dims], dtype=np.int32
    )
    cdef ndarray[cnp.float64_t, ndim=1] amp = np.zeros(dim_sizes, dtype=np.float64).flatten()
    cdef ndarray[cnp.float64_t, ndim=1] temp_amp = np.zeros(dim_sizes, dtype=np.float64).flatten()

    cdef int n_dims = dim_sizes.size
    cdef int n_points = points.size / n_dims

    print(points)

    linear_interpolation(
        &points[0, 0],
        &nearest_points[0, 0],
        &nearest_distance[0, 0],
        &amp[0],
        &temp_amp[0],
        &dim_sizes[0],
        n_dims,
        n_points,
    )

    return amp.reshape(dim_sizes)