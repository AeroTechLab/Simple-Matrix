# Simple Matrix

A set of basic C routines to abstract vector/matrix storage and operations, offering:

- Matrix memory management (creation, deletion, copy, resizing, etc.)
- Reading/writing matrix values for single elements or as a whole through raw buffers ([row-major order](https://en.wikipedia.org/wiki/Row-_and_column-major_order))
- Matrices/vectors sum and multiplication
- Transpose of a matrix
- Inverse and determinant of a square matrix
- Matrix formatted printing

Internally, the library uses [BLAS/LAPACK](https://en.wikipedia.org/wiki/LAPACK) routines, so the library must be linked to one of its available implementations, like the [reference BLAS/LAPACK](http://www.netlib.org/lapack/lug/node11.html), [OpenBLAS](http://www.openblas.net/), [ATLAS](http://math-atlas.sourceforge.net/), [Intel's MKL](https://software.intel.com/en-us/intel-mkl), etc.

### Building example

For instance, building this library with [GCC](https://gcc.gnu.org/) as a shared object, using reference **BLAS/LAPACK**, would require the shell command (from root directory):

>$ gcc matrix.c -I. -shared -fPIC -o matrix.so -lblas -llapack

### Documentation

Descriptions of how the functions and data structures work are available at the [Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html)-generated [documentation pages](https://labdin.github.io/Simple-C-Matrix-Library/)
