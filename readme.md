# Fourier Transform

This is the source code for lab 01 on DPoSaI.
Transforms are implemented as function templates under different namespaces to
distinguish between trivial and fast implementations.

See [Fourier.h](Fourier.h) for details.

Fast transform is implemented in three different ways.

The first version uses trivially implemented recursive function under anonymous
namespace. This version uses a lot of valarray instances with duplicated data
and function calls; moves data several times across result vector during
computation.

The second version uses valarray slices and does not create redundant instances,
all computations are made in place. Returns generalized slice that refers to
properly ordered output vector elements.

The third version is not recursive. It uses prebuilt recursion tree to process
over data in vector. Supposed to be the fastest implementation.
