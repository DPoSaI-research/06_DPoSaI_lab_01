# Fourier Transform

This is the source code for lab 01 on DPoSaI.
Transforms are implemented as function templates under different namespaces to
distinguish between trivial and fast implementations.

See [Fourier.h](Fourier.h) for details.

Fast transform uses trivially implemented recursive functions under anonymous
namespace. This functions use a lot of valarray instances with duplicated data
and function calls; move data several times across result vector during
computation.
