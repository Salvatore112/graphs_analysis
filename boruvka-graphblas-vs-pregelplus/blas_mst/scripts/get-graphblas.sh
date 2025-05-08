export GRAPHBLAS_INCLUDE_DIR=`pwd`/graphblas-binaries/include
export GRAPHBLAS_LIBRARY=`pwd`/graphblas-binaries/lib/libgraphblas.so
mkdir graphblas-binaries
cd graphblas-binaries
wget --quiet https://anaconda.org/conda-forge/graphblas/7.4.1/download/linux-64/graphblas-7.4.1-hcb278e6_0.conda
unzip graphblas-7.4.1-hcb278e6_0.conda
tar xf pkg-graphblas-7.4.1-hcb278e6_0.tar.zst
cd ..

