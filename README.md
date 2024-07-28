# Parallel_Aligner
A naive progressive alignment algorithm for protein sequences which applies vectorisation (AVX), cluster computing (MPI) and threading (OpenMP)

For an explanation of the algorithm and benchmarking results, watch the following [video](https://youtu.be/DkpAvOKyZMg).

# Compilation 

If you wish to use the parallel implementation you will need Math Kernel Library (MKL) and OpenMPI installed on your machine. Refer to the Makefile for details. 

There is also a serial implementation in the `serial` directory. 

```
cd <path-to-download> Parallel_Aligner/
```

```
make
```

# Usage

```
./msaParallel <path-to-data>/example.fasta
```

In the current implementation, you must execute the program with mpiexec. An example command would look like this: 
```
mpiexec -n 2 -map-by node -bind-to none ./msaParallel ./data/globin/100_seqs_globin.fasta 
```
