
MPI_INC = /usr/include/openmpi-x86_64
MPI_LIB = /usr/lib64/openmpi/lib

CC = g++
CFLAGS = -std=c++11 -Wall -O3 -pg -Wpedantic -ffast-math -mavx2 -fopenmp -I$(MPI_INC) 
SRC = msa.cpp matrix.cpp
HEADERS = msa.h matrix.h
TARGET = msaAvx


$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) -L$(MPI_LIB) -lmpi_cxx -lmpi

clean:
	rm -f msa
	
.PHONY: clean


