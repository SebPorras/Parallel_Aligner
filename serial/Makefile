CC = g++
CFLAGS = -std=c++11 -Wall -O3 -pg -Wpedantic -ffast-math
SRC = msa.cpp matrix.cpp
HEADERS = msa.h matrix.h
TARGET = msa

$(TARGET): $(SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f msa
	
.PHONY: clean


