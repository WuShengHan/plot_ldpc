
This repo open some sources
這裡開放一些parity matrix提供驗證比較

 spec: ALL QC64 matrix fits in one-soft-in-with-hw-3bit Decoder 
       AVX-C bit-true to 50Gbps-380K-gates RTL-IP
       
 mode 0: H8832_8192_QC64_DEG4
 mode 1: H5120_4096_QC64_DEG4
 mode 2: H2944_2048_QC64_DEG4
 mode 3: H2048_1024_QC64_DEG4
 mode 4: H1152_ 512_QC64_DEG4
 mode 5: H 768_ 256_QC64_DEG4
 mode 6: H1280_1024_QC64_DEG4

compile command
g++ -o run_wsh -mavx2 -O3 -march=native main_wsh.cpp libwsh.a -lpthread

execute
time ./run_wsh

1. encoding throughput 

this result is so high-speed due to fully utilize width=64 and deg-4.

2. decoding throughput

this allowes to detect error floor and verify coding family over 1000 gigabytes(1 terabyte) in 10 hours.
Welcome to compare with 5g-ldpc and 5g-polar coding.

3. more, decoding over hard-input.
   This also allowes you to change ch_idx into binary channel model through bpsk_awgn function.
   The outstanding parity matrix of WSH coding family are selected by this criterions.
   Althrough, the decoder needs another LLR table for hard calculation.





