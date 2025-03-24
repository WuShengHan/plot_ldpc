# WSH Coding Family - Parity Matrix Verification

This repository provides open-source parity matrices for verification and comparison.

## Specification

- **All QC64 matrices** fit within a **one-soft-input-with-hw-3bit decoder**.
- **Implementation** AVX-C code bit-true with optimized RTL IP  .
- **Cost-Performance** 380K gate count achieves over **40Gbps** throughput.

## Supported Modes

| Mode | Matrix Configuration |
|------|----------------------|
| 0    | H8832_8192_QC64_DEG4 |
| 1    | H5120_4096_QC64_DEG4 |
| 2    | H2944_2048_QC64_DEG4 |
| 3    | H2048_1024_QC64_DEG4 |
| 4    | H1152_512_QC64_DEG4  |
| 5    | H768_256_QC64_DEG4   |
| 6    | H1280_1024_QC64_DEG4 |

## Updates (2025-03-24)

### a.  Verification Adjustment  
For those with experience in ECC/FEC verification, it is known that in `main_wsh.cpp`, modifying `line 228` and `line 237` from `±0.5` to `±0.72` can significantly reduce the FER.  

- For `SNR = 3 dB`, the FER improves from `7.69231e-05` to `2.7027e-05` when using the code `N2048_K1024`.  

However, this adjustment only modifies the channel behavior, altering the soft value distribution without actually improving the decoder's capability or error resilience.  

### b.  Biased Codewords  
As far as I know, another approach is to introduce noise into codewords that contain more `bit-0`s to decode, summarize and obtain a pretending-better scenario.  

- A particularly effective method is to use the all-zero codeword and leverage calculations where `LLR = 0` to determine `bit-0` values, allowing for more accurate guesses that lead to a converged codeword.  

This technique is beneficial when handling biased data with a high occurrence of `bit-0`s, yielding improved results. However, in real-world applications, data distributions may vary, leading to differences in practical performance and user experience.  


## Compilation & Execution

Compile using:  

```sh
g++ -o run_wsh -mavx2 -O3 -march=native main_wsh.cpp libwsh.a -lpthread
```

Executing using:  

```sh
./run_wsh 
```

## Encoding Throughput
Testing encoding for 10,000,000 codewords: 

```
real    0m1.730s
user    0m1.706s
sys     0m0.012s
```

Extremely high-speed performance achieved by fully utilizing width = 64 and degree = 4.  
  
## Decoding Throughput  
Perform your own tests.
This is under one-soft-bit decoding if you review the `bpsk_awgn` function.  
This implementation could bit-true match verilog-RTL ip.  
The `main_wsh.cpp` C code allows error floor detection and verification over 1 terabyte of data in 10 hours.  
Welcome to compare against 5G LDPC and 5G Polar coding.  

## Features
Decoding over hard-inputs: Supports binary channel modeling using the `bpsk_awgn` function.   
I use another LLR table to run hard and refine matrix sieving.  

- **WSH Coding Family**: These are partial released. If interested, please feel free to contact me via linkedin.






