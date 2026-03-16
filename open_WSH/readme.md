# WSH Coding Family - Parity Matrix Verification

This repository provides open-source parity matrices for verification and comparison.

## Specification

- **All QC64 matrices** fit within a **one-soft-input-with-hw-3bit decoder**.
- **Implementation** AVX-C code bit-true with optimized RTL IP  .
- **Cost-Performance** 380K gate count achieves over **40Gbps** throughput.

## Supported Modes

|Mode|Matrix Configuration|Description|
|-|-|-|
| 0    | H8832_8192_QC64_DEG4  | ^ |
| 1    | H8960_8192_QC64_DEG4  | ^ |
| 2    | H9088_8192_QC64_DEG4  | ^ longer info length |
| 3    | H9216_8192_QC64_DEG4  | ^ higher rate |
| 4    | H4736_4096_QC64_DEG4  | ^ to gain TxRx info efficiency |
| 5    | H4864_4096_QC64_DEG4  | ^ |
| 6    | H4992_4096_QC64_DEG4  | ^ |
| 7    | H5120_4096_QC64_DEG4  | ^ |
| 8    | H2688_2048_QC64_DEG4  |   |
| 9    | H2816_2048_QC64_DEG4  | V |
| 10   | H2944_2048_QC64_DEG4  | V |
| 11   | H3072_2048_QC64_DEG4  | V shorter info length |
| 12   | H1664_1024_QC64_DEG4  | V lower rate |
| 13   | H1792_1024_QC64_DEG4  | V to gain reliability |
| 14   | H1920_1024_QC64_DEG4  | V |
| 15   | H2048_1024_QC64_DEG4  |   |
| ...  | more modees are ok    | stay tuned  |


## Release Verilog RTL IP & Cosim with Verilator and AVX-C (2026-03)

Open to announce the latest verification env of our LDPC IP core, and featuring high-performance simulation capabilities.

### A. Browser-Based Simulation
You can still conveniently launch LDPC simulations directly in your browser. 
*   **Deep Analysis**: Run long-term simulations by local run or browser to generate a detailed `report.txt`.
*   **Design Coherence**: ensure the 100% match from Verilog IP to AVX-C by its results and timing cycles
*   **Data Visualization**: Use the exported report data to plot BER/FER curves for your personal comparison.

### B. Comprehensive IP Evaluation
Users are encouraged to independently verify the following metrics to ensure the IP meets your design requirements:
*   **Performance**: Throughput and Latency analysis.
*   **Power**: Estimated power consumption profile.
*   **Area Cost**: Resource utilization (Gate count / LUTs / Flip-Flops).
*   **Error Floor**: Correction performance and error floor verification.
*   **FPGA Behavior**: Validated RTL behavior on target FPGA platforms.

---

### License & Terms of Commercial Use

By using this Verilog RTL IP, you agree to the following **Dual-Licensing** terms:

#### Any use of this IP by a **company, corporation, or for-profit entity** requires a formal **Commercial License Agreement**. This includes:
*   Integrating this IP into a commercial SoC, ASIC, or FPGA product.
*   Using this IP for internal corporate R&D, prototyping, or evaluation.
*   **Unauthorized commercial use will be subject to legal action.**

#### Licensing Inquiry
To obtain the **Full RTL Source Code**, **Technical Support**, or a **Commercial License**, please contact the author: ShengHan Wu, from Taiwan, WuShengHan@outlook.com

---|

## Directly run code on Google Colab (2025-04-24)

### a. Launch the LDPC simulation in your browser
This Colab notebook allows you to reproduce the LDPC decoding results with high efficiency — reaching FER down to $10^{-4}$ within approximately 60 seconds.
https://colab.research.google.com/github/WuShengHan/plot_ldpc/blob/main/open_WSH/colab_run.ipynb

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






