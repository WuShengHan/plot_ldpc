# 
#
This python code would demostrate my LDPC code performance,
which include
    1. WSH code family for codeword length from 128 to 8320
    2. 5G NR NTU code proposed family (On-going)
    3. IP 1KB code
    4. IP 2KB code
    5. IP 4KB code

It is noticeable that
    The result of 1,2   is shown in soft-input correction
    The result of 3,4,5 is shown in hard-input correction

Here provide some figure property
    a. raw bit error rate vs. frame error rate
    b. signal noise ratio vs frame error rate, containing energy-bit signal and energy-symbol signal
    c. exact error bit vs eaxct frame error rate, which means bit-err is directly compared with BCH correction
    d. exact error bit vs throughput: demo hard-input throughput
