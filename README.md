<pre> 
#
This python code would demostrate my LDPC code performance,
which include
    1. WSH code family for codeword info = 128, 256, 512, 1024, 2048, 4096 and 8320
    2. 5G NR NTU code proposed family 
       http://www.3gpp.org/ftp/TSG_RAN/WG1_RL1/TSGR1_AH/NR_AH_1701/Docs/R1-1700645.zip
       http://www.3gpp.org/ftp/tsg_ran/WG1_RL1/TSGR1_88/Docs/R1-1703208.zip
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
<pre>
