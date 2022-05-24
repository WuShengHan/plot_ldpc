#
Here provides comparison material to implement the optimal FEC in B5G/6G system

X-axis is noticed that Es/N0 is direct linked to raw bit error rate (channel condition)
Eb/N0 is for fair comparison among the same info length with different parity size

Y-axis also must be noticed more if next-gen communication. 10^3 CWs(or frames) means around 100KB - 1MB transmission testing.
10^6 CWs(or frames) means around 100MB - 1GB transmission testing.
10^8 or more CWs(or frames) will reach over-100GB-data reliable and high-quality communication to match every use case.

The parameter of code rate is suitable for turbo code, convolution code but not for ldpc code.
LDPC design focuses on parity length (number of check node), qc-size(number of parallel processing), codeword length (time cost) 

5G-RANGE: Remote Area Access Network for the 5th Generation
Deliverable 3.1 Physical layer of the 5G-RANGE
https://ec.europa.eu/research/participants/documents/downloadPublic?documentIds=080166e5c0722c76&appId=PPGMS

Comparison of Polar Decoders with Existing Low-Density Parity-Check and Turbo Decoders
https://arxiv.org/pdf/1702.04707.pdf

A Comparison of Channel Coding Schemes for 5G Short Message Transmission
https://www.semanticscholar.org/paper/A-Comparison-of-Channel-Coding-Schemes-for-5G-Short-Iscan-Lentner/87c7417fe10b2efb4b936270cad1e8216ac57d2b

BER Comparison Between Convolutional, Turbo, LDPC, and Polar Codes
https://publik.tuwien.ac.at/files/publik_262129.pdf

On the Evaluation of the Polyanskiy-Poor-Verdu Converse Bound for Finite Block-length Coding in AWGN
https://arxiv.org/pdf/1401.7169.pdf

Physical Layer Latency Management Mechanisms: A Study for Millimeter-Wave Wi-Fi
https://www.mdpi.com/2079-9292/10/13/1599/htm

Link Budget Analysis: Error Control & Detection by Atlanta RF
https://atlantarf.com/Error_Control.php

#
