%High speed MATLAB codes
 
%This program use SCL with small list size as a "filter", i.e., under the
%same noise realization, If CRC-SCL with smaller L is correct, then CRC-SCL
%with larger L must be correct. You may not believe in this conjecture, but you can have
%a try, This is true with high probability (almost 1)
 
%This program has another accelerator. If in lower snr, CRC-SCL is correct, then at higher snr, 
%the same CRC-SCL must be correct.  You may not believe in this conjecture, but you can have
%a try, This is true with high probability (almost 1).
 
%I know above methods sound dangerous, but it is safe under following 4
%conditions
 
%1. The same code word
%2. The same AWGN noise realization with distribution N(0, 1) 
%3. CRC must be used (Only SCL is not permitted)
%4. The same code construction in all SNR range.
 
%This program satisfies above 4 conditions. You can verify the bler performance by comparisons with existing results. 
 
%Since MATLAB is not good at recursive function (too many parameters to be passed)
%I cancel the well-known recursiveCalcP() and recursiveCalcB() proposed by
%I. Tal
%Instead, 'For' function is used. You may argue that we can use objective
%oriented (OO) style. However, OO in matlab is also slow.
 
%Besides, the algorithms of following papers are provided.
 
%How to Construct Polar Codes
%Fast Successive-Cancellation Decoding of Polar Codes: Identification and Decoding of New Nodes
%beta-expansion: A Theoretical Framework for Fast and Recursive Construction of Polar Codes
 
%Above three algorithms are made by myself so the correctness is not guaranteed. 

%Modify by zwy in 07/27/2019
%Come from yu yong run

clear
clc

Ng = 12;
poly = [1 1 1 1 1 0 0 0 1 0 0 1 1];
n = 9;  % you should change the channel parameter simutaneously
N = 2^n;
K = 256+12;
Kp = 42;
ebno_vec = [1 2 3 3.5 4]; %row vec, you can write it like [1 1.5 2 2.5 3] 
list_vec = [1 2 4 8];  %row vec, you can write it like [1 4 16 32 ...]. The first element is always 1 for acceleration purpose. The ramaining elements are power of two.
max_runs = 1e7;
max_err = 100;
resolution = 1e4;%the results are shown per max_runs/resolution.

load('Pe_N512_snr3_r5.mat');
[bler, ber] = simulation(N, K, Kp, max_runs, max_err, P, resolution, ebno_vec, list_vec,Ng, poly);

