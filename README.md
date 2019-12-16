# GPU-SFFT (Repository under construction, i.e., no ready yet, thanks)
The Sparse Fast Fourier Transform (MIT-SFFT) is an algorithm to compute the discrete Fourier transform of a signal with a sublinear time complexity,  i.e. algorithms with runtime complexity proportional to the sparsity level k,  where k is the number of non-zero coefficients of the signal in the frequency domain. The GPU-SFFT software is a  highly scalable GPU-based parallel algorithm for computing the SFFT of k-sparse signals. Our implementation of GPU-SFFT is based on parallel optimizations that leads to enormous speedups. These include carefully crafting parallel regions in the sequential MIT-SFFT code to exploit parallelism, and minimizing data movement between the CPU and the GPU. This allows us to exploit extreme parallelism for the CPU-GPU architectures and to maximize the number of concurrent threads executing instructions. Our experiments  show that our designed CPU-GPU specific optimizations lead to enormous decrease in the run times needed for computing the SFFT. Further we show that GPU-SFFT is 38x times faster than the MIT-SFFT and 5x faster than cuFFT, the NVIDIA CUDA Fast Fourier Transform (FFT) library. More details about GPU-SFFT in our paper cited in the Publications section that follows.
# Prerequisites
This software has been tested on the following dependences:
* CUDA 10.1
* gcc 5.4.0 
* Ubuntu 16.04.5

# Install
Series of instructions to compile and run your code

1. Download the software

git clone https://github.com/pcdslab/GPU-SFFT.git

2. Compile the code with the Makefile. Please set the library and include directories paths on the Makefile available for your particular machine

$ cd GPU-SFFT

$ make
# Run 
$ ./exec -q -k -Bcst -loops_loc -loops_thresh -loops_est -RLFC_cst -RLFC_loops 

|argument|description|
|--------|-----------|
|q |n = 2^q, n: size of the input signal|
|k |sparsity of the input signal, i.e., the number of frequency coefficients different than zero|
|Bsct|Parameter to compute the number of bins used to estimate the locations and the values of largest frequencies coefficients.|
|loops_loc|Number of loops to estimate the locations of the largest frequencies coefficients.|
|loops_thresh|Number of times a frequency coefficient has to appear in the location bins.|
|loops_est|Number of loops to estimate the values of the largest frequencies coefficients.|
|RLFC_cst|Constant to compute the number of bins used to compute the restriction of largest frequency coefficients.|
|RLFC_loops|Number of loops to run the function to compute the restriction of largest frequency coefficients.|
# Experiments

All the experiments to test the GPU-SFFT software were performed on a Linux server with Ubuntu operating system version 16.04.5, 44 Intel Xeon Gold processors, clock speed 2.1 GHz, and 125 GB of RAM. The GPU in this server is a NVIDIA Titan Xp, with 30 SM, 128 cores/SM, maximum clock rate of 1.58 GHz, 12196 MB of  global memory, L2 cache size of                             3145728 bytes, 65536 bytes of constant memory, 49152 bytes of shared memory per block, 65536 registers available per block, 
a maximum number of 2048 threads per multiprocessor, a maximum number of 1024 threads per block, maximum dimension size of a thread block (x,y,z): (1024, 1024, 64), maximum dimension size of a grid size (x,y,z): (2147483647, 65535, 65535),
 and CUDA version 10.1 with CUDA capability of 6.1.


# Parameters for experiments
k= 1000

|q |19|20|21|22|23|24|25|26|
|--|--|--|--|--|--|--|--|--|
|Bcst |3.85|1.40|2.16|0.683|0.99|0.663|0.662|0.478|

q = 27, n = 2^27

|k |1000|4000|7000|10000|13000|16000|19000|22000|
|--|--|--|--|--|--|--|--|--|
|Bcst|0.68|0.468|0.77|0.7|0.665|0.65|0.66|0.8|
|RLFC_cst|128|1024|4096|4096|4096|4096|4096|4096|

|k |25000|28000|31000|34000|37000|40000|43000|
|--|--|--|--|--|--|--|--|
|Bcst|0.79|0.71|0.69|0.71|0.7|0.7|0.8|
|RLFC_cst|8192|16384|16384|16384|32768|32768|32768|



# Ouput sample
oarti001@raptor:~/documents/FIU_cudaSFFT/FIU_cudaSFFT_SaeedLab$ ./exec 27 16000 0.65 2 2 7 4096 2

RUNNING EXPERIMENT: n:134217728 k:16000 loops_loc=2 loops_est:7

GPU-SFFT filter parameters for: n=134217728, k=16000.

******************************************************************************

 RLFC Filter: RLFC_loops: 2  B_thresh:32000 W_RLFC: 4194304
 
 Filter: (numlobes=183314.0, tolerance=1e-08, b_f=966) B: 131072 

 computing filter starts
 
 Window size: Filter : 2230607
 
 Time for computing  filter:  3.709905s
 
******************************************************************************

Simulation starting:

GPU-SFFT Results

******************************************************************************

Total GPU-SFFT time: 0.172697

Time distribution: 
                   
         votetable  restr-loc-index  perm+filter reverse-hash  estimation FFT+Cutoff    other    total

           0.0007       0.0473          0.0227      0.0112        0.0714     0.0163      0.0031   0.1727
                    
             0.4%         27.4%           13.2%       6.5%         41.3%       9.5%        1.8%   100.0%
                     

HtoD time = 0.232404s, HtoD time/total time  = 134.573197%, total time + HtoD time = 0.405101s

ERROR:

K=16000; MISSED (estimation, result) = (0, 0), MSE*k = 1.19875, Mean Square Error (MSE) = 7.49222e-05

time to extract the values of the largest frequencies: 0.042122s

******************************************************************************

cuFFT(NVIDIA-CUDA) library results

******************************************************************************

Memory transfer host to device execution time 233.19milisec

Plan creation execution time 2.05ms

cuFFT execution time    91.65ms

Memory transfer device to host execution time 1066.93ms

total time 1391.77ms

******************************************************************************

# Publications

If you use this software please cite our paper:

Oswaldo Artiles, and Fahad Saeed, “GPU-SFFT: A GPU based parallel algorithm for computing the Sparse Fast Fourier Transform (SFFT) of k-sparse signals”, Workshop on Performance Engineering with Advances in Software and Hardware for Big Data Sciences (PEASH), Proceedings of IEEE Conference on Big Data (IEEE BigData 2019), Los Angeles, CA, USA, Dec. 09-12, 2019 

# Acknowledgements
This research was supported by the National Science Foundations (NSF) under the Award Numbers CAREER OAC-1925960. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Science Foundation. We would also like to acknowledge the donation of a K-40c Tesla GPU and a TITAN Xp GPU from NVIDIA which was used for all the GPU-based experiments performed with the software.




