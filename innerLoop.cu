/* 
 *  Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *  Eric Price, Massachusetts Institute of Technology.
 */
/*
 *  Copyright (C) 2019 Oswaldo Artiles and Fahad Saeed 
 *  Florida International University, Florida, USA.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE   
 * Please refer to the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>

//includes CUDA project
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <cuComplex.h>

// includes Thrust project
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>

extern "C"{
           #include "utils.h"
           #include "innerLoop.h"
           #include "outerLoop.h"
           #include "cudaFft.h"
           #include "timer.h"
}

#define vprintf(...) if (VERBOSE){ printf(__VA_ARGS__);}

////////////////////////prototype functions/////////////////////////////////////
void cutoff(int *output,int num,thrust::device_ptr<double> &d_samples_ptr,
	    double *h_samples,int n);
double nth_element_immutable_cuda(thrust::device_ptr<double> &d_samples_ptr,int n,
                                  int num);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
__global__ void printComDxKernel(complex_t *d_x,int n);
__global__ void printIntDxKernel(int *d_x,int n);
__global__ void printDouDxKernel(double *d_x,int n);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
__global__ void restLocFreqCoeffKernel(complex_t *d_x_filt_t, complex_t *d_x,
				       int W_RLFC,int offset,int sigma);
__global__ void samplesKernel(complex_t *d_samp,double *d_samples,
			      double *d_samples_u,int samp_size);
__global__ void permFilterTillKernel(complex_t *d_bins_t,complex_t *d_x,
				     complex_t *d_filter_t,int n,int B,int ai,
				     int tiles,int residue);
__global__ void permFilterKernel(complex_t *d_bins_t, complex_t *d_x,
				 complex_t *d_filter_t,int n,int B,int ai,
				 int d_filter_size);
__global__ void reverseHashKernel(int *d_I,int *d_vote,int *d_J,Pair *d_J2_sigma,
				  int n,int B_thresh,int B, int a,int loop_threshold,
				  int *d_I_F,int num_RLFC, int W_RLFC);
__global__ void estimateValuesKernel(Node *d_ans,int *d_I,int my_I_F,
				     complex_t *d_bins_f,int loops,int n,
				     int *d_permute,int B,complex_t *d_filter_f);
__global__ void cutoffKernel(int *d_output,double *d_cutoff,double *d_samples,
			     int samp_size,int *d_count,int B);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
__device__ int timesMod(const int x, const int a, const int n);
__device__ int findUpperBound(Pair *d_J, int n, int val);
__device__ void insertionSortDev(double *d_x, int n);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * This function implements a heuristic used in SFFT 2.0
 * as  a preprocessing stage to restrict the location
 * of the largest frequency coefficients in the DFT of
 * the input signal. The heuristic is based on a aliasing 
 * filter that has not leakage at all.
 * This stage  is implemented with the following steps:
 *
 * 1. Allocate device  memory  for the device filtered vector
 *    in the time domain, d_x_filt_t, with a size W_RLFC.
 * 2. Filter the device input vector d_x  and stores 
 *    the filtered output  in d_x_filt_t,
 * 6. Allocate device  memory for  the filtered vector in the
 *    frequency domain,d_x_filt_f,  with a size W_RLFC.    
 * 7. Create a plan to compute the DFT of d_x_filt_t and stores
 *    it in d_x_filt_f.
 * 8. Execute the plan.
 * 9. Initiate the  cutoff stage. 
 *
 */
int restr_loc_freq_coeff(complex_t *d_x, int n, int B_thresh, int W_RLFC, int* d_J2, double *RLFC_Time){

  /*timing variables */
  float t_filt = 0.0;
  float t_CUFFT_ex = 0.0;
  float t_total = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  int sigma = n/W_RLFC;
  int offset = (unsigned) floor( drand48() * sigma);

  /*Allocate device  memory  for the time domain sampling vector: d_x_filt_t*/
  cudaEventRecord(start);
  complex_t *d_x_filt_t;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_x_filt_t), sizeof(*d_x_filt_t)*W_RLFC));

  /*filtering input signal */
  dim3 dimBlock(256);
  dim3 dimGrid(W_RLFC/dimBlock.x);
  restLocFreqCoeffKernel<<<dimGrid, dimBlock>>>(d_x_filt_t, d_x, W_RLFC,offset, sigma);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_filt,start, stop);

  /*CUFFT plan creation*/
  cudaEventRecord(start);
  cufftHandle plan;
  int batch = 1;
  cudaEventRecord(start);
  checkCudaErrors(cufftPlan1d(&plan,W_RLFC,CUFFT_Z2Z,batch));

  /*Allocate device  memory for the frequency domain filtered vector:d_x_filt_f */
  complex_t *d_x_filt_f;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_x_filt_f), sizeof(*d_x_filt_f)*W_RLFC));

  /*CUFFT plan  execution*/
  checkCudaErrors(cufftExecZ2Z(plan,reinterpret_cast<cufftDoubleComplex *>(d_x_filt_t),
			       reinterpret_cast<cufftDoubleComplex *>(d_x_filt_f),CUFFT_FORWARD));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_CUFFT_ex,start, stop);
  checkCudaErrors(cudaFree(d_x_filt_t));

  t_total =  t_filt +t_CUFFT_ex;
  *RLFC_Time  += t_total/1.0e3;
  samples_indices(d_x_filt_f,B_thresh,W_RLFC,d_J2,RLFC_Time);

  int print_t = 0;
  if (print_t){
    printf("\ninnerLoop:restr_loc_freq_coeff:time filtering = %lfs \n", t_filt/1.0e3);
    printf("innerLoop:restr_loc_freq_coeff:time CUDAFFT execution = %lfs \n", t_CUFFT_ex/1.0e3);
    printf("innerLoop:restr_loc_freq_coeff:time restr_loc_freq_coeff =  %lfs, RLFC_Time = %lfs \n\n",t_total/1.0e3, *RLFC_Time);
  }

  /*cleanup memory*/
  checkCudaErrors(cufftDestroy(plan));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));
  checkCudaErrors(cudaFree(d_x_filt_f));

  return 0;
}// end restr_loc_freq_coeff
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * This function implements the steps: permutation and filtering,and  FFT 
 * and cutoff, of the inner location loop   of the GPU-SFFT algorithm,
 * as follows:
 * 1. Allocates device memory for the   permuted and filtered  vector, 
 *    d_bins_t, with the size = B. 
 * 4. Random spectrum permutation and apply filter in one loop:
 *    Randomly permutes x, computes the dot product of d_x and
 *    d_filter_t, and stores the result  of the dot product in 
 *    d_bins_t.
 * 6. Compute the DFT of d_bins_t and stores it in d_bins_f.
 * 7. Initiate the  cutoff stage. 
 *
 */
int perm_filter_cutoff(complex_t *d_x,int n,complex_t *d_filter_time,
		       int d_filter_size, int B_thresh,int B,int a,
		       int ai, complex_t *d_bins_f,int *d_J,
		       double *PF_T,double *BC_T){

  /*timing variables */
  float t_set_binst = 0.0;
  float  t_perm_filter = 0.0;
  float t_CUFFT_ex = 0.0;
  float t_total = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  if (n % B){
    fprintf(stderr, "Warning: n is not divisible by B, which algorithm expects.\n");
  }

  /*Allocate device  memory  for the time domain bins vector: d_bins_t*/
  cudaEventRecord(start);
  complex_t *d_bins_t;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_bins_t), sizeof(*d_bins_t)*B));

  /*Set to zero the components of  d_bins_t*/
  checkCudaErrors(cudaMemset(d_bins_t, 0, sizeof(*d_bins_t)*B));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_set_binst,start, stop);

  /*permute and filter the input signal d_x and store it in d_bins_t*/
  int tiles = d_filter_size/B;
  int residue = d_filter_size%B;
  cudaEventRecord(start);
  dim3 dimBlock(512);
  dim3 dimGrid((d_filter_size+dimBlock.x-1)/dimBlock.x);
  if (n<pow(2,27)){
    permFilterTillKernel<<<dimGrid, dimBlock>>>(d_bins_t,d_x,d_filter_time,n,
						B,ai,tiles,residue);
  }
  else{
    permFilterKernel<<<dimGrid, dimBlock>>>(d_bins_t,d_x,d_filter_time,n,B,ai,
					    d_filter_size);
  }
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_perm_filter,start, stop);

  if (TIMING) {
    *PF_T = (t_set_binst + t_perm_filter)/1.0e3 ;
    vprintf("inner_loop:Step 1.A (PERM + FILTER):------------------------ %lf\n", *PF_T);
  }

  /*CUFFT plan creation*/
  cudaEventRecord(start);
  cufftHandle plan;
  int batch = 1;
  cudaEventRecord(start);
  checkCudaErrors(cufftPlan1d(&plan,B,CUFFT_Z2Z,batch));

  /*CUFFT plan  execution*/
  checkCudaErrors(cufftExecZ2Z(plan,reinterpret_cast<cufftDoubleComplex *>(d_bins_t),
			       reinterpret_cast<cufftDoubleComplex *>(d_bins_f),CUFFT_FORWARD));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_CUFFT_ex,start, stop);

  samples_indices(d_bins_f,B_thresh,B,d_J,BC_T);

  *BC_T += t_CUFFT_ex/1.0e3;
  if (TIMING) {
    vprintf("innerLoop::perm_filter_cutoff:: Step 1.B (Sub-sampling and CUDAFFT)---------: %lf\n",*BC_T);
  }

  t_total =  *PF_T*1.0e3 + *BC_T*1.0e3;

  int print_t = 0;
  if (print_t){
    printf("innerLoop:perm_filter_cutoff::  time set bins on time domain =%lfs\n",t_set_binst/1.0e3);
    printf("innerLoop:perm_filter_cutoff::time perm+filter=%lfs \n", *PF_T );
    printf("innerLoop:perm_filter_cutoff:: time CUDAFFT execution =%lfs\n",t_CUFFT_ex/1.0e3);
    printf("innerLoop:perm_filter_cutoff::time perm_filter_cutoff=%lfs,time CUDAFFT + Cutoff = %lfs,\n\n",t_total/1.0e3,*BC_T);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));
  checkCudaErrors(cufftDestroy(plan));
  checkCudaErrors(cudaFree(d_bins_t));

  return 0;
}// end perm_filter_cutoff
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * This function implements the stage of  reversing the hash function for
 * recovery of the real indices of the largest frequency coefficients
 * as follows:
 * 1. For each  i = 0 to i = num-1, find indices that (1)map to J,
 *    i.e., lie within n/(2B) of (J * n/B) after permutation and
 *    (2) lie in d_J2 mod W_RLFC
 * 2. For each such index = j, increment d_vote[j] and if d_vote[j]
 *    is equal to  loop_threshold, increment I_F, and append
 *    to the vector d_I
 */
int reverse_hash(int *d_J,int n,int B_thresh,int B,int a,int ai,int loop_thresh,
		 int *d_vote, int *d_I,int *I_F,double *G_T,
		 int *d_J2,int num_RLFC,int W_RLFC,int RLFC_loops){

  /*timing variables */
  double t_initial_GT = 0.0;
  double t_initial = 0.0;
  double t_h_perm = 0.0;
  double t_qsort =  0.0;
  float  t_H_to_D_perm = 0.0;
  float  t_true_loc = 0.0;
  float t_total = 0.0;

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  t_initial_GT =  get_time();
  t_initial =  get_time();
  Pair *h_J2_sigma = (Pair *)malloc(num_RLFC*sizeof(*h_J2_sigma));
  for(int m =0; m< num_RLFC; m++){
    int prev = timesmod(d_J2[m], ai, W_RLFC);
    h_J2_sigma[m].first  = prev;
    h_J2_sigma[m].second = timesmod(prev, a, n);
  }
  t_h_perm = get_time() - t_initial;

  t_initial =  get_time();
  qsort(h_J2_sigma, num_RLFC, sizeof(h_J2_sigma[0]), comp_struct3);
  t_qsort = get_time() - t_initial;

  /*Allocate device  memory  for d_I_F*/
  int *d_I_F;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_I_F),sizeof(*d_I_F)));
  /*Copy host memory I_F to device memory d_I_F */
  checkCudaErrors(cudaMemcpy(d_I_F,I_F,sizeof(*d_I_F),cudaMemcpyHostToDevice));

  /*Allocate device  memory  for the device vector d_J2_sigma*/
  cudaEventRecord(start);
  Pair *d_J2_sigma;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_J2_sigma),num_RLFC*sizeof(*d_J2_sigma)));
  /* Copy host memory h_J2_sigma to device memory d_J2_sigma */
  checkCudaErrors(cudaMemcpy(d_J2_sigma, h_J2_sigma,num_RLFC*sizeof(*d_J2_sigma),cudaMemcpyHostToDevice));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_H_to_D_perm,start, stop);

  /* Compute the  intersection of d_J2_sigma and indices close to d_J * n/B,
     then invert to get true locations */
  dim3 dimBlock(512);
  dim3 dimGrid((B_thresh+dimBlock.x)/dimBlock.x);
  reverseHashKernel<<<dimGrid, dimBlock>>>(d_I,d_vote,d_J,d_J2_sigma,n,B_thresh,B,a,
					   loop_thresh,d_I_F,num_RLFC,W_RLFC);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_true_loc,start, stop);

  /*Copy device memory d_I_F  to host memory I_F */
  checkCudaErrors(cudaMemcpy(I_F,d_I_F,sizeof(*d_I_F),cudaMemcpyDeviceToHost));

  if (TIMING) {
    *G_T = get_time()-t_initial_GT;
    vprintf("Step 1.D (REVERSE HASH):----------------------------- %lf\n\n", get_time()-t_initial_GT);
    vprintf("#####################################################################\n\n");
  }

  t_total = t_h_perm  + t_qsort + (t_H_to_D_perm +  t_true_loc)/1.0e3;
  int print_t = 0;
  if (print_t){
    printf("innerLoop:reverse_hash: time h_permute fillin = %lfs,  time quick sort =  %lfs \n",
	   t_h_perm, t_qsort);
    printf("innerLoop:reverse_hash: time HtoD permute = %lfs,time true locations = %lfs \n",
	   t_H_to_D_perm/1.0e3,t_true_loc/1.0e3);
    printf("innerLoop:reverse_hash: total time = %lfs \n", t_total);
    printf("innerLoop:reverse_hash:time REVERSE HASH = %lfs \n", *G_T );
  }

  /*cleanup memory*/
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));
  checkCudaErrors(cudaFree(d_J2_sigma));

  return 0;
}// end reverse_hash
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * d_I contains the indices that we want to estimate.
 * d_bins_f contains a B-dimensional array for each of the `loops`
 * iterations of the outer loop.  Every coordinate i of x "hashes to" a
 * corresponding coordinate (permute[j] * i) mod B of d_bins_f[j], which
 * gives an estimate of x[i].
 * We estimate each coordinate as the median (independently in real and
 * imaginary axes) of its images in the rows of x_samp.
 */
Node *estimate_values(int *d_I,int *I_F,complex_t *d_bins_f,int loops,int n,
     	              int *h_permute,int  B,complex_t *d_filter_freq,int loops_loc){

  int my_I_F = *I_F;

  /*timing variables */
  float  t_H_to_D_permute = 0.0;
  float  t_estimate = 0.0;
  float  t_H_to_D_hat_x = 0.0;
  float t_total = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device  memory  for the vector d_permute and copy host memory to device memory*/
  cudaEventRecord(start);
  int *d_permute;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_permute),loops*sizeof(*d_permute)));
  checkCudaErrors(cudaMemcpy(d_permute,h_permute,loops*sizeof(*d_permute),cudaMemcpyHostToDevice));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_H_to_D_permute,start, stop);

  /*Allocate device  memory  for the  d_hat_x vector*/
  Node *d_hat_x;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_hat_x),sizeof(*d_hat_x)*my_I_F));

  /* Compute the  estimated values */
  dim3 dimBlockEst(512);
  int gridDimxEst = (my_I_F+dimBlockEst.x-1)/dimBlockEst.x;
  dim3 dimGridEst(gridDimxEst);
  cudaEventRecord(start);
  estimateValuesKernel<<<dimGridEst,dimBlockEst>>>(d_hat_x,d_I,my_I_F,d_bins_f,
						   loops,n,d_permute,B,d_filter_freq);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_estimate,start, stop);

  /*Allocate host memory for the h_hat_x vector and copy device memory to host memory*/
  cudaEventRecord(start);
  Node *h_hat_x = (Node *) malloc(my_I_F*sizeof(*h_hat_x));
  checkCudaErrors(cudaMemcpy(h_hat_x,d_hat_x,my_I_F*sizeof(*d_hat_x),cudaMemcpyDeviceToHost));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_H_to_D_hat_x,start, stop);

  /*timing*/
  t_total =  (t_H_to_D_permute + t_estimate + t_H_to_D_hat_x)/1.0e3;
  int print_t = 0;
  if (print_t){
    printf("innerLoop:estimate_values: time HtoD permute=%lfs\n",t_H_to_D_permute/1.0e3);
    printf("innerLoop:estimate_values: time estimate values kernel=%lfs,time HtoD hat_x=%lfs\n",t_estimate/1.0e3,t_H_to_D_hat_x/1.0e3);
    printf("innerLoop:estimate_values: total time estimate values=%lfs\n",t_total);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));
  checkCudaErrors(cudaFree(d_permute));
  checkCudaErrors(cudaFree(d_hat_x));

  int print_hat_x = 0;
  if (print_hat_x){
    for (int i=0;i<my_I_F;i+=16){
      printf("inner_loop:estimate_values h_hat_x[%d].value.x = %lf, h_hat_x[%d].value.y = %lf,i = %d\n",i,h_hat_x[i].value.x,i,h_hat_x[i].value.y,i);
    }
  }

  return h_hat_x;
}// end estimate_values
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Compute the samples vector as the  square of the
 * absolute value of a device sampling vector.
 * Stages:
 * 1. Allocate device memory for the device vectors  
 *    d_samples and d_samples_u to hold the magnitudes  of
 *    the components of the sampling vector
 *    d_samp. 
 * 2. Compute  the d_samples and d_samples_u   vectors.
 * 3. Allocate device memory for the host vector
 *     h_samples to copy  the d_samples vector.  
 * 4. Wrap the h_samples vector with a  thrust pointer 
 * 5. Start cutoff stage.
 *
 */
int samples_indices(complex_t *d_samp,int B,int samp_size,int* d_output,double *CTime){

  /*timing variables */
  float t_fill_samples_kernel = 0.0;
  float t_fill_samples = 0.0;
  float t_index = 0.0;
  float t_total = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device memory for the vector d_samples to be sorted */
  cudaEventRecord(start);
  double *d_samples;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_samples), sizeof(*d_samples)*samp_size));
  /*Allocate unified memory for the vector d_samples_u to compute largest indices*/
  double *d_samples_u;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_samples_u), sizeof(*d_samples_u)*samp_size));

  /*compute  d_samples and d_samples_u */
  dim3 dimBlock(256);
  dim3 dimGrid(samp_size/dimBlock.x);
  samplesKernel<<<dimGrid, dimBlock>>>(d_samp,d_samples,d_samples_u,samp_size);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_fill_samples_kernel,start, stop);

  /* wrap d_samples with a thrust device pointer*/
  cudaEventRecord(start);
  thrust::device_ptr<double> d_samples_ptr(d_samples);

  /*cutoff*/
  cutoff(d_output,B,d_samples_ptr,d_samples_u,samp_size);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_index,start, stop);

  t_fill_samples = t_fill_samples_kernel;
  t_total =  t_fill_samples + t_index;
  *CTime += t_total/1.0e3;
  int print_t = 0;
  if (print_t){
    printf("\ninnerLoop::samples_indices:time fill samples kernel = %lfs,time fill samples = %lfs \n",
	   t_fill_samples_kernel/1.0e3,t_fill_samples/1.0e3);
    printf("innerLoop::samples_indices:time indices = %lfs, time (fill samples + indices) =  %lf \n", t_index/1.0e3,t_total/1.0e3);
  }

  /*cleanup memory*/
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));
  checkCudaErrors(cudaFree(d_samples));
  checkCudaErrors(cudaFree(d_samples_u));

  return 0;
}// end samples_indices
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
 * Output the indices corresponding to the num largest elements of samples.
 * Output is sorted.
 * int *output -- output array
 * int B -- size of output, equals to the cutoff number
 * real_t *samples -- input array
 * int samp_size -- size of d_samples array
 * Cut B into num
 */
void cutoff(int *d_output,int B, thrust::device_ptr<double> &d_samples_ptr, double *d_samples,int samp_size){

  float  t_fst_loop = 0.0;
  double  t_sec_loop = 0.0;
  float  t_sort = 0.0;
  float  t_total = 0.0;
  double nth_time;
  double initial_time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  assert(samp_size >= B + 1);

  initial_time = get_time();
  double *d_cutoff;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_cutoff),sizeof(*d_cutoff)));
  *d_cutoff = nth_element_immutable_cuda(d_samples_ptr,samp_size,samp_size-B-1);
  nth_time= get_time()-initial_time;

  /* allocate device memory for d_count and set to zero*/
  int *d_count;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_count),sizeof(*d_count)));
  checkCudaErrors(cudaMemset(d_count,0,sizeof(*d_count)));
  checkCudaErrors(cudaDeviceSynchronize());

  /*compute d_output  with first for loop*/
  cudaEventRecord(start);
  dim3 dimBlock1(512);
  dim3 dimGrid1((samp_size+dimBlock1.x-1)/dimBlock1.x);
  cutoffKernel<<<dimGrid1,dimBlock1>>>(d_output,d_cutoff,d_samples,samp_size,d_count,B);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_fst_loop,start, stop);

  /*copy device memory d_count to host memory h_count*/
  int *h_count = (int*)malloc(sizeof(*h_count));
  checkCudaErrors(cudaMemcpy(h_count,d_count,sizeof(*d_count),cudaMemcpyDeviceToHost));
  initial_time = get_time();
  if (*h_count < B){
    /*compute d_output  with second for loop*/
    for(int i = 0; i < samp_size && *h_count < B; i++){
      if (d_samples[i] == *d_cutoff) {
	d_output[(*h_count)++] = i;
      }
    }
    t_sec_loop = get_time()-initial_time;
  }
  assert(*h_count == B);

  /*wrap and sort d_output,extract d_output*/
  cudaEventRecord(start);
  thrust::device_ptr<int> d_output_ptr(d_output);
  thrust::stable_sort(d_output_ptr,d_output_ptr+*h_count);
  d_output = thrust::raw_pointer_cast(d_output_ptr);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_sort,start, stop);
  
  t_total = nth_time*1.0e3 + t_fst_loop + t_sec_loop*1.0e3 + t_sort;
  int print_t = 0;
  if(print_t){printf("innerLoop::cutoff: nth_time=%lfs,time loop 1=%lfs,time loop 2=%lfs,time sort=%lfs,total time=%lfs\n",
		     nth_time,t_fst_loop/1.0e3,t_sec_loop,t_sort/1.0e3,t_total/1.0e3);}

}// end cutoff
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
 * This function returns the num'th smallest element of the length n input.
 * The input is a thrust pointer to a vector of length n.
 * 
 */
double nth_element_immutable_cuda(thrust::device_ptr<double> &d_samples_ptr, int n, int num){

  /*time measurement variables*/
  float t_sort = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*sort the device pointer*/
  cudaEventRecord(start);
  thrust::stable_sort(d_samples_ptr, d_samples_ptr+n);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_sort,start, stop);

  //printf("Sort time   %4.2fms\n", t_sort);

  return d_samples_ptr[num];
}//end nth_element_immutable_cuda
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* computes the mod n of a product of integers */
int timesmod(const int x, const int a, const int n) {
  return (int)((((long long int)(x))*a)%n);
}//end timesmod

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a device vector d_x */
__global__
void printComDxKernel(complex_t *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    printf("innerLoop:PrintComDxKernel::d_x[%d].x = %lf \n", i, d_x[i].x);
    printf("innerLoop:PrintComDxKernel::d_x[%d].y = %lf \n", i, d_x[i].y);
  }
}//end printComDxKernel

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a int device vector d_x */
__global__
void printIntDxKernel(int *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    printf("innerLoop:PrintIntDxKernel::d_x[%d] = %d \n", i, d_x[i]);
  }
}//end printIntDxKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a double device vector d_x */
__global__
void printDouDxKernel(double *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    printf("innerLoop:PrintDouDxKernel::d_x[%d] = %6.25lf \n", i, d_x[i]);
  }
}//end printDouDxKernel
////////////////////////////////////////////////////////////////////////////////
////////////// Kernel for restr_loc_freq_coeff function.
////////////////////////////////////////////////////////////////////////////////
/*Filtering  and storing the input signal*/
__global__
void restLocFreqCoeffKernel(complex_t *d_x_filt_t, complex_t *d_x,
			    int W_RLFC,int offset,int sigma){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < W_RLFC){
    d_x_filt_t[i].x = d_x[offset + i*sigma].x;
    d_x_filt_t[i].y = d_x[offset + i*sigma].y;
  }
}//end restLocFreqCoeffKernel
////////////////////////////////////////////////////////////////////////////////
//////////////Kernel for the samples_indices function
////////////////////////////////////////////////////////////////////////////////
/*Computes the samples vectors needed  to compute largest indices*/
__global__
void samplesKernel(complex_t *d_samp,double *d_samples,double *d_samples_u,
		   int samp_size){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < samp_size){
    d_samples[i] = pow(d_samp[i].x,2) + pow(d_samp[i].y,2);
    d_samples_u[i] = d_samples[i];
  }
}// end samplesKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////Kernels for the perm_filter_cutoff function
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*permute,filter and bin  the input signal d_x */
__global__
void permFilterTillKernel(complex_t *d_bins_t,complex_t *d_x,complex_t *d_filter_t,
			  int n,int B,int ai,int tiles,int residue){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < B || i == B){

    if(i < B){

      for(int j=0; j<tiles; j++){

	int id = j*B + i;
	int index = (id*ai)&(n-1);

	d_bins_t[i].x += d_x[index].x*d_filter_t[id].x - d_x[index].y*d_filter_t[id].y;
	d_bins_t[i].y += d_x[index].x*d_filter_t[id].y + d_x[index].y*d_filter_t[id].x;
      }
    }

    if(i == B){

      for(int j=0; j<residue; j++){

	int id = j*tiles + i;
	int index = (id*ai)&(n-1);

	d_bins_t[j].x += d_x[index].x*d_filter_t[id].x - d_x[index].y*d_filter_t[id].y;
	d_bins_t[j].y += d_x[index].x*d_filter_t[id].y + d_x[index].y*d_filter_t[id].x;
      }
    }
  }
}//end permFilterTillKernel

////////////////////////////////////////////////////////////////////////////////
/*permute and filter the input signal d_x */
__global__
void permFilterKernel(complex_t *d_bins_t, complex_t *d_x, complex_t *d_filter_t,
		      int n,int B,int ai,int d_filter_size){

  int i = threadIdx.x + blockIdx.x * blockDim.x;
  if(i < d_filter_size){

    int i_mod_B = i&(B-1);
    int index = (i*ai)&(n-1);
    d_bins_t[i_mod_B].x  += d_x[index].x*d_filter_t[i].x - d_x[index].y*d_filter_t[i].y;
    d_bins_t[i_mod_B].y  += d_x[index].x*d_filter_t[i].y + d_x[index].y*d_filter_t[i].x;
  }
}//end permFilterKernel
////////////////////////////////////////////////////////////////////////////////
//////////////Kernel for the reverse_hash function
////////////////////////////////////////////////////////////////////////////////

/*
 *  Compute the  intersection of d_J2_sigma and indices close to J * n/B,
 *  then invert to get true locations.
 * 
 */
__global__
void reverseHashKernel(int *d_I,int *d_vote,int *d_J,Pair *d_J2_sigma,int n,
		       int B_thresh,int B, int a,int loop_threshold,int *d_I_F,
		       int num_RLFC, int W_RLFC){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i <  B_thresh){
    int low, high;
    low = ((int)(ceil((d_J[i] - 0.5) * n / B)) + n)&(n-1);
    high = ((int)(ceil((d_J[i] + 0.5) * n / B)) + n)&(n-1);
    int key = low&(W_RLFC-1);
    int index = findUpperBound(d_J2_sigma, num_RLFC, key);
    int location = low - (low&(W_RLFC-1));
    int locinv = timesMod(location, a, n);
    for(int j = index; ; j++){
      if (j == num_RLFC){
	j -= num_RLFC;
	location = (location + W_RLFC)&(n-1);
	locinv = timesMod(location, a, n);
      }
      int approved_loc = location + d_J2_sigma[j].first;
      if((low < high && (approved_loc >= high || approved_loc < low)) || (low > high && (approved_loc >= high && approved_loc < low))){
	break;
      }
      int loc = (locinv + d_J2_sigma[j].second)&(n-1);

      atomicAdd(&d_vote[loc],1);

      if(d_vote[loc]==loop_threshold){
	d_I[atomicAdd(d_I_F,1)]=loc;
      }
    }
  }
}//end reverseHashKernel
////////////////////////////////////////////////////////////////////////////////
//////////////Kernel for estimate_values function
////////////////////////////////////////////////////////////////////////////////
/*
 *  Compute the  estimate values of the largest indices contained in d_I.
 *
 */
__global__
void estimateValuesKernel(Node *d_hat_x,int *d_I,int my_I_F,
			  complex_t *d_bins_f,int loops,int n,
			  int *d_permute,int B,complex_t *d_filter_f){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i <  my_I_F){
    double *values_real = (double*)malloc(loops*sizeof(*values_real));
    double *values_img  = (double*)malloc(loops*sizeof(*values_img));

    int position = 0;
    for(int j = 0; j < loops; j++){
      int permuted_index= timesMod(d_permute[j], d_I[i],n);
      int hashed_to = permuted_index / (n/B);
      int dist = permuted_index & ((n/B)-1);
      if (dist > (n/B)/2) {
	hashed_to = (hashed_to + 1)&(B-1);
	dist -= n/B;
      }
      dist = (n - dist) & (n-1);
      complex_t filter_value = d_filter_f[dist];
      int k = j*B + hashed_to;
      double abs2 = pow(filter_value.x,2) + pow(filter_value.y,2);
      double mult_real = d_bins_f[k].x*filter_value.x + d_bins_f[k].y*filter_value.y;
      double mult_img  = d_bins_f[k].y*filter_value.x - d_bins_f[k].x*filter_value.y;
      values_real[position] = mult_real/abs2;
      values_img[position] =  mult_img/abs2;
      position++;
    }
    int location = (loops - 1) / 2;

    insertionSortDev(values_real,loops);
    insertionSortDev(values_img,loops);

    double realv = values_real[location];
    double imgv = values_img[location];
    d_hat_x[i].key = d_I[i];
    d_hat_x[i].value.x = realv;
    d_hat_x[i].value.y = imgv;
    /*clean memory*/
    free(values_real);
    free(values_img);

  }
}//end estimateValuesKernel
/////////////////////////////////////////////////////////////////////////////
///////////// Kernel for cutoff function
////////////////////////////////////////////////////////////////////////////////
/*Output the indices corresponding to the num largest elements of samples. */
__global__
void cutoffKernel(int *d_output,double *d_cutoff,double *d_samples,int samp_size,
		  int *d_count,int B){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < samp_size){
    if (d_samples[i]>*d_cutoff && *d_count < B){
      d_output[atomicAdd(d_count,1)] = i;
      int m  = 0;
      int ii = 100;
      if(blockIdx.x<m && i<ii){printf("threadIdx.x=%d,blockIdx.x=%d,i=%d,B=%d,d_output[%d] = %d,*d_count-1=%d,d_samples[%d]=%6.25lf,*d_cutoff=%6.25lf\n",
				      threadIdx.x,blockIdx.x,i,B,*d_count-1,d_output[i],*d_count-1,i,d_samples[i],*d_cutoff);}
    }
  }
}//end cutoffKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* computes  a product of integers mod n*/
__device__
int timesMod(const int x, const int a, const int n) {
  return (int)((((long long int)(x))*a)%n);
}//end timesMod
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* binary search to find the upper bound of an array*/
__device__
int findUpperBound(Pair *d_J, int n, int val){
  int mid, low = 0;
  int high = n - 1;

  if(val >= d_J[high].first){
    return high;
  }

  mid = (low + high) / 2;
  while (high > low) {
    if (d_J[mid].first >=val)
      high  = mid;
    else
      low = mid + 1;

    mid = (low + high) / 2;
  }
  return mid;
}//end findUpperBound
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*returns a sorted array  */
__device__
void insertionSortDev(double *d_x, int n) {
  double key;
  int j;
  for (int i = 1; i < n; i++) {
    key = d_x[i];
    j = i - 1;
    while (j >= 0 && d_x[j] > key) {
      d_x[j + 1] = d_x[j];
      j = j - 1;
    }
    d_x[j + 1] = key;
  }
}//end insertionSortDev
