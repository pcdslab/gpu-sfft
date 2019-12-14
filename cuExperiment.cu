/* 
 *  Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *  Eric Price, Massachusetts Institute of Technology.
 *
 */
/*
 *  Copyright (C) 2019 Oswaldo Artiles and Fahad Saeed
 *  Florida International University, Florida, USA.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
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
 
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

// includes CUDA project
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>


extern "C"{
      #include "utils.h"
      #include "experiment.h"
      #include "filters.h"
      #include "cudaFft.h"
      #include "timer.h"
      #include "outerLoop.h"
}
////////////////////////prototype functions/////////////////////////////////////
__global__ void printComDxExKernel(complex_t *d_x,int n);
__global__ void makeDolphchebyshevKernel1(complex_t *d_filter_time,int w,double t0,
					  double tolerance);
__global__ void makeDolphchebyshevKernel2(complex_t *d_filter_time,int w);
__global__ void makeMultipleKernel1(complex_t *d_filter_freq,double *d_max,int n);
__global__ void makeMultipleKernel2(complex_t *d_filter_time,int w,int n);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
__device__ double cuCheb(double m, double x);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* Compute the input signal and the DFT of the input signal.
 * Compute the filter in both time and frequency domain. 
 * Starts running the experiment.
 */

  int cudafft_experiment(int n,double lobefrac,double tolerance,int b_f,
			 int B,int B_thresh,int loops_loc,int loops_thresh,
			 int loops_est,int W_RLFC,int RLFC_loops,
			 int repetitions,int k,double std_noise){

  float t_H_to_D;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*The host input signal vector*/
  complex_t *h_x=(complex_t *)malloc(n*sizeof(*h_x));

  /*The host DFT output signal vector */
  complex_t *x_f = (complex_t *)calloc(n, sizeof(*x_f));

  /* indices of the large frequencies in the DFT output signal vector x_f */
  int *LARGE_FREQ = (int *)malloc(k*sizeof(*LARGE_FREQ));

  /*Initialize the large bins in freq domain */
  for (int i=0; i<k; i++){
    LARGE_FREQ[i] = (unsigned)floor(drand48() * n);
    x_f[LARGE_FREQ[i]].x = 1.0;
    x_f[LARGE_FREQ[i]].y = 0.0;
  }
  
  /*Reverse cuFFT to generate input signal, h_x, from DFT output signal, x_f. */
  int forward = 0;
  int print_y = 0;
  cuda_fft_dft_d(x_f,h_x,n,forward,print_y);
  
  /*Add  noise to the input signal*/
  double snr_achieved;
  snr_achieved = AWGN(h_x,n,std_noise);
  if(std_noise != 0){
    printf("SNR = %g / %.2f dB \n\n", snr_achieved, 10 * log10(snr_achieved));
  }

  /*Forward cuFFT to generate the host DFT output signal  x_f from input signal h_x */
  forward = 1;
  print_y = 0;
  cuda_fft_dft_d(h_x,x_f,n,forward,print_y);
  
  /*normalization of DFT output signal x_f */
  for(int i = 0; i < n; i++){
    x_f[i].x /= n;
    x_f[i].y /= n;
  }

  /*Allocate device  memory  for the device input signal d_x*/
  complex_t *d_x;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_x), sizeof(cufftDoubleComplex)*n));
  
  /* Copy host memory to device memory*/
  cudaEventRecord(start);
  checkCudaErrors(cudaMemcpy(d_x, h_x, n*sizeof(complex_t),cudaMemcpyHostToDevice));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_H_to_D,start, stop);
  
  /*compute the filter with window size filter = w */
  printf("\n computing filter starts\n" );
  make_multiple_t(lobefrac,tolerance,n,b_f,d_x,B,B_thresh,loops_loc,
		  loops_thresh,loops_est,W_RLFC,RLFC_loops,repetitions,
		  LARGE_FREQ,k,x_f,t_H_to_D);

  /*cleanup memory*/
  checkCudaErrors(cudaFree(d_x));
  free(h_x);
  free(x_f);
  free(LARGE_FREQ);

  return 0;
}// end cudafft_experiment
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
 * Run an experiment.  Parameters are:
 * d_x:            the device input signal, 
 * n:              the length of d_x
 * d_filter_time:  the time component of the filter
 * d_filter_freq   the frequency component of the filter
 * d_filter_size   the size of the filter
 * B:              number of samples in subsampling during location
 * B_thresh:       number of samples considered "heavy"
 * loops_loc:      number of location loops
 * loops_thresh:   number of times a coordinate must seem heavy.
 * loops_est:      number of estimation loops
 * repetitions:    repeat the experiment this many times for timing
 * LARGE_FREQ:     locations of largest coefficients.
 * k:              number of HHs, used for evaluation only.
 * x_f:            true DFT computed with the CUDAFFTW package.
 * t_H_to_D:       time for memory copy of host input signal h_x to 
 *                 device input signal d_x
 */
void run_experiment(complex_t *d_x,int n,complex_t *d_filter_time,
		    complex_t *d_filter_freq,int d_filter_size,int B,
		    int B_thresh,int loops_loc,int loops_thresh,
		    int loops_est,int W_RLFC,int RLFC_loops,
		    int repetitions,int *LARGE_FREQ,int k,
		    complex_t *x_f,float t_H_to_D){

  printf("GPU-SFFT Results\n");
  printf("******************************************************************************\n");
  //////////////////////////////////run_experiment///////////////////////////////////////
  Node *hat_x = NULL;
  int I_F = 0;  //the number of indices of  largest frequencies coefficients
  int loops = loops_loc + loops_est;
  
  for(int i = 0; i < repetitions; i++){
    reset_timer();
    hat_x = outer_loop(d_x,n,d_filter_time,d_filter_freq,d_filter_size,B,B_thresh,
		     W_RLFC,RLFC_loops,loops_thresh,loops_loc,loops,&I_F, t_H_to_D);

  }
  //////////////////////////////////Mean Absolute Error Computation///////////////////////////////////////
  
  /* computation of the error of the SFFT computed with GPU-SFFT as
     compared to the DFT computed with cuFFT (NVIDIA-CUDA) library*/

  double initial_time = get_time();
  complex_t *hat_x_Large = (complex_t *)calloc(n,sizeof(*hat_x_Large));

  /*sort the hat_x according to the cabs value of second element*/
  qsort(hat_x, I_F, sizeof(*hat_x), comp_struct);
 
  for(int i=0; i< k; i++){
    hat_x_Large[hat_x[i].key].x = hat_x[i].value.x;
    hat_x_Large[hat_x[i].key].y = hat_x[i].value.y;
  }
  double t_hat_x_large = get_time() - initial_time;

  int large_found = 0;
  int FOUND=0;
 
  qsort(hat_x, I_F, sizeof(*hat_x), comp_struct2);

  for(int i = 0; i < k; i++){
    Node *item = (Node*)bsearch(&LARGE_FREQ[i], hat_x, I_F, sizeof(*hat_x), comp_struct2);
    if(item != NULL ){
       FOUND++;
    }
     large_found += (hat_x_Large[LARGE_FREQ[i]].x != 0.0 || hat_x_Large[LARGE_FREQ[i]].y != 0.0);
  }

  /* compute Mean Absolute Error (MAE) */
  double ERROR = 0.0;
  for(int i=0; i< n ; i++){
    ERROR += cabs_t2(complexSubstr(hat_x_Large[i],x_f[i]));
  }
  double MSE = ERROR/k;
  free(hat_x_Large);

  printf("ERROR:\n");
  printf("K=%d; MISSED (estimation, result) = (%d, %d), MSE*k = %lg, Mean Square Error (MSE) = %lg\n",k,
	 k-FOUND, k-large_found, ERROR, MSE);
  printf("time to extract the values of the largest frequencies: %lfs \n",t_hat_x_large);
  printf("******************************************************************************\n\n");

} // end run_experiment
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* Compute the times for the computation of the cuFFT of the input signal.*/

int cudafft_results (complex_t *d_x,int n,int forward,int print_y){

  /*Allocate memory for the host DFT output signal vector */
  complex_t *x_ft = (complex_t *)malloc(n*sizeof(*x_ft));
  /*Allocate memory  for the host input signal vector*/
  complex_t *h_x=(complex_t *)malloc(n*sizeof(*h_x));

  /* Copy device memory to host memory*/
  checkCudaErrors(cudaMemcpy(h_x,d_x,n*sizeof(complex_t),cudaMemcpyDeviceToHost));

  /*cuFFT to generate the DFT output signal  x_ft from input signal x */
  cuda_fft_dft_d(h_x,x_ft,n,forward,print_y);

  free(x_ft);
  free(h_x);

  return 0;
}// end cudafft_results
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void make_multiple_t(double lobefrac,double tolerance,int n,int b_f,
		     complex_t *d_x,int B,int B_thresh,int loops_loc,
		     int loops_thresh,int loops_est,int W_RLFC,int RLFC_loops,
		     int repetitions,int *LARGE_FREQ,int k,complex_t *x_f,
		     float t_H_to_D){

  /*timing variables */
  float  t_comp_filt_time = 0.0;
  float  t_CUDAFFT_time = 0.0;
  double  t_shift_time = 0.0;
  float  t_set_zero_time = 0.0;
  double t_total_filt_time = 0.0;
  float  t_set_zero_dg = 0.0;
  float  t_copy_dg = 0.0;
  float  t_CUDAFFT_g = 0.0;
  double  t_comp_s = 0.0;
  double t_filtf_ds_dmax = 0.0;
  float  t_norm_filtf = 0.0;
  double t_comp_filtf = 0.0;
  float  t_CUDAFFT_g_filtf = 0.0;
  float  t_copy_dg_filt = 0.0;
  float  t_norm_filtt = 0.0;
  float  t_total = 0.0;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  double initial_time_f = get_time();

  /*compute the  d_filter_time*/
  double initial_filt_time = get_time();
  int w = (int)((1 / M_PI) * (1/lobefrac) * acosh(1./tolerance));
  if (!(w % 2)){
    w--;
  }
  assert(w <= n);
  double t0 = cosh(acosh(1/tolerance) / (w-1));

  /*allocate device memory for d_filter_time, size=w */
  complex_t *d_filter_time;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_filter_time),w*sizeof(*d_filter_time)));
  checkCudaErrors(cudaDeviceSynchronize());

  /*compute d_filter_time*/
  cudaEventRecord(start);
  dim3 dimBlockt1(512);
  dim3 dimGridt1((w+dimBlockt1.x-1)/dimBlockt1.x);
  makeDolphchebyshevKernel1<<<dimGridt1,dimBlockt1>>>(d_filter_time,w,t0,tolerance);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_comp_filt_time,start, stop);

  /*CUFFT plan (time) creation*/
  cudaEventRecord(start);
  cufftHandle plan_t;
  int batch = 1;
  checkCudaErrors(cufftPlan1d(&plan_t,w,CUFFT_Z2Z,batch));
  /*CUFFT plan forward execution*/
  checkCudaErrors(cufftExecZ2Z(plan_t,reinterpret_cast<cufftDoubleComplex *>(d_filter_time),
			       reinterpret_cast<cufftDoubleComplex *>(d_filter_time),CUFFT_FORWARD));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_CUDAFFT_time,start, stop);

  /*shift d_filter_time*/
  double initial_time = get_time();
  cuShift(d_filter_time,w,w/2);
  t_shift_time  = get_time()-initial_time;

  /*set to zero the imaginary components of  d_filter.time*/
  cudaEventRecord(start);
  dim3 dimBlockt2(512);
  dim3 dimGridt2((w+dimBlockt2.x-1)/dimBlockt2.x);
  makeDolphchebyshevKernel2<<<dimGridt2,dimBlockt2>>>(d_filter_time,w);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_set_zero_time,start, stop);
  double t_filt_time = get_time()-initial_filt_time;

  /*allocate unified memory for d_filter_freq, size = n */
  complex_t *d_filter_freq;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_filter_freq),n*sizeof(*d_filter_freq)));
  checkCudaErrors(cudaDeviceSynchronize());

  /*allocate unified memory for d_g, size =n*/
  cudaEventRecord(start);
  complex_t *d_g;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_g),n*sizeof(*d_g)));
  /*set to zero the components of  d_g*/
  checkCudaErrors(cudaMemset(d_g,0.0,sizeof(*d_g)*n));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_set_zero_dg,start, stop);

  /*copy device memory d_filter_time+w/2 to device memory d_g,size=w/2*/
  cudaEventRecord(start);
  checkCudaErrors(cudaMemcpy(d_g,d_filter_time+w/2,(w/2)*sizeof(*d_g),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaDeviceSynchronize());

  /*copy device memory d_filter.time to device memory d_g+n-(w/2),size=w/2*/
  checkCudaErrors(cudaMemcpy(d_g+n-(w/2),d_filter_time,(w/2)*sizeof(*d_g),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_copy_dg,start, stop);

  /*CUFFT plan (frequency) creation,size=n*/
  cudaEventRecord(start);
  cufftHandle plan_f;
  batch = 1;
  checkCudaErrors(cufftPlan1d(&plan_f,n,CUFFT_Z2Z,batch));
  /*CUFFT forward plan  execution over d_g*/
  checkCudaErrors(cufftExecZ2Z(plan_f,reinterpret_cast<cufftDoubleComplex *>(d_g),
			       reinterpret_cast<cufftDoubleComplex *>(d_g),CUFFT_FORWARD));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_CUDAFFT_g,start, stop);

  /*compute s */
  initial_time = get_time();
  complex_t s;
  s.x = 0.0;
  s.y = 0.0;
  /*compute d_s*/
  for(int i = 0; i < b_f; i++){
    s.x += d_g[i].x;
    s.y += d_g[i].y;
  }
  t_comp_s = get_time() - initial_time;

  /*allocate unified memory for d_max and set to zero*/
  initial_time = get_time();
  double *d_max;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_max),sizeof(*d_max)));
  checkCudaErrors(cudaMemset(d_max,0.0,sizeof(*d_max)));
  /*compute  d_filter_freq,*d_max and s*/
  int offset = b_f/2;
  for(int i = 0; i < n; i++){
    d_filter_freq[(i+n+offset)&(n-1)].x = s.x;
    d_filter_freq[(i+n+offset)&(n-1)].y = s.y;
    *d_max = *d_max>cabs_t(s)?*d_max:cabs_t(s);
    s.x +=  d_g[(i + b_f)&(n-1)].x - d_g[i].x;
    s.y +=  d_g[(i + b_f)&(n-1)].y - d_g[i].y;
  }
  t_filtf_ds_dmax = get_time()-initial_time;

  /*normalize d_filter_freq*/
  cudaEventRecord(start);
  dim3 dimBlock1(512);
  dim3 dimGrid1((n+dimBlock1.x-1)/dimBlock1.x);
  makeMultipleKernel1<<<dimGrid1,dimBlock1>>> (d_filter_freq,d_max,n);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_norm_filtf,start, stop);

  /*compute d_filter_freq with offset*/
  initial_time = get_time();
  complex_t offsetc;
  offsetc.x = 1.0;
  offsetc.y = 0.0;
  complex_t step;
  step.x = cos(-2*M_PI*(w/2)/n);
  step.y = sin(-2*M_PI*(w/2)/n);
  for(int i = 0; i < n; i++){
    d_filter_freq[i] = complexMul(d_filter_freq[i],offsetc);
    offsetc = complexMul(offsetc,step);
  }
  t_comp_filtf = get_time()-initial_time;

  /*CUFFT inverse plan  execution to compute d_g as the inverse DFT of d_filter_freq*/
  cudaEventRecord(start);
  checkCudaErrors(cufftExecZ2Z(plan_f,reinterpret_cast<cufftDoubleComplex *>(d_filter_freq),
			       reinterpret_cast<cufftDoubleComplex *>(d_g),CUFFT_INVERSE));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_CUDAFFT_g_filtf,start, stop);

  /*copy device memory d_g  to device memory d_filter_time,size=w*/
  cudaEventRecord(start);
  checkCudaErrors(cudaMemcpy(d_filter_time,d_g,w*sizeof(*d_g),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_copy_dg_filt,start, stop);

  /*normalize d_filter_time,size=w */
  dim3 dimBlock2(512);
  dim3 dimGrid2((w+dimBlock2.x-1)/dimBlock2.x);
  makeMultipleKernel2<<<dimGrid2,dimBlock2>>> (d_filter_time,w,n);
  checkCudaErrors(cudaDeviceSynchronize());
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&t_norm_filtt,start, stop);
  double t_filter = get_time() - initial_time_f;

  /*timing filter_time*/
  t_total_filt_time = t_comp_filt_time  + t_CUDAFFT_time + t_shift_time*1.0e3 + t_set_zero_time;
  /*timing filter_freq + filter_time*/
      t_total = t_filt_time*1.0e3 + t_set_zero_dg + t_copy_dg + t_CUDAFFT_g + t_comp_s*1.0e3
	                    + t_filtf_ds_dmax*1.0e3 + t_norm_filtf + t_comp_filtf*1.0e3 + t_CUDAFFT_g_filtf
	+ t_copy_dg_filt + t_norm_filtt;

      int print_t = 0;
      if (print_t){
	printf("make_dolphchebyshev_t: time filter computation=%lfs,time CUDAFFT=%lfs,time shift=%lfs,time set zero=%lf,total time=%lfs \n",
	       t_comp_filt_time/1.0e3,t_CUDAFFT_time/1.0e3,t_shift_time,t_set_zero_time/1.0e3,t_total_filt_time/1.0e3);
	printf("make_multiple_t: time computation filter_time=%lfs,time set zero d_g=%lfs,time copy d_g=%lfs,time CUDAFFT d_g=%lfs\n",
	       t_filt_time,t_set_zero_dg/1.0e3,t_copy_dg/1.0e3,t_CUDAFFT_g/1.0e3);
	printf("make_multiple_t: time compute d_s=%lfs,time computation filt freq,d_max and d_s=%lfs,time normalization filter freq=%lfs,\n",
	       t_comp_s,t_filtf_ds_dmax,t_norm_filtf/1.0e3);
	printf("make_multiple_t:time computation filter freq=%lfs,time CUDAFFT d_g filter freq=%lfs,time copy dg filt time =%lfs\n",
	       t_comp_filtf,t_CUDAFFT_g_filtf/1.0e3,t_copy_dg_filt/1.0e3);
	printf("make_multiple_t:normalization filter time=%lfs,total time computation filter=%lfs\n",t_norm_filtt/1.0e3,t_total/1.0e3);
      }

      int d_filter_size = w;
      printf(" Window size: Filter : %d\n", w);
      printf(" Time for computing  filter:  %lfs\n", t_filter );
      printf("******************************************************************************\n\n");

      /*cleanup memory*/
      checkCudaErrors(cudaEventDestroy(start));
      checkCudaErrors(cudaEventDestroy(stop));
      checkCudaErrors(cufftDestroy(plan_t));
      checkCudaErrors(cufftDestroy(plan_f));
      checkCudaErrors(cudaFree(d_g));

      /*run the experiment*/
      printf("Simulation starting:\n");

      run_experiment (d_x,n,d_filter_time,d_filter_freq,d_filter_size,B,B_thresh,
		      loops_loc,loops_thresh,loops_est,W_RLFC,RLFC_loops,
		      repetitions,LARGE_FREQ,k,x_f,t_H_to_D);

      /*cleanup memory*/
      checkCudaErrors(cudaFree(d_filter_time));
      checkCudaErrors(cudaFree(d_filter_freq));


      printf("cuFFT(NVIDIA-CUDA) library results\n");
      printf("******************************************************************************\n");

      int forward = 1;
      int print_y = 1;
      cudafft_results(d_x,n,forward,print_y);

      printf("******************************************************************************\n\n");

}//end make_multiple_t
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*shift r elements of d_x */
void cuShift(complex_t *d_x, int n, int r){
  r = (n + r)%n;
  assert(n >= r);

  complex_t *d_tmp1;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_tmp1),r*sizeof(*d_tmp1)));
  complex_t *d_tmp2;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_tmp2),(n-r)*sizeof(*d_tmp1)));

  checkCudaErrors(cudaMemcpy(d_tmp1,d_x+n-r,r*sizeof(*d_tmp1),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaMemcpy(d_tmp2,d_x,(n-r)*sizeof(*d_tmp2),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaMemcpy(d_x,d_tmp1,r*sizeof(*d_tmp1),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaMemcpy(d_x+r,d_tmp2,(n-r)*sizeof(*d_tmp2),cudaMemcpyDeviceToDevice));
  checkCudaErrors(cudaFree(d_tmp1));
  checkCudaErrors(cudaFree(d_tmp2));
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a device vector d_x */
__global__
void printComDxExKernel(complex_t *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    if (isnan(d_x[i].x)){printf("cuExperiment:d_x[%d].x = %lf \n", i, d_x[i].x);}
  }
}//end printComDxExKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*compute d_filter_time*/
__global__
void makeDolphchebyshevKernel1(complex_t *d_filter_time,int w,double t0,double tolerance){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < w){
    if (isnan(cuCheb(w-1,t0*cos(M_PI*i/(w)))*tolerance)){
      d_filter_time[i].x = 0.0;
    }
    else{
      d_filter_time[i].x = cuCheb(w-1,t0*cos(M_PI*i/(w)))*tolerance;
    }
    d_filter_time[i].y = 0.0;
  }
}// end makeDolphchebyshevKernel1
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*set to zero the imaginary components of  d_filter_time*/
__global__
void makeDolphchebyshevKernel2(complex_t *d_filter_time,int w){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < w){
    d_filter_time[i].y = 0.0;
  }
}// end makeDolphchebyshevKernel2
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*normalize d_filter_freq */
__global__
void makeMultipleKernel1(complex_t *d_filter_freq,double *d_max,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    d_filter_freq[i].x /= *d_max;
    d_filter_freq[i].y /= *d_max;
  }
}// end makeMultipleKernel1
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*normalize d_filter_time */
__global__
void makeMultipleKernel2(complex_t *d_filter_time,int w,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < w){
    d_filter_time[i].x /= n;
    d_filter_time[i].y /= n;
  }
}// end makeMultipleKernel2
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
__device__
double cuCheb(double m, double x){
  if (fabs(x) <= 1){
    return cos(m * acos(x));
  }
  else{
    return cosh(m * acosh(x));
  }
}//end cuCheb
