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

//includes CUDA project
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

extern "C"{
         #include "cudaFft.h"
         #include "utils.h"
}
/* Compute and returns the FFT of the input signal
 * using cuFFT(NVIDIA-CUDA) library.
 * Prints the timing of the FFT transformation. 
 */ 
int cuda_fft_dft_d(complex_t *h_input,  complex_t *h_output, int nx, int forward, int print_y){

  int batch = 1;

  /*timing variables */
  float time_create_plan;
  float time_transf_H_to_D;
  float time_transf_D_to_H;
  float time_execution;
  float time_total;

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  /*Allocate device  memory  for the input signal*/
  complex_t *d_input;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_input), sizeof(cufftDoubleComplex)*nx*batch));

  /* Copy host memory to device memory*/
  cudaEventRecord(start);
  checkCudaErrors(cudaMemcpy(d_input, h_input, nx*batch*sizeof(complex_t),cudaMemcpyHostToDevice));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time_transf_H_to_D,start, stop);

  /*cuFFT plan creation*/
  cufftHandle plan;
  cudaEventRecord(start);
  checkCudaErrors(cufftPlan1d(&plan,nx,CUFFT_Z2Z,batch));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time_create_plan,start, stop);

  /*cuFFT plan  execution*/
  cudaEventRecord(start);
  checkCudaErrors(cufftExecZ2Z(plan,reinterpret_cast<cufftDoubleComplex *>(d_input),
			       reinterpret_cast<cufftDoubleComplex *>(d_input),forward ? CUFFT_FORWARD:CUFFT_INVERSE));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time_execution,start, stop);

  /*Copy device memory to host memory*/
  cudaEventRecord(start);
  checkCudaErrors(cudaMemcpy(h_output, d_input, nx*batch*sizeof(complex_t),cudaMemcpyDeviceToHost));
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time_transf_D_to_H,start, stop);

  /* Destroy cuFFT context*/
  checkCudaErrors(cufftDestroy(plan));

  /*cleanup memory*/
  checkCudaErrors(cudaFree(d_input));
  checkCudaErrors(cudaEventDestroy(start));
  checkCudaErrors(cudaEventDestroy(stop));

  time_total = time_transf_H_to_D + time_execution + time_transf_D_to_H;
  print_timing (time_create_plan,  time_transf_H_to_D, time_transf_D_to_H,
		time_execution, time_total, print_y);


  return 0;

}// end cuda_fft_dft_d


void print_timing (float time_create_plan, float time_transf_H_to_D,
		   float time_transf_D_to_H, float time_execution,
		   float time_total, int print_y){

  if (print_y) {
    printf("Memory transfer host to device execution time %4.2fmilisec\n",time_transf_H_to_D);
    printf("Plan creation execution time %4.2fms\n", time_create_plan);
    printf("cuFFT execution time    %4.2fms\n", time_execution);
    printf("Memory transfer device to host execution time %4.2fms\n",time_transf_D_to_H);
    printf("total time %4.2fms\n",time_total);
  }

}//end print_timing
