/*                 \
 *  Copyright (c) 2012-2013 Haitham Hassanieh, Piotr Indyk, Dina Katabi,
 *  Eric Price, Massachusetts Institute of Technology.
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
 * SOFTWARE.
 * Please refer to the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
	       
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

//includes CUDA project
#include <cuda.h>
#include <cufft.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

extern "C"{
       #include "utils.h"
       #include "timer.h"
       #include "outerLoop.h"
       #include "innerLoop.h"
}

bool_t ALGORITHM1   = true;
bool_t VERBOSE      = false;
bool_t TIMING       = true;

#define vprintf(...) if (VERBOSE){ printf(__VA_ARGS__);}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/* prototype functions */
__global__
void printDxOutKernel(int *d_x,int n);

__global__
void printComDxOutKernel(complex_t *d_x,int n);
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Node *outer_loop(complex_t *d_x,int n,complex_t *d_filter_time,
		 complex_t *d_filter_freq,int d_filter_size,int B,
		 int B_thresh,int W_RLFC,int RLFC_loops,int loop_threshold,
		 int loops_loc,int loops,int *I_F,float t_H_to_D){

  /*Variables used for timing*/
  double VOTE_T = 0;        // vote vector timing
  double PF_T    = 0;       // permutation and filtering partial timing
  double G_T     = 0;       // reverse hash function partial timing
  double BC_T    = 0;       // FFT and cutoff partial timing
  double PF_ALL  = 0;       // permutation and filtering total  timing
  double G_ALL   = 0;       // reverse hash function  total timing
  double BC_ALL  = 0;       // FFT and cutoff total timing
  double RLFC_Time = 0.0;   // restrict location coefficients timing

  double DDD = get_time();

  /*Allocate device  memory for the vector d_bins_f,size = B*loops */
  complex_t *d_bins_f;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_bins_f), sizeof(*d_bins_f)*B*loops));

  /*Allocate device  memory  for vector d_vote */
  double t_initial = get_time();
  int *d_vote;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_vote), sizeof(*d_vote)*n));

  /*Set to zero the components of  d_vote*/
  t_initial = get_time();
  checkCudaErrors(cudaMemset(d_vote, 0, sizeof(*d_vote)*n));
  VOTE_T = get_time() - DDD;

  /*Allocate device  memory  for vector d_I to store the real indices of the largest frequency coefficients */
  int *d_I;
  checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_I), sizeof(*d_I)*MAXSIZE_HITS));

  DDD = get_time();
  /*Begin function for the restriction of indices of largest frequency coefficients
    (SFFT version 2.0)*/
  /*Allocate unified   memory  for vector d_J2 to store the indices of largest frequency coefficients */
  int *d_J2;
  checkCudaErrors(cudaMallocManaged(reinterpret_cast<void **>(&d_J2), sizeof(*d_J2)*RLFC_loops*B_thresh));
  int num_RLFC = B_thresh;
  for(int i = 0; i < RLFC_loops; i++){
    restr_loc_freq_coeff(d_x, n, B_thresh, W_RLFC, d_J2+i*B_thresh, &RLFC_Time);
  }

  if(RLFC_loops > 1){
    radix_sort_t(d_J2, RLFC_loops * B_thresh);
    int Last =0;
    for(int i = 1; i < RLFC_loops * B_thresh; i++){
      if(d_J2[i]!=d_J2[Last]){
	d_J2[++Last]=d_J2[i];
      }
    }
    num_RLFC=Last+1; // number of elements with non repeated  values  in d_J2.
    vprintf("RLFC:%d----->%d\n\n", B_thresh*RLFC_loops,num_RLFC);
  }
  double RLFC_time = get_time()-DDD;
  /*End function for the restriction of  indices of largest frequency coefficients
    (SFFT version 2.0)*/

  /*Begin inner location and estimation Loops*/
  int *permute =(int*)malloc(loops * sizeof(*permute));
  for(int i = 0; i < loops; i++){
    int a = 0;
    while(gcd(a, n) != 1){
      a = (int)(random() % n);
    }
    int ai = mod_inverse(a, n);
    permute[i]=ai;
    int perform_location = (i < loops_loc);
    int *d_J;
    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&d_J), sizeof(*d_J)*B_thresh));
    PF_T = 0.0;
    BC_T = 0.0;
    G_T = 0.0;
    perm_filter_cutoff(d_x,n,d_filter_time,d_filter_size,B_thresh,B,a,ai,d_bins_f+i*B,
		       d_J,&PF_T,&BC_T);
    if (perform_location) {
      reverse_hash(d_J,n, B_thresh,B,a,ai,loop_threshold,d_vote,d_I,I_F,&G_T,
		   d_J2,num_RLFC,W_RLFC,RLFC_loops);
    }
    checkCudaErrors(cudaFree(d_J));
    PF_ALL += PF_T;
    BC_ALL += BC_T;
    if (perform_location){
      G_ALL += G_T;
    }
  }
  checkCudaErrors(cudaFree(d_vote));
  checkCudaErrors(cudaFree(d_J2));
  /*end inner Location Loops*/

  /*Begin estimation*/
  DDD = get_time();

  Node *hat_x = NULL;
  hat_x = estimate_values(d_I,I_F,d_bins_f,loops,n,permute,B,d_filter_freq,loops_loc);

  /*End estimation*/
  double E_T = get_time() - DDD;

  DDD = get_time();

  /*cleanup memory*/
  checkCudaErrors(cudaFree(d_I));
  checkCudaErrors(cudaFree(d_bins_f));
  free(permute);

  if(TIMING){
    printf("Total GPU-SFFT time: %lf\n", DDD);
    printf("Time distribution: votetable  restr-loc-index  perm+filter reverse-hash  estimation FFT+Cutoff    other    total\n");
    printf("                    %2.4lf       %2.4lf          %2.4lf      %2.4lf        %2.4lf     %2.4lf      %2.4lf   %2.4lf\n",
	   VOTE_T, RLFC_time, PF_ALL, G_ALL, E_T, BC_ALL, DDD-PF_ALL-G_ALL-E_T-RLFC_time-BC_ALL-VOTE_T, DDD);

    double tott = (DDD)/100;
    printf("                    %4.1f%%         %4.1f%%           %4.1f%%      %4.1f%%         %4.1f%%      %4.1f%%       %4.1f%%   %5.1f%%\n",
	   VOTE_T / tott, RLFC_time/tott, PF_ALL/tott, G_ALL/tott, E_T/tott, BC_ALL/tott, (DDD-PF_ALL-G_ALL-E_T-RLFC_time-BC_ALL-VOTE_T)/tott, (DDD)/tott);

    printf("\n");
    t_H_to_D = t_H_to_D/1.0e3;
    printf("HtoD time = %lfs, HtoD time/total time  = %lf%%, total time + HtoD time = %lfs \n", t_H_to_D,t_H_to_D/tott,DDD+t_H_to_D);
  }

  return hat_x;
}// end outer_loop

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a int device vector d_x */
__global__
void printDxOutKernel(int *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    printf("outerLoop:PrintDxKernel::d_x[%d] = %d \n", i, d_x[i]);
  }
}//end printDxOutKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*print a complex device vector d_x */
__global__
void printComDxOutKernel(complex_t *d_x,int n){

  int i = threadIdx.x + blockIdx.x * blockDim.x;

  if(i < n){
    printf("innerLoop:PrintComDxKernel::d_x.x[%d] = %lf \n", i, d_x[i].x);
    printf("innerLoop:PrintComDxKernel::d_x.y[%d] = %lf \n", i, d_x[i].y);
  }
}//end printComOutDxKernel
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////