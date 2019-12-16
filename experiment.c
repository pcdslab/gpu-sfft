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
/*
 *
 * Compile with $ make
 * Run with $ ./exec -q -k -Bcst -loops_loc -loops_thresh -loops_est -RLFC_cst -RLFC_loops
 *
 */

#include <getopt.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "cudaFft.h"
#include "filters.h"
#include "outerLoop.h"
#include "timer.h"
#include "utils.h"
#include "experiment.h"

int main(int argc, char **argv){
  int arg = 9;
  if (argc < arg){
      printf("number of arguments less than %d.\n", arg);
      return 0;
  }
  int q = atoi(argv[1]);                         // n = 2^q
  int n = pow(2,q);                              // Input signal size: total number of time samples in the input signal.
  int k = atoi(argv[2]);                         // Input  signal sparsity,i.e., the  number of large frequencies coefficients.
  int repetitions = 1;                           // Number of times the experiment is repeated.
  double Bcst= strtod(argv[3], 0);               // Constant to compute the number of bins to estimate the locations and the values of the largest frequencies coefficients.
  int loops_loc = atoi(argv[4]);                 // Number of loops to estimate the locations of the largest frequencies coefficients.
  int loops_thresh = atoi(argv[5]);              // Number of times a frequency coefficient has to appear in the location bins.
  int loops_est  = atoi(argv[6]);                // Number of loops to estimate values of the  largest frequencies coefficients.
  double RLFC_cst= atoi(argv[7]);                // Constant to compute the number of bins for the function to compute the restriction of largest frequency coefficients.
  int RLFC_loops = atoi(argv[8]);                // Number of loops to run the function to compute the restriction of largest frequency coefficients.
  double snr=1.0e+12;                            // Specifies the signal to noise ratio.
  double tolerance = 1.0e-8;                     // Sets the noise (Leakage) in the filter.
  VERBOSE = false;                               // When true, prints detailed timing of each step of the code.

  double STD_NOISE_Y = 0;
  double std_noise = sqrt(k/(2*snr))*STD_NOISE_Y;

  n = floor_to_pow2(n);

  double BB = (unsigned) (Bcst*sqrt((double)n*k/(log2(n))));
  double lobefrac = 0.5 / (BB);

  int b_f = (int)(1.2*1.1*((double) n/BB));
  assert(b_f <= n);

  int B = floor_to_pow2(BB);
  int B_thresh = 2*k;
  int W_RLFC = floor_to_pow2(RLFC_cst*n/B);
  assert(B_thresh < W_RLFC);
  assert(B_thresh < B);
  assert(loops_thresh <= loops_loc);


  printf("\n\nRUNNING EXPERIMENT: n:%d k:%d loops_loc=%d loops_est:%d\n", n, k,loops_loc,loops_est);

  printf("GPU-SFFT filter parameters for: n=%d, k=%d.\n", n, k);
  printf("******************************************************************************\n");
  printf(" RLFC Filter: RLFC_loops: %d  B_thresh:%d W_RLFC: %d\n",
	 RLFC_loops,B_thresh, W_RLFC);
  printf(" Filter: (numlobes=%.1lf, tolerance=%lg, b_f=%d) B: %d \n",
	 0.5/lobefrac,tolerance,b_f,B);
  printf("\n");

  cudafft_experiment(n,lobefrac,tolerance,b_f,B,B_thresh,loops_loc,
		     loops_thresh,loops_est,W_RLFC,RLFC_loops,
		     repetitions,k,std_noise);

  return 0;
}// end main
