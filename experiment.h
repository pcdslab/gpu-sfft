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

#ifndef EXPERIMENT_H_
#define EXPERIMENT_H_

#include "filters.h"

int cudafft_experiment(int n,double lobefrac,double tolerance,int b_f,
		       int B,int B_thresh,int loops_loc,int loops_thresh,
		       int loops_est,int W_RLFC,int RLFC_loops,
		       int repetitions,int k,double std_noise);

int cudafft_results (complex_t *d_x,int n,int forward,int print_y);

void run_experiment(complex_t *d_x,int n,complex_t *d_filter_time,
		    complex_t *d_filter_freq,int d_filter_size,int B,
		    int B_thresh,int loops_loc,int loops_thresh,
		    int loops_est,int W_RLFC,int RLFC_loops,
		    int repetitions,int *LARGE_FREQ,int k,
		    complex_t *x_f,float t_H_to_D);

void cuShift(complex_t *d_x, int n, int r);

#endif 
