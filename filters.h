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

#ifndef FILTERS_H_
#define FILTERS_H_

#include "utils.h"


/**
 * Create a window function such that:
 * the main lobe has width 2 * n * filterfrac
 * outside the main lobe, the linf residual is tolerance
 * Computes the required w.
 * Allocates and computes the time component of the filter.
 * Modifies a w-dimensional window function to have n-dimensional FFT
 * the sum of b adjacent ones previously.
 * Allocates and computes the frequency component of the filter.
 * Start running the experiment.
 */
void make_multiple_t(double lobefrac,double tolerance,int n,int b_f,
		     complex_t *d_x,int B,int B_thresh,int loops_loc,
		     int loops_thresh,int loops_est,int W_Comb,int Comb_loops,
		     int repetitions,int *LARGE_FREQ,int k,complex_t *x_f,
		     float t_H_to_D);

#endif 
