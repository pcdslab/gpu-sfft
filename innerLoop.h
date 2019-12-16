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
#ifndef INNERLOOP_H_
#define INNERLOOP_H_


/* prototype functions */

int restr_loc_freq_coeff(complex_t *d_x,int n,int B_thresh,int W_RLFC,
			 int *d_J2,double *RLFC_Time);

int samples_indices(complex_t *d_samp,int B,int samp_size,int* d_output,
		    double *CTime);

int perm_filter_cutoff(complex_t *d_x,int n,complex_t *d_filter_time,
		       int d_filter_size, int B_thresh,int B,int a,
		       int ai,complex_t *d_bins_f,int *d_J,
		       double *PF_T,double *BC_T);

int reverse_hash(int *d_J,int n,int B_thresh,int B,int a,int ai,int loop_thresh,
		 int *d_vote, int *d_I,int *I_F,double *G_T,
		 int *d_J2,int num_RLFC,int W_RLFC,int RLFC_loops);

Node *estimate_values(int *d_I,int *I_F,complex_t *d_bins_f,int loops,int n,
		      int *permute,int  B,complex_t *d_filter_freq,int loops_loc);

int timesmod(const int x,const int a,const int n);

#endif
