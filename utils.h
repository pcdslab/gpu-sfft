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

#ifndef UTILS_H
#define UTILS_H

#include <math.h>

typedef struct double2C{
  double x;
  double y;
}complex_t;

typedef int bool_t;
#define true 1
#define false 0
#define MAXSIZE_HITS 524288

typedef struct nodes{
  unsigned int key;
  complex_t value;
}Node;


typedef struct pairs{
  int first;
  int second;
}Pair;

/* functions prototypes */
double cabs_t(complex_t x);

double cabs_t2(complex_t x);

complex_t complexMul(complex_t a, complex_t b);

complex_t complexDiv(complex_t a, complex_t b);

complex_t complexSubstr(complex_t a, complex_t b);

int gcd(int a, int b);

double phase(complex_t x);

int mod_inverse(int a, int n);

void radix(int byte, int size, int *A, int *TEMP);

void radix_sort_t(int *A, int size);

void radix_filt(int byte, int size, int *A, int *TEMP, complex_t* Filter, complex_t* TMP_F);

void radix_sort_filt(int *A, complex_t* Filter,int size);

int floor_to_pow2(double x);

double AWGN(complex_t *x, int n, double std_noise);

double binomial_cdf(double prob, int n, int needed);

int comp_struct(const void * a, const void * b);

int comp_struct2(const void * a, const void * b);

int comp_struct3(const void * a, const void * b);

#endif /* UTILS_H_ */
