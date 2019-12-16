/*									\
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
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "timer.h"


double cabs_t2(complex_t x){
  return (pow(x.x,2) + pow(x.y,2));
}

double cabs_t(complex_t x){
  return sqrt (pow(x.x,2) + pow(x.y,2));
}

complex_t complexMul(complex_t a, complex_t b) {
  complex_t c;
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.x * b.y + a.y * b.x;
  return c;
}

complex_t complexDiv(complex_t a, complex_t b){

  complex_t c;
  if (b.x == 0.0 && b.y == 0.0){
    printf("utils.c: division by zero is not allowed \n");
    c.x = 1.0e10;
    return c;
  }

  complex_t b_conj;
  b_conj.x = b.x;
  b_conj.y = -b.y;

  c = complexMul(a,b_conj);
  float div = cabs_t2 (b);
  c.x = c.x/div;
  c.y = c.y/div;
  return c;
}

complex_t complexSubstr(complex_t a, complex_t b){
  complex_t c;
  c.x = a.x - b.x;
  c.y = a.y - b.y;
  return c;
}

/* Compute the gcd of a and b
   assumes a, b > 0*/
int gcd(int a, int b){
  if (a%b == 0) return b;
  return gcd(b, a%b);
}

double phase(complex_t x){
  return atan2(x.y, x.x);
}

/*crappy inversion code I stole from elsewhere
  Undefined if gcd(a, n) > 1*/
int mod_inverse(int a, int n) {
  int i = n, v = 0, d = 1;
  while (a>0) {
    int t = i/a, x = a;
    a = i % x;
    i = x;
    x = d;
    d = v - t*x;
    v = x;
  }
  v %= n;
  if (v<0) v = (v+n)%n;
  return v;
}

void radix(int byte, int size, int *A, int *TEMP) {

  int* COUNT = (int*)calloc(256,sizeof(*COUNT));

  byte = byte << 3;

  for (int i = 0; i < size; ++i)
    ++COUNT[((A[i]) >> (byte)) & 0xFF];
  for (int i = 1; i < 256; ++i)
    COUNT[i] += COUNT[i - 1];
  for (int i = size - 1; i >= 0; --i) {
    TEMP[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = A[i];
    --COUNT[(A[i] >> (byte)) & 0xFF];
  }

  free(COUNT);

}

void radix_sort_t(int *A, int size) {
  int* TEMP = (int*)malloc(size*sizeof(*TEMP));

  for (unsigned int i = 0; i < sizeof(int); i += 2) {

    // even byte
    radix(i, size, A, TEMP);

    // odd byte
    radix(i + 1, size, TEMP, A);
  }

  free(TEMP);
}

void radix_filt(int byte, int size, int *A, int *TEMP, complex_t* Filter, complex_t* TMP_F) {

  int* COUNT = (int*)calloc(256,sizeof(*COUNT));

  byte = byte << 3;

  for (int i = 0; i < size; ++i)
    ++COUNT[((A[i]) >> (byte)) & 0xFF];
  for (int i = 1; i < 256; ++i)
    COUNT[i] += COUNT[i - 1];
  for (int i = size - 1; i >= 0; --i) {
    TEMP[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = A[i];
    TMP_F[COUNT[(A[i] >> (byte)) & 0xFF] - 1] = Filter[i];
    --COUNT[(A[i] >> (byte)) & 0xFF];
  }

  free(COUNT);

}

void radix_sort_filt(int *A, complex_t* Filter,int size) {

  int *TEMP = (int*)malloc(size*sizeof(*TEMP));

  complex_t *TMP_F = (complex_t*)malloc(size*sizeof(*TMP_F));

  for (unsigned int i = 0; i < sizeof(int); i += 2) {

    // even byte
    radix_filt(i, size, A, TEMP, Filter, TMP_F);

    // odd byte
    radix_filt(i + 1, size, TEMP, A, TMP_F, Filter);
  }

  free(TEMP);
  free(TMP_F);
}

int floor_to_pow2(double x){
  unsigned int ans;
  for(ans = 1; ans <= x; ans <<= 1)
    ;
  return ans / 2;
}

double AWGN(complex_t *x, int n, double std_noise){

  if(std_noise==0)
    return 1000000000;

  complex_t gn;
  gn.x = 0.0;
  gn.y = 0.0;

  double sig_power =0;
  double noise_power =0;
  double snr;
  double u, v;

  for(int h = 0 ; h < n; h++){
    sig_power += cabs_t(x[h])*cabs_t(x[h]);

    u=drand48();
    v=drand48();
    gn.x = std_noise * sqrt(-2*log(u)) * cos(2*M_PI*v);
    gn.y = std_noise * sqrt(-2*log(u)) * sin(2*M_PI*v);

    noise_power += -2*log(u);

    x[h].x += gn.x;
    x[h].y += gn.y;
  }

  noise_power = noise_power * std_noise * std_noise;
  snr = sig_power/noise_power;

  return snr;
}

double binomial_cdf(double prob, int n, int needed){
  double ans = 0;
  double choose = 1;
  for(int i = n; i >= needed; i--){
    ans += choose * pow(prob, i) * pow(1-prob, n-i);
    choose = choose * i / (n-i+1);
  }
  return ans;
}

int comp_struct(const void *a,const void *b){
  Node *aa=(Node *)a;
  Node *bb=(Node *)b;
  return(((cabs_t(aa->value))<(cabs_t(bb->value)))?1:-1);
}

int comp_struct2(const void *a,const void *b){
  Node *aa=(Node *)a;
  Node *bb=(Node *)b;
  if((aa->key)==(bb->key))
    return 0;

  return((aa->key)>((bb->key))?1:-1);
}

int comp_struct3(const void *a,const void *b){
  Pair *aa=(Pair *)a;
  Pair *bb=(Pair *)b;
  if((aa->first)==(bb->first))
    return 0;

  return((aa->first)>((bb->first))?1:-1);
}
