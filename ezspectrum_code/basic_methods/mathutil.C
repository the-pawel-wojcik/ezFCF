/*! \file mathutil.C
\ingroup BASIC_METHODS
*/

#include "mathutil.h"

//! Factorial function
double Factorial(const int N)
{
  double tmp=double(abs(N));
  
  if(tmp==0)
    return 1.0;
  
  for(int j=abs(N)-1;j>1;j--)
    tmp*=j;
  
  return tmp;
}

//! Combination function (both parts should be less than max-unsigned-int)
// nChoosek(n, k) gives the binomial coeff
// $$ \binom{n}{k} = \frac{n!}{k!(n-k)!} $$
// This implementation avoids division of large numbers
// https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
unsigned long nChoosek(int n, int k) {
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;

  unsigned long result = n;
  for (int i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }
  return result;
}

//! provide next combination once previous is given in j (total number is C_n^k)
bool enumerateCombinations(int n, int k, std::vector <int>& j)
/*-------------------------------------------------------------------
NOTE: k>0!!

 Description:
   Enumerates all possible combinations of choosing k objects from n
   distint objects.  Initialize the enumeration by setting j[0] to a
   negative value.  Then, each call to enumerateCombinations will
   generate the next combinations and place it in j[0..k-1].  A
   return value of false indicates there are no more combinations to
   generate.  j needs to be allocated with a size for at least k
   elements.

 Usage:
   int comb[10] = {-1};
   while(enumerateCombinations(10, 5, comb)){
       // do something with comb[0..4]
   }

 Reference:
   Tucker, Allen.  Applied Combinatorics.  3rd Ed.  1994.
-------------------------------------------------------------------*/
{
    int i;
    if(j[0] < 0){
        for(i = 0; i < k; ++i){
            j[i] = i;
        }
        return true;
    } else {
        for(i = k - 1; i >= 0 && j[i] >= n - k + i; --i){}
        if(i >= 0){
            j[i]++;
            int m;
            for(m = i + 1; m < k; ++m){
                j[m] = j[m-1] + 1;
            }
            return true;
        } else {
            return false;
        }
    }
}
