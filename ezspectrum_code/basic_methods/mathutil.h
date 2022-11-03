#ifndef _mathutil_h_
#define _mathutil_h_

/*! \file mathutil.h
\brief General math utilities.
*/

#include "genincludes.h"

//! Factorial; N=(-inf..0..+inf)  
double Factorial(const int N); 

//! Combination, C_n^k
unsigned long nChoosek(int n, int k);

//! provide next combination once previous is given in j (total number is C_n^k) -- see comments in  .C
bool enumerateCombinations(int n, int k, std::vector <int>& j);

#endif


