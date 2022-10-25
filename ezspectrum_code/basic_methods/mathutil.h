#ifndef _mathutil_h_
#define _mathutil_h_

/*! \file mathutil.h
\brief General math utilities: fact, MIN, MAX, SIGN, RANDOM etc.
*/

#include "genincludes.h"

//! Factorial; N=(-inf..0..+inf)  
double Factorial(const int N); 

//! Combination, C_n^k
unsigned long nChoosek(int n, int k);

//! provide next combination once previous is given in j (total number is C_n^k) -- see comments in  .C
bool enumerateCombinations(int n, int k, std::vector <int>& j);

//Max and Min f-s
//! A macro that returns the maximum of a and b
#define MAX(a,b) ( (a>=b) ? (a) : (b) )
//! A macro that returns the minimum of a and b
#define MIN(a,b) ( (a<b) ? (a) : (b) ) 
//! A macro that returns the sign of as a (double)+/-1.
#define SIGN(a)  ( (a<0.) ? (-1.) : (1.) )
//! random number in [0,1] 

#define RANDOM ((double)random()/RANDMAX) 
//! Randomizer (use system time)
void Randomize(void );


//! Making sure bool, TRUE and FALSE are defined
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (1)
#endif


#endif


