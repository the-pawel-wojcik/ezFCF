#ifndef _genutil_h_
#define _genutil_h_

/*! \file genutil.h
\brief General utilities (get Time; error messages )
\ingroup BASIC_METHODS
*/

#include "genincludes.h"
#include <time.h>

//! Error handling (print msg, exit(1))
void error(const std::string msg);
void error(const char *const msg); // old one

//! if expr -> error(msg);
inline void check(const bool expr, const char *const msg) {
  if (expr)
    error(msg);
}

//! Returns string with the current time:
std::string GetTime();
time_t GetRawTime();

// Whitespace trimming functions from
// https://stackoverflow.com/questions/1798112/removing-leading-and-trailing-spaces-from-a-string
// trim whitespaces from left
inline std::string &ltrim(std::string &s, const char *t = " \t\n\r\f\v") {
  s.erase(0, s.find_first_not_of(t));
  return s;
}

// trim whitespaces from right
inline std::string &rtrim(std::string &s, const char *t = " \t\n\r\f\v") {
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim whitespaces from left & right
inline std::string &trim(std::string &s, const char *t = " \t\n\r\f\v") {
  return ltrim(rtrim(s, t), t);
}

void print_qchem_style_vector(arma::Col<double> v, std::string header);

#endif
