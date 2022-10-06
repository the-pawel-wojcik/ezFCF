/*! \file genutil.C
\ingroup BASIC_FUNCTIONS
*/

#include "genutil.h"
#include "genincludes.h"

// Return string with the current time:
std::string GetTime(){
  time_t rawtime;
  time ( &rawtime );
  return  ctime ( &rawtime );
}

time_t GetRawTime(){
  time_t rawtime;
  time ( &rawtime );
  return  rawtime;
}

void error(const char *const msg)
{
  std::cout << "\nError! " << msg <<"\n\n";
  exit(EXIT_FAILURE);
}
