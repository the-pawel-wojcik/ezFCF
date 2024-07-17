#ifndef _genincludes_h_
#define _genincludes_h_

/*! \file genincludes.h
\brief General C/C++/STL includes
\ingroup (BASIC_METHODS)
*/

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>

#include <limits> // used for std::numeric_limits<double>::max()

#include <queue>
#include <set>
#include <string>
#include <vector>

#include <algorithm>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <armadillo>

const std::string ATOMIC_MASSES_FILE("atomicMasses.xml");
const std::string GLOBAL_DATA_PATH("/usr/local/share/ezFCF");
const std::string ENVIRONMENT_VAR_NAME("EZFCF_ROOT");

#endif
