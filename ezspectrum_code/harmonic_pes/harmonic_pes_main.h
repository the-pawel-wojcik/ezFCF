#ifndef _harmonic_pes_main_
#define _harmonic_pes_main_

/*! \file harmonic_pes__main.h
\brief  the "main" for analytic-harmonic Franck-Condon Factors calculation
*/


//////////////////////////////////////////////////////////////////////////////
//                                                                          // 
//  calculates photoelectron spectrum in 1D Franck-Condon approximation      // 
//                                                                          // 
//////////////////////////////////////////////////////////////////////////////

//2DO: 
// 1 ADD DEGENERATE NORMAL MODES
// 2 several target states
// cut into nice picies

//! .xml file name with atomName<->atomMass table
#define ATOMIC_MASSES_FILE ("atomicMasses.xml")

#include "genincludes.h"
#include <vector>
#include <queue>
#include <set>

#include "molstate.h"
#include "kmatrix.h"
#include "vector3d.h"

#include "parallel_approximation.h"
#include "dushinsky.h"
#include "vibrational_indexing.h"
#include "aik_xml_parser.h"

#include <limits>

//! program itself: 
bool harmonic_pes_main (const char *InputFileName, xml_node& node_input, xml_node& node_amu_table);


#endif


