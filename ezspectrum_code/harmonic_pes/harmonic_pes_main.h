#ifndef _harmonic_pes_main_
#define _harmonic_pes_main_

/*! \file harmonic_pes__main.h
  \brief  the "main" for analytic-harmonic Franck-Condon Factors calculation
  */

#include "genincludes.h"

#include "molstate.h"
#include "parallel_approximation.h"
#include "dushinsky.h"
#include "vibrational_indexing.h"
#include "aik_xml_parser.h"
#include "vib_state_parser.h"

#include <limits>

//! program itself:
bool harmonic_pes_main(const char *InputFileName, xml_node &node_input,
                       xml_node &node_amu_table);

#endif


