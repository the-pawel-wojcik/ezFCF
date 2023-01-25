#ifndef _vib_state_parser_h_
#define _vib_state_parser_h_

/*! \file vib_state_parser.h
\brief Handle parsing of the strings describing vibrational state
\ingroup BASIC_METHODS
*/

#include "genincludes.h"
#include "genutil.h"

// TODO: unite it either with VibronicState or with TheOnlyInitialState

//! splits string of type "3v21" into two integers 3 and 21
//! qnt = 3; nm = 21
void get_qnt_nm(std::string& ex_str, int& qnt, int& nm );

//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state 
//! stored as a vector of integers. This is an alternative to the
//! fillVibronicState which stores the vibrational state as 
//! an object of VibronicState
//! v_state should be initilized to zeros. Each index of v_state is the number
//! of a normal mode in use
void vib_state_to_vec_int(std::string& text, std::vector<int>& v_state);

#endif
