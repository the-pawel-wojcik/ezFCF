#ifndef _constants_h_
#define _constants_h_

/*! \file constants.h
  \brief Math constants 
  \ingroup METHODS
  */

#include <string>
#include <algorithm>


const double PI = 3.14159265358979323846;

//! Dimension of Cartesian 3D world
const int CARTDIM = 3;

//! c (speed of light) cm/second
const double SPEEDOFLIGHT_CM_SEC =  2.99792458e10;
//! c (speed of light) bohr/a.time.u
const double SPEEDOFLIGHT_AU = 137.0359996;

const double BOHR_2_CM = 0.00000000529177249;
//! h, erg-seconds 
const double PLANKCONSTANT_ERGxSEC = 6.626176e-27;

const double KELVINS2EV = 8.61735e-05;

const double AMU2GRAM = 1.6605655e-24;
const double AU2ANGSTROM =  0.52918;
const double ANGSTROM2AU = 1.889716;
const double WAVENUMBERS2EV = (1.0/8065.479);
const double EV2HARTREE = (1.0/27.211386245988);

const double HBAR = 1.054572669e-34;
const double C_IN_SI =  2.99792458e8;

const double AMU_2_ELECTRONMASS = (1822.8884843);

//! if less than threshold, considered to be zero
const double COORDINATE_THRESHOLD = 0.00000001; //Angstrom
//! if less than threshold, considered to be zero
// threshold for moment of intertia tensor
const double MOI_THRESHOLD = 0.00000001;


//! Convert energy to eV
bool covert_energy_to_eV(double& ene, const std::string& units);

#endif
