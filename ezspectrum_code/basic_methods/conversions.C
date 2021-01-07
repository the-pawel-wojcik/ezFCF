#include "constants.h"

bool covert_energy_to_eV(double& ene, const std::string& units) {

  bool all_is_good=true;

  //If units were specified .....
  if (units!= "") {   
    if (units=="cm-1")
      ene*=WAVENUMBERS2EV;
    else if (units=="K")
      ene*=KELVINS2EV;
    else if (units!="eV") 
      all_is_good=false;
  }
  return all_is_good;
}

