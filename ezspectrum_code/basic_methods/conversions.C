#include "constants.h"

bool covert_energy_to_eV(double& ene, const std::string& units) {

  //SG: modified to be case insensitive
  bool all_is_good=true;
  std::string tmp_unts = units ;
  std::transform(tmp_unts.begin(), tmp_unts.end(), tmp_unts.begin(), ::tolower) ;

  //If units were specified .....
  if (tmp_unts!= "") {   
    if (tmp_unts=="cm-1")
      ene*=WAVENUMBERS2EV;
    else if (tmp_unts=="k")
      ene*=KELVINS2EV;
    else if (tmp_unts!="ev") 
      all_is_good=false;
  }
  return all_is_good;
}
