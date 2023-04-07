#include "genutil.h"
#include "vib_state_parser.h"
#include "vibronic_state.h"
#include <algorithm>

// copy of the void vib_state_to_vec_int(std::string &text, std::vector<int>
// &v_state) function
//
// TODO: Turn vib_state_to_vec_int into a funciton that processes text
// and turns it into a VibronicState
//
// TODO: Make VibronicState return the state as a vector of ints.
VibronicState::VibronicState(std::string &input,
                             const int n_molecular_normal_modes,
                             const int el_st_idx)
    : elStateIndex(el_st_idx), vibrQuanta(n_molecular_normal_modes, 0) {
  // TODO: once transition to pair happens this should not be needed
  // Make sure that the case when no excitates are requested: i.e., "0"
  // is also covered. It needs to be convered as the no excitations vibronic
  // state is almost always active. At the same time it should be implemented
  // so that a serires of ifs for this special case is not needed.
  for (int normal_mode = 0; normal_mode < n_molecular_normal_modes;
       ++normal_mode) {
    excite_subspace.push_back(normal_mode);
  }

  if (input.empty()) {
    std::string msg("Vibrational state represeted by an empty string.");
    error(msg);
  }

  // Special case that differes from the standard format
  // TODO: possible remove leading and trailsing whitespaces
  if (input == "0") {
    return;
  }

  // slice the input string into individual excitations
  std::queue<std::string> requested_nms_with_excitations;
  std::istringstream iss(input);
  std::string single_mode_excitation;
  while (!iss.eof()) {
    iss >> single_mode_excitation;
    requested_nms_with_excitations.push(single_mode_excitation);
  }

  // Extract data from each excitation
  while (requested_nms_with_excitations.size()) {
    std::string excitation = requested_nms_with_excitations.front();
    requested_nms_with_excitations.pop();
    // parse "10v3" into
    // how_excited = 10
    // mode_number = 3
    int how_exicted = -1;
    int mode_number = -1;
    get_qnt_nm(excitation, how_exicted, mode_number);

    if (mode_number >= n_molecular_normal_modes || mode_number < 0) {
      std::stringstream msg;
      msg << " Normal mode (" << mode_number << ") out of bounds." << std::endl;
      msg << "Should be no less than 0 and no larger than "
          << n_molecular_normal_modes - 1 << ".";
      error(msg);
    }

    if (how_exicted < 0) {
      std::stringstream msg;
      msg << "Negative number of vibrational quanta (" << how_exicted
          << ") requested for mode #" << mode_number << "." << std::endl;
      error(msg);
    }

    // TODO: this can only be checked once the objects of this class are
    // not intialized by default to the full space.
    /* auto pos = */
    /*     std::find(excite_subspace.begin(), excite_subspace.end(),
     * mode_number); */
    /* if (pos != excite_subspace.end()) { */
    /*   std::stringstream msg; */
    /*   msg << " Repeated specification of the number of vibrational quanta for
     * " */
    /*          "mode #" */
    /*       << mode_number << "."; */
    /*   error(msg); */
    /* } */

    setVibrQuanta(mode_number, how_exicted);
  }
};

int VibronicState::getTotalQuantaCount() {
  int count = 0;
  for (int nm = 0; nm < getV().size(); nm++)
    count += getV()[nm];
  return count;
}

int VibronicState::getV_full_dim(const int nm) {
  int return_quanta = 0;
  // check if nm in the "excite space"
  for (int i = 0; i < getEx().size(); i++)
    if (nm == getEx()[i])
      return_quanta = getV()[i];

  return return_quanta;
}

std::ostream &operator<<(std::ostream &os, const VibronicState &vibst) {
  int quanta_printed = 0;
  os << vibst.elStateIndex << '(';

  for (int nm = 0; nm < vibst.excite_subspace.size(); nm++) {
    if (vibst.vibrQuanta[nm] > 0) {
      if (quanta_printed > 0) // not the first normal mode in the vector
        os << ',';
      quanta_printed++;
      os << vibst.vibrQuanta[nm] << 'v' << vibst.excite_subspace[nm];
    }
  }

  // no excitations stands for the ground vibrational state, i.e., "(0)"
  if (quanta_printed == 0)
    os << '0';
  os << ")";
  return os;
}
