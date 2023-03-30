#include "vib_state_parser.h"
#include <ostream>

//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string &ex_str, int &qnt, int &nm) {
  if (ex_str.find("v") == std::string::npos) {
    std::string msg(" Format error in [");
    msg += ex_str;
    msg += std::string("] excitation: should contain symbol \'v\'.");
    error(msg);
  }

  std::string ex_str_backup(ex_str);

  ex_str.replace(ex_str.find("v"), 1, " ");
  std::istringstream ex_strs(ex_str);
  ex_strs >> qnt >> nm;

  if (ex_strs.fail()) {
    std::stringstream msg;
    msg << "Format error in \"" << ex_str_backup
        << "\".\n Expected two integers separated by 'v'.";
    error(msg);
  }
}

//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state
//! stored as a vector of integers. This is an alternative to the
//! fillVibronicState which stores the vibrational state as
//! an object of VibronicState
//! v_state should be initilized to zeros. Each index of v_state is the number
//! of a normal mode in use
void vib_state_to_vec_int(std::string &text, std::vector<int> &v_state) {

  if (text.empty()) {
    std::string msg("Vibrational state represeted by an empty string.");
    error(msg);
  }

  // number of molecule's normal modes: 3N - 5/6
  const int max_nm = v_state.size();

  std::queue<std::string> non_zero_modes;

  std::string non_zero_mode;
  size_t pos = text.find(",");
  while (true) {
    non_zero_mode = text.substr(0, pos);
    non_zero_modes.push(non_zero_mode);
    if (pos == std::string::npos)
      break;
    text.erase(0, pos + 1);
    pos = text.find(",");
  }

  while (non_zero_modes.size()) {
    std::string excitation = non_zero_modes.front();
    non_zero_modes.pop();
    // parse "10v3" into
    // how_excited = 10
    // mode_number = 3
    int how_exicted = -1;
    int mode_number = -1;
    get_qnt_nm(excitation, how_exicted, mode_number);

    if (mode_number >= max_nm || mode_number < 0) {
      std::stringstream msg;
      msg << " Normal mode (" << mode_number << ") out of bounds." << std::endl;
      msg << "Should be no less than 0 and no larger than " << max_nm - 1
          << ".";
      error(msg);
    }

    if (how_exicted < 0) {
      std::stringstream msg;
      msg << "Negative number of vibrational quanta (" << how_exicted
          << ") requested for mode #" << mode_number << "." << std::endl;
      error(msg);
    }

    if (v_state[mode_number] != 0) {
      std::stringstream msg;
      msg << " Repeated specification of the number of vibrational quanta for "
             "mode #"
          << mode_number << ".";
      error(msg);
    }

    v_state[mode_number] = how_exicted;
  }
}
