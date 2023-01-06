#include "vib_state_parser.h"

//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string &ex_str, int &qnt, int &nm) {
  if (ex_str.find("v") == std::string::npos) {
    std::string msg("\nFormat error in [");
    msg += ex_str;
    msg += std::string("] excitation: should contain symbol \'v\'.");
    error(msg);
  }

  ex_str.replace(ex_str.find("v"), 1, " ");
  std::istringstream ex_strs(ex_str);
  ex_strs >> qnt >> nm;

  if (ex_strs.fail()) {
    ex_str.replace(ex_str.find(" "), 1, "v");
    std::string msg("Format error in [");
    msg += ex_str;
    msg += std::string("] excitation. Expected two integers separated by 'v'.");
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
    std::string msg("Empty string represeting a vibrational state.");
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
    // no_of_quanta = 10 
    // mode_number = 3
    int no_of_quanta = -1;
    int mode_number = -1;
    get_qnt_nm(excitation, no_of_quanta, mode_number);

    if (mode_number >= max_nm || mode_number < 0) {
      std::string msg("the_only_initial_state: "
                      " Normal mode out of the range 0 to ");
      msg += std::to_string(max_nm - 1);
      msg += std::string(".");
      error(msg);
    }

    if (no_of_quanta < 0) {
      std::string msg(
          "the_only_initial_state: "
          " Negative number of vibrational quanta in a vibrational mode.");
      error(msg);
    }

    if (v_state[mode_number] != 0) {
      std::string msg("the_only_initial_state: "
                      " The number of vibrational quanta in mode #");
      msg += std::to_string(mode_number);
      msg += std::string(" is defined more than once.");
      error(msg);
    }

    v_state[mode_number] = no_of_quanta;
  }
}
