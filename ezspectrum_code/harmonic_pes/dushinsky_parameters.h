#ifndef _dushinsky_rotation_h_
#define _dushinsky_rotation_h_

#include "genincludes.h"
#include "genutil.h"
#include "job_parameters.h"

/* Class for the "dushinsky_rotations" input node. This class is used only for
 * processing the node input and offering access to collected data. This class
 * is not used for computing FCFs.*/
class DushinskyParameters {
private:
  // Dushinsky works with only one target state
  int targN;
  int max_quanta_ini;
  int max_quanta_targ;
  int Kp_max_to_save;

  /* Helpers */
  int n_el_states;
  void test_targN() {
    if ((targN < 1) or (targN > n_el_states - 1)) {
      std::stringstream msg;
      msg << "\"dushinsky_rotations\"->\"target_state\" value = " << targN
          << " is incorrect.\n"
          << " Pick no less than 1 and at most " << n_el_states - 1 << ".";
      error(msg);
    }
  }

public:
  /* ```n_el_states```: number of electronic states in the input, i.e., if the
   * input has only the initial state and one target state ```n_el_states```
   * should be set as equal to 2. */
  DushinskyParameters(const xml_node &node_input, int n_el_states,
                    const JobParameters &jp);
  int get_targN() const { return targN; }
  int get_max_quanta_init() const { return max_quanta_ini; }
  int get_max_quanta_targ() const { return max_quanta_targ; }
  int get_Kp_max_to_save() const { return Kp_max_to_save; }
};

#endif
