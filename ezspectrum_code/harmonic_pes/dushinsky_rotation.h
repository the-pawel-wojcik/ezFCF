#ifndef _dushinsky_rotation_h_
#define _dushinsky_rotation_h_

#include "genincludes.h"
#include "genutil.h"
#include "job_parameters.h"

/* Class for the "dushinsky_rotations" input node. This class is used only for
 * processing the node input and offering access to collected data. This class
 * is not used for computing FCFs.*/
class DushinskyRotation {
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
  DushinskyRotation(const xml_node &node_input, int n_el_states,
                    const JobParameters &jp)
      : targN(1), max_quanta_ini(3), max_quanta_targ(3), Kp_max_to_save(32000),
        n_el_states(n_el_states) {
    std::cout
        << "\n\n=== Reading the Dushinsky rotations job parameters ===\n\n"
        << std::flush;
    // TODO: check that there is just one as in JobParameters:
    // more than one make sense as this would allow to find spectrum
    // for transitions to more than just one target state
    // although it sounds smarter to allow more target states as it
    // is allowed in the parallel approximation
    xml_node node_dushinsky_rotations(node_input, "dushinsky_rotations", 0);
    targN = node_dushinsky_rotations.read_int_value("target_state");
    max_quanta_ini = node_dushinsky_rotations.read_int_value(
        "max_vibr_excitations_in_initial_el_state");
    max_quanta_targ = node_dushinsky_rotations.read_int_value(
        "max_vibr_excitations_in_target_el_state");
    // TODO: test that the input numbers make sense
    if (jp.get_temp() == 0) {
      max_quanta_ini = 0;
      std::cout << "\nSince temperature=0, "
                   "\"max_vibr_excitations_in_initial_el_state\" has been "
                   "set to 0.\n"
                << std::flush;
    }

    // TODO: find and document what is Kp_max_to_save
    if (node_dushinsky_rotations.find_subnode("max_vibr_to_store")) {

      xml_node node_max_vibr_to_store(node_dushinsky_rotations,
                                      "max_vibr_to_store", 0);
      Kp_max_to_save = node_max_vibr_to_store.read_int_value("target_el_state");
    }

    test_targN();
    std::cout << "Duschinsky calculations will work with the input target "
                 "state number "
              << targN << ".\n\n";
  }
  int get_targN() const { return targN; }
  int get_max_quanta_init() const { return max_quanta_ini; }
  int get_max_quanta_targ() const { return max_quanta_targ; }
  int get_Kp_max_to_save() const { return Kp_max_to_save; }
};

#endif
