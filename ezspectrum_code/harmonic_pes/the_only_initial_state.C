#include "the_only_initial_state.h"

VibronicState TheOnlyInitialState::get_vibronic_state(int elState_idx) const {

  VibronicState vib_state;
  // initialization
  // TODO: make it a ctor of the VibronicState
  for (int mode_number = 0; mode_number < n_molecular_nms; mode_number++) {
    int no_of_excitations = the_only_state[mode_number];
    vib_state.addVibrQuanta(no_of_excitations, mode_number);
  }

  vib_state.setElStateIndex(elState_idx);

  return vib_state;
}

TheOnlyInitialState::TheOnlyInitialState(const xml_node &node, int n_mol_nms)
    : n_molecular_nms(n_mol_nms), node_present(false),
      the_only_state(n_mol_nms, 0) {

  int no_toiss = node.find_subnode("the_only_initial_state");

  if (no_toiss != 1) {
    // TODO: error if more than one node exists
    // See similar issue in the JobParamerets ctor
    return;
  }

  node_present = true;

  xml_node node_the_only_initial_state(node, "the_only_initial_state", 0);
  std::string text = node_the_only_initial_state.read_string_value("text");
  vib_state_to_vec_int(text, the_only_state);
}
