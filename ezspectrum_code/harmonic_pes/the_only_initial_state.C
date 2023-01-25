#include "the_only_initial_state.h"

TheOnlyInitialState::TheOnlyInitialState(const xml_node &node, int n_mol_nms)
    : n_molecular_nms(n_mol_nms), node_present(false),
      the_only_state(n_mol_nms, 0) {

  int no_toiss = node.find_subnode("the_only_initial_state");

  if (no_toiss != 1) {
    // TODO: scream error if more than one node exists
    return;
  }

  node_present = true;

  xml_node node_the_only_initial_state(node, "the_only_initial_state", 0);
  std::string text = node_the_only_initial_state.read_string_value("text");
  vib_state_to_vec_int(text, the_only_state);
}
