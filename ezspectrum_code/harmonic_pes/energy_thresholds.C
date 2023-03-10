#include "energy_thresholds.h"

/* Parses the input node and extracts the energy thresholds for initial and
 * target states.
 * `node_parent` can contain "energy_thresholds" as its subnode. If no
 * "energy_thresholds" subnode is found the `present()` function returns false 
 * and the energy thresholds are set to DBL_MAX.
 */
EnergyThresholds::EnergyThresholds(const xml_node &node_parent)
    : energy_thresholds_detected(false), thresh_initial_eV(DBL_MAX),
      thresh_target_eV(DBL_MAX) {
  if (node_parent.find_subnode("energy_thresholds")) {
    read_energy_tresholds(node_parent);
    energy_thresholds_detected = true;
  }
}

/* A helper function used for parsing of the 'energy_thresholds' subnode
 * from either 'parallel_approximation' or 'dushinsky_rotations' nodes. */
void EnergyThresholds::read_energy_tresholds(
    const xml_node &node_appox_params) {

  std::cout << "Reading energy thresholds.\n\n" << std::flush;
  xml_node node_energy_thresh(node_appox_params, "energy_thresholds", 0);

  if (node_energy_thresh.find_subnode("initial_state")) {
    thresh_initial_eV = get_energy_thresh(node_energy_thresh, "initial_state");
  }

  if (node_energy_thresh.find_subnode("target_state")) {
    thresh_target_eV = get_energy_thresh(node_energy_thresh, "target_state");
  }
}

/* A helper function used in parsing of the "energy_thresholds" node of the
 * input xml file. `label` is one of "initial_state" or or "target_state".
 * Returns the threshold value in eV. */
double EnergyThresholds::get_energy_thresh(const xml_node &node_energy_thresh,
                                           const std::string label) {

  if (label != std::string("initial_state") &&
      label != std::string("target_state")) {
    std::stringstream msg;
    msg << "get_energy_thresh: unknown state label: ";
    msg << label;
    error(msg);
  }

  xml_node node_state_thresh(node_energy_thresh, label, 0);
  std::string units = node_state_thresh.read_string_value("units");
  double energy_thresh = node_state_thresh.read_node_double_value();

  if (!covert_energy_to_eV(energy_thresh, units)) {
    std::stringstream msg;
    msg << "get_energy_thresh: Unknown units of the " << label
        << " energy threshold: \"" << units << "\"\n"
        << "  (accepted units: \"eV\", \"K\", or \"cm-1\").";
    error(msg);
  }

  return energy_thresh;
}
