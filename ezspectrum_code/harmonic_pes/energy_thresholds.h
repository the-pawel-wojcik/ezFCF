#ifndef _energy_thresholds_h_
#define _energy_thresholds_h_

#include "aik_xml_parser.h"
#include "constants.h"
#include "genincludes.h"
#include "genutil.h"

/* Class for handlig information given in the input "energy_thresholds" nodes.*/
class EnergyThresholds {
private:
  bool energy_thresholds_detected;
  double thresh_initial_eV;
  double thresh_target_eV;

  void read_energy_tresholds(const xml_node &node_appox_params);
  double get_energy_thresh(const xml_node &node_energy_thresh,
                           const std::string label);

public:
  EnergyThresholds(const xml_node &node_parent);

  /* Returns true if the energy thresholds were present in the input. */
  bool present() const { return energy_thresholds_detected; }

  double initial_eV() const { return thresh_initial_eV; }
  double target_eV() const { return thresh_target_eV; }
};

#endif
