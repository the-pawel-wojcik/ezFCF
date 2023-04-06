#ifndef _single_excitations_h_
#define _single_excitations_h_

#include "aik_xml_parser.h"
#include "molstate.h"
#include "spectralpoint.h"

/* Container for the single excitations.
 * Parses the "single_excitation" nodes.
 */
class SingleExcitations {
public:
  std::vector<SpectralPoint> single_excitations;

  SingleExcitations(xml_node &node_head, const MolState &targ_elst,
                    const int n_molecular_nms, const int iniN);

  bool empty() const { return single_excitations.empty(); }
};

#endif
