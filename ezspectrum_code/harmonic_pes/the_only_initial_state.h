#ifndef _the_only_initial_state_h_
#define _the_only_initial_state_h_

#include "genutil.h"
#include "genincludes.h"
#include "aik_xml_parser.h"
#include "vib_state_parser.h"

/* Class for handing the access to data from the "the_only_initial_state" input
 * node. */
// TODO: add input test 
// TODO: unite it with the Duschinsky version VibronicState class
class TheOnlyInitialState {
private:
  /* Number of molecular normal modes, i.e., 3 N_{atoms} - 5/6. */
  const int n_molecular_nms;
  /* A subnode "the_only_initial_state" exist in the input. */
  bool node_present;

  /* Containers for values of arguments passed to the input node */
  std::vector<int> the_only_state;

public:
  /* ```node``` points to an input node that might contain the
   * "the_only_initial_state" as its subnode. If there is no
   * "do_not_excite_subspace" subnode the object will still work: use the
   * member function ```empty()``` to check the status of the object.
   *
   * ```n_mol_nms``` is the number of molecular normal modes, i.e.,
   * 3 (# atoms) - 5 or 6. */
  TheOnlyInitialState(const xml_node &node, int n_mol_nms);

  std::vector<int> get_state() const { return the_only_state; }
  bool non_empty() const { return node_present; }
};

#endif
