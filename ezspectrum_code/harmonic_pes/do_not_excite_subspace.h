#ifndef _do_not_excite_subspace_h_
#define _do_not_excite_subspace_h_

#include "genutil.h"
#include "genincludes.h"
#include "aik_xml_parser.h"

/* Class for handing the access to data from the "do_not_excite_subspace" input
 * node. */
class DoNotExcite {
private:
  /* Two containers for arguments passed to the input node */
  int size;
  std::set<int> subspace;

  /* Number of molecular normal modes, i.e., 3 N_{atoms} - 5/6. */
  const int n_molecular_nms;

  /* Helpers */
  void parse_normal_modes_stream(std::stringstream & nmodes_stream);

  /* Functions testing user's input */
  void run_tests() const;
  void check_the_largest(const int n_mol_nms) const;
  void check_the_smallest() const;
  void check_for_duplicates() const;

public:
  DoNotExcite(const xml_node &node, int n_mol_nms);

  /* Get number of modes that will be excluded, i.e., the size of the do not
   * excite subspace. */
  int get_size() const { return size; }

  std::set<int> get_as_int_set() const { return subspace; }
  std::vector<int> get_active_subspace() const;

  /* Check if there are any modes to be excluded in the calculations. */
  bool empty() const { return subspace.empty(); }
  bool non_empty() const { return !empty(); }

  void print_summary(const bool nm_were_reordered) const;
};

#endif
