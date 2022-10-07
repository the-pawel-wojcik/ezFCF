#ifndef _parallel_approximation_h_
#define _parallel_approximation_h_

/*! \file parallel_approximation.h
\ingroup METHODS
\brief calculates and stores Franck-Condon factors in parallel approximation

To initialize provide:
1. vector of molStates; ([0]-initial, [1],[2],... --target states)
2. threshold for FCFs to accept (in_fcf_threshold)
3. target state index =1..Nmax for spectral lines assignment (in_targN)

*/

#include "constants.h"
#include "genincludes.h"
#include "mathutil.h"
#include "vibrational_indexing.h"
#include "vibronic_state.h"

#include "franck_condon.h"

#include "spectralpoint.h"
#include "spectrum.h"

#include "molstate.h"

class Parallel {
  //! spectrum which stores all the points above the intensity threshold
  Spectrum spectrum;

  // index of the initial state in the list of electronic states
  const int iniN;
  // number of atoms in the molecule
  const int n_atoms;
  // number of normal modes in the molecule, i.e., 3*n_atoms-6 or -5 for linear
  const int n_molecule_nm;
  // size of the "excite subspace" = n_molecule_nm - "do_not_excite_subsp.size"
  const int n_active_nm;

public:
  Parallel(std::vector<MolState> &molStates, std::vector<int> &active_nm,
           double fcf_threshold, double temperature, int max_n_initial,
           int max_n_target, bool if_the_only_initial_state,
           std::vector<int> the_only_initial_state, bool if_comb_bands,
           bool if_use_target_nm, bool if_print_fcfs, bool if_web_version,
           const char *nmoverlapFName, double energy_threshold_initial,
           double energy_threshold_target);

  //! returns the up to date spectrum
  Spectrum &getSpectrum() { return spectrum; };
};

#endif
