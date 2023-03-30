#include "harmonic_pes_main.h"
#include "molstate.h"

void harmonic_pes_parallel(xml_node &node_input,
                           std::vector<MolState> &elStates,
                           const std::string InputFileName);
void harmonic_pes_dushinksy(xml_node &node_input,
                            std::vector<MolState> &elStates,
                            const std::string InputFileName);

/* *** Check functions *** */

/* Don't continue without at least one target state. */
void require_target_states(std::vector<MolState> &elStates) {
  if (elStates.size() <= 1) {
    error("\nError! No target states found in the input.\n\n");
  }
}

/* VGA doesn't make sense for the initial states. */
void disallow_initial_state_gradient(std::vector<MolState> &elStates) {
  if (elStates[0].IfGradient()) {
    error("\nError! Use of the vertical gradient method is allowed only "
          "in target states.\n\n");
  }
}

/* Normal modes reordering is allowed in target states. */
void disallow_initial_state_nm_reordering(std::vector<MolState> &elStates) {
  if (elStates[0].ifNMReorderedManually()) {
    std::string msg("Manual reordering of normal modes is not "
                    "allowed in the initial state.");
    error(msg);
  }
}

/*
  Checks for:
  - the same number of atoms,
  - the same order of the atomic names,
  - the same "linearity"
*/
void require_similarity(std::vector<MolState> &elStates, int state_i) {
  if (not(elStates[state_i].ifSimilar(elStates[0]))) {
    std::string msg = std::string("target state #");
    msg += std::to_string(state_i);
    msg += std::string(
        " has different # of atoms, their order or \"linear\" flag.");
    error(msg);
  }
}

/* VGA is allowed in all or none of the target states */
void require_consistent_VGA(std::vector<MolState> &elStates, int state_i) {
  // ^ is a binary XOR: returns True if the two bools are different
  bool gradient_used_in_the_first_target = elStates[1].IfGradient();
  if (elStates[state_i].IfGradient() ^ gradient_used_in_the_first_target) {
    std::string msg = std::string("Error: target state #");
    msg += std::to_string(state_i);
    msg += std::string(
        " breaks the consistent use of the VGA in target states.\n\n"
        " Make sure to use vertical gradient in all target states.");
    error(msg);
  }
}

/* Input states have to satisfy some constraints. */
void run_checks(std::vector<MolState> &elStates) {
  require_target_states(elStates);
  disallow_initial_state_gradient(elStates);
  disallow_initial_state_nm_reordering(elStates);

  // run target state checks
  for (int state_i = 1; state_i < elStates.size(); state_i++) {
    // TODO: make sure that the program aborts if any of geometry, nm, or freqs
    // were not read in the non-VG node
    require_similarity(elStates, state_i);
    require_consistent_VGA(elStates, state_i);
  }
}

/* *** End of check functions *** */

/* Reorient and shift molecular geometries. Print new ones.
 * The print_nms variable prints transformed normal modes. */
void run_transformations(std::vector<MolState> &elStates, bool print_nms) {
  for (int state_i = 0; state_i < elStates.size(); state_i++) {

    // unless manual transformations were requested, run the default ones
    if (not(elStates[state_i].ifAlignedManually())) {
      elStates[state_i].align();
      // align each target state with the initial one
      if (state_i > 0)
        elStates[state_i].align(elStates[0]);
      // get clean zeros
      elStates[state_i].applyCoordinateThreshold(COORDINATE_THRESHOLD);
    }

    std::cout << "\nNew molecular geometry:\n";
    elStates[state_i].printGeometry();
    std::cout << std::fixed << std::setprecision(6) << std::setw(13);
    // Armadillo's `raw_print` respects cout flags
    elStates[state_i].getMomentOfInertiaTensor().raw_print("\nMOI tensor:");

    if (print_nms) {
      std::cout << "Normal modes after the geometry transformations:\n\n";
      elStates[state_i].printNormalModes();
    }
  }

  std::cout << "Done with the transformations" << std::endl;
  std::string line(80, '-');
  std::cout << line << "\n";
}

/* *** Constructor *** */
bool harmonic_pes_main(const std::string InputFileName, xml_node &node_input,
                       xml_node &node_amu_table) {
  //============================================================================
  // Read initial state and N target states; i.e. (N+1) electronic states total
  std::vector<MolState> elStates;

  // TODO: separate Reading from processing
  // TODO: Read what's recomended depending on if VG requested or not
  std::cout << "\n====== Reading the initial state ======\n";
  xml_node node_istate(node_input, "initial_state", 0);
  MolState elSt;
  elSt.Read(node_istate);
  elSt.ApplyKeyWords(node_amu_table, elSt);
  elStates.push_back(elSt);

  size_t n_target_states = node_input.find_subnode("target_state");
  for (int state_i = 0; state_i < n_target_states; state_i++) {
    MolState elSt_t;
    xml_node node_t_state(node_input, "target_state", state_i);
    std::cout << "===== Reading the target state #" << state_i << " =====\n";
    elSt_t.Read(node_t_state);
    elSt_t.ApplyKeyWords(node_amu_table, elStates[0]);
    elStates.push_back(elSt_t);
  }
  std::cout << "===== Done reading states =====\n\n";

  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes = node_input.read_flag_value("print_normal_modes");

  // Perform various checks and transformations
  run_checks(elStates);
  run_transformations(elStates, if_print_normal_modes);

  // if parallel or dushinsky
  bool if_something_to_do = false;

  if (node_input.find_subnode("parallel_approximation")) {
    if_something_to_do = true;
    harmonic_pes_parallel(node_input, elStates, InputFileName);
  }

  if (node_input.find_subnode("dushinsky_rotations")) {
    if_something_to_do = true;
    harmonic_pes_dushinksy(node_input, elStates, InputFileName);
  }

  if (!if_something_to_do) {
    error("Input is missing both \"parallel_approximation\" and "
          "\"dushinsky_rotations\" sections.\n Nothing to do.");
  }

  return true;
}

/* Returns true if normal modes reordering was requested for at least one state.
 */
bool normal_modes_reordered(std::vector<MolState> &elStates) {
  return std::any_of(
      elStates.begin(), elStates.end(),
      [](const MolState &state) { return state.ifNMReorderedManually(); });
}

/* print the overlap matrix with the initial state for each target states
 * TODO: getNormalModeOverlapWithOtherState should be const.
 * TODO: should use const DoNotExcite & insted of the set. */
void print_overlap_matrix(std::vector<MolState> &elStates,
                          std::set<int> do_not_excite_subspace,
                          bool if_print_normal_modes) {
  for (int state = 1; state < elStates.size(); state++) {
    std::cout << "\n===== Overlap matrix of the target state #" << state
              << " with the initial state =====\n";

    std::vector<int> normal_modes_list;
    arma::Mat<double> NMoverlap; // normal modes overlap matrix (for each
                                 // target state the same matrix is used)
    bool if_overlap_diagonal;

    // select nondiagonal submatrix of the overlap matrix:
    if_overlap_diagonal = elStates[state].getNormalModeOverlapWithOtherState(
        elStates[0], NMoverlap, normal_modes_list);
    // rows -- norm modes of the target state; colums norm modes of the
    // initial state;

    // remove normal modes from normal_modes_list that are in the
    // do_not_excite_subspace:
    std::vector<int> new_normal_modes_list;

    for (int nm = 0; nm < normal_modes_list.size(); nm++) {
      // if nm is not in the do_not_excite set:
      auto iter_set = do_not_excite_subspace.find(normal_modes_list[nm]);
      if (iter_set == do_not_excite_subspace.end())
        // then copy it to the new list:
        new_normal_modes_list.push_back(normal_modes_list[nm]);
    }

    // Create an overlap submatrix:
    if ((if_overlap_diagonal) or (new_normal_modes_list.size() <= 1)) {
      std::cout << "The normal modes overlap matrix with the initial state "
                   "is diagonal.\n";
      // TODO: This message can be interpreted both ways. Make it clearer
      // and then add it back.
      /* if (!do_not_excite_subspace.empty()) { */
      /*   std::cout */
      /*       << "  (do_not_excite_subspace is excluded)\n"; */
      /* } */
      std::cout << "\n";
    } else {
      std::cout << "WARNING! The normal modes overlap matrix with the "
                   "initial state\n"
                << "         is non-diagonal! Consider reordering the normal "
                   "modes.\n\n";
      // create a normal mode submatrix:
      arma::Mat<double> overlap_submatrix(new_normal_modes_list.size(),
                                          new_normal_modes_list.size(),
                                          arma::fill::zeros);
      for (int nm1 = 0; nm1 < new_normal_modes_list.size(); nm1++)
        for (int nm2 = 0; nm2 < new_normal_modes_list.size(); nm2++)
          overlap_submatrix(nm1, nm2) =
              NMoverlap(new_normal_modes_list[nm1], new_normal_modes_list[nm2]);

      // print the overlap_submatrix (with correct column/row labels):
      std::cout << "  The non-diagonal part of the normal modes overlap matrix "
                   "(do_not_excite_subspace is excluded):";
      std::cout << "\n     ";

      // print header row: column numbers
      for (int j = 0; j < new_normal_modes_list.size(); j++)
        std::cout << std::fixed << std::setprecision(0) << std::setw(8)
                  << new_normal_modes_list[j];

      // print rows
      for (int i = 0; i < new_normal_modes_list.size(); i++) {
        // print row number
        std::cout << "\n  " << std::fixed << std::setprecision(0)
                  << std::setw(3) << new_normal_modes_list[i];

        // print overlap values
        for (int j = 0; j < new_normal_modes_list.size(); j++)
          if (fabs(overlap_submatrix(i, j)) >= 0.001)
            std::cout << std::fixed << std::setprecision(3) << std::setw(8)
                      << overlap_submatrix(i, j);
          else
            std::cout << "      --";
      }
      std::cout << "\n\n";
    }

    // print in a "fit 80 chars wide terminal" form
    if (if_print_normal_modes)
      NMoverlap.print("Normal modes overlap matrix with the initial state "
                      "\n(if significantly non diagonal, please consider "
                      "normal modes reordering)");
  }
}

void print_warnings(bool ifAnyNormalModesReordered,
                    DoNotExcite no_excite_subspace) {

  if (ifAnyNormalModesReordered)
    std::cout
        << "\n"
        << "WARNING! The normal modes of one of the target states were "
           "reordered!\n"
        << "         New order is used for the target state assignment.\n";

  if (no_excite_subspace.non_empty()) {

    std::cout << "\nNOTE: only the following normal modes were excited: "
                 "(\"excite subspace\"):\n  ";

    for (int nm : no_excite_subspace.get_active_subspace())
      std::cout << nm << ' ';
    std::cout << "\n";

    if (ifAnyNormalModesReordered)
      std::cout << "\nWARNING! The normal modes of one of the target states "
                   "were reordered!\n"
                << "         New order is used for the \"excite subspace\"\n";
  }
  std::cout << "\n";
}

//======================================================================
// Parallel approximation
//======================================================================
void harmonic_pes_parallel(xml_node &node_input,
                           std::vector<MolState> &elStates,
                           const std::string InputFileName) {

  xml_node node_jobparams(node_input, "job_parameters", 0);

  // read global paramters
  double temperature = node_jobparams.read_double_value("temperature");
  // fcf threshold (from the <job_parameters> tag)
  double fcf_threshold =
      sqrt(node_jobparams.read_double_value("spectrum_intensity_threshold"));
  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes = node_input.read_flag_value("print_normal_modes");
  // check if the web version format of the output (do not print the input
  // file & create a ".nmoverlap" file)
  bool if_web_version = node_input.read_flag_value("if_web_version");

  // TODO: this should be a call to the pes_main once it turns into a class
  bool ifAnyNormalModesReordered = normal_modes_reordered(elStates);

  // total number of the normal modes (in the initial state)
  int n_molecular_nms = elStates[0].NNormModes();

  std::cout << "\n=== Reading the parallel approximation job parameters ===\n"
            << std::flush;

  xml_node node_parallel_approx(node_input, "parallel_approximation", 0);

  // maximum number of vibrational levels for initial and target state
  // i.e. 2 means three vibr. states: ground and two excited states
  int max_n_initial = node_parallel_approx.read_int_value(
      "max_vibr_excitations_in_initial_el_state");
  int max_n_target = node_parallel_approx.read_int_value(
      "max_vibr_excitations_in_target_el_state");

  // states with excitations in more than one mode
  bool if_comb_bands =
      node_parallel_approx.read_bool_value("combination_bands");

  // Use normal modes of the target state
  bool if_use_target_nm = node_parallel_approx.read_bool_value(
      "use_normal_coordinates_of_target_states");

  // TODO: move it to processing
  if (temperature == 0) {
    max_n_initial = 0;
    std::cout
        << "\nSince temperature=0, "
           "\"max_vibr_excitations_in_initial_el_state\" has been set to 0.\n"
        << std::flush;
  }

  // TODO: Add a use example to Samples/
  bool if_print_fcfs =
      node_parallel_approx.read_flag_value("print_franck_condon_matrices");

  EnergyThresholds thresholds(node_parallel_approx);

  // read normal modes do_not_excite subspace
  DoNotExcite no_excite_subspace(node_parallel_approx, n_molecular_nms);
  no_excite_subspace.print_summary(ifAnyNormalModesReordered);

  std::set<int> do_not_excite_subspace = no_excite_subspace.get_subspace();

  // create active_nm -- "excite subspace" (full_space-do_not_excite_subspace)
  std::vector<int> active_nm_parallel =
      no_excite_subspace.get_active_subspace();

  print_overlap_matrix(elStates, do_not_excite_subspace, if_print_normal_modes);

  // for the web version: save the overlap matrix (with displacements) in an
  // xml file
  std::string nmoverlapFName = InputFileName + std::string(".nmoverlap");

  std::cout << HorizontalLine << "\n"
            << " Beginning the parallel mode approximation computations.\n"
            << HorizontalLine << "\n\n"
            << std::flush;

  // Read "the_only_initial_state" node from the input
  TheOnlyInitialState initial_vibrational_state(node_parallel_approx,
                                                n_molecular_nms);

  Parallel parallel(elStates, active_nm_parallel, fcf_threshold, temperature,
                    max_n_initial, max_n_target, initial_vibrational_state,
                    if_comb_bands, if_use_target_nm, if_print_fcfs,
                    if_web_version, nmoverlapFName, thresholds.initial_eV(),
                    thresholds.target_eV());

  //--------------------------------------------------------------------------------
  // Print the updated spectrum:
  parallel.getSpectrum().Sort();
  std::cout << HorizontalLine << "\n";
  std::cout
      << "           Stick photoelectron spectrum (parallel approximation)\n";
  std::cout << HorizontalLine << "\n";

  print_warnings(ifAnyNormalModesReordered, no_excite_subspace);

  parallel.getSpectrum().PrintStickTable();

  // save this spectrum to the file
  std::string spectrumFName = InputFileName + std::string(".spectrum_parallel");
  parallel.getSpectrum().PrintStickTable(spectrumFName);
  std::cout << "\nStick spectrum was also saved in \"" << spectrumFName
            << "\" file \n";
  if (no_excite_subspace.non_empty())
    std::cout << " (Full list of the normal modes was used for assigning "
                 "transitions)\n";

  std::cout << HorizontalLine << "\n\n";
}

/*
 * =============================================================================
 *  Dushinski rotation (reach exact solution within harmonic approximation)
 * =============================================================================
 *  Notation and equations are from [Berger et al. JPCA 102:7157(1998)]
 * =============================================================================
 */
void harmonic_pes_dushinksy(xml_node &node_input,
                            std::vector<MolState> &elStates,
                            const std::string InputFileName) {

  JobParameters job_config(node_input);
  DushinskyParameters dush_config(node_input, elStates.size(), job_config);

  const int iniN = 0;
  const int targN = dush_config.get_targN();

  // total number of the normal modes (in the molecule)
  // TODO: this is likely the most used variable throughout the program it
  // should be treated in a more general, unified way
  const int n_molecular_normal_modes = elStates[0].NNormModes();

  xml_node node_dushinsky_rotations(node_input, "dushinsky_rotations", 0);
  DoNotExcite no_excite_subspace(node_dushinsky_rotations,
                                 n_molecular_normal_modes);
  no_excite_subspace.print_summary(elStates[targN].ifNMReorderedManually());
  EnergyThresholds thresholds(node_dushinsky_rotations);
  TheOnlyInitialState the_only_init_state(node_dushinsky_rotations,
                                          n_molecular_normal_modes);

  std::cout << HorizontalLine << "\n"
            << " Beginning computations with an inclusion of the Duschinsky "
               "effect.\n"
            << HorizontalLine << "\n\n"
            << std::flush;

  // TODO: move iniN to arg of dushinsky and remove targN from the list as it
  // is already in the dush_config
  Dushinsky dushinsky(elStates, targN, thresholds, dush_config, job_config,
                      no_excite_subspace, the_only_init_state);

  //----------------------------------------------------------------------
  // now load the list of single transitions to evaluate FCFs recursively
  // single_transition is in the "full space"; do_not_excite_subspace does not
  // apply;

  size_t n_single_ex =
      node_dushinsky_rotations.find_subnode("single_excitation");
  if (n_single_ex > 0) {
    if (elStates[targN].ifNMReorderedManually()) {
      std::cout
          << "WARNING! The normal modes of the target state were reordered!\n"
          << "         New order is used for the single transitions.\n\n";
    }

    std::cout
        << "The following single transitions were added to the spectrum:\n"
        << std::flush;

    for (size_t nsex = 0; nsex < n_single_ex; nsex++) {

      xml_node node_single_ex(node_dushinsky_rotations, "single_excitation",
                              nsex);

      std::string initial_state_str = node_single_ex.read_string_value("ini");
      VibronicState init_vibronic_st(initial_state_str,
                                     n_molecular_normal_modes, iniN);

      std::string target_state_str = node_single_ex.read_string_value("targ");
      VibronicState targ_vibronic_st(target_state_str, n_molecular_normal_modes,
                                     iniN);

      SpectralPoint single_transition(init_vibronic_st, targ_vibronic_st);

      // evaluate FCF for each transition and add to the spectrum:
      int K = single_transition.getVibrState1().getTotalQuantaCount();
      int Kp = single_transition.getVibrState2().getTotalQuantaCount();
      // std::cout << "K=" << K << " Kp=" << Kp << std::endl;

      double s_fcf = dushinsky.evalSingleFCF_full_space(
          single_transition.getVibrState1(), K,
          single_transition.getVibrState2(), Kp);
      dushinsky.addSpectralPoint(s_fcf, single_transition.getVibrState1(),
                                 single_transition.getVibrState2());

      std::cout << "FCF=" << std::scientific << std::setprecision(6) << s_fcf
                << " ";
      single_transition.getVibrState1().print();
      std::cout << "->";
      single_transition.getVibrState2().print();
      std::cout << "\n" << std::flush;
    }
  }

  std::cout
      << "\nUpdating the energies and applying the Boltzmann distribution..."
      << std::flush;
  //--------------------------------------------------------------------------------
  // update (fill) energies for every point in the spectrum and add the
  // Boltzmann distribution:
  int points_removed = 0;
  double temperature = job_config.get_temp();
  for (int pt = 0; pt < dushinsky.getSpectrum().getNSpectralPoints(); pt++) {
    double energy = -elStates[targN].Energy();
    double E_prime_prime = 0; // no hot bands

    // run it over the full space, if nm not in the nms_dushinsky subspace,
    // getV_full_dim() returns zero (no excitations):
    for (int nm = 0; nm < elStates[iniN].NNormModes(); nm++) {
      energy += elStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV *
                dushinsky.getSpectrum()
                    .getSpectralPoint(pt)
                    .getVibrState1()
                    .getV_full_dim(nm);
      energy -= elStates[targN].getNormMode(nm).getFreq() * WAVENUMBERS2EV *
                dushinsky.getSpectrum()
                    .getSpectralPoint(pt)
                    .getVibrState2()
                    .getV_full_dim(nm);

      E_prime_prime += elStates[iniN].getNormMode(nm).getFreq() *
                       WAVENUMBERS2EV *
                       dushinsky.getSpectrum()
                           .getSpectralPoint(pt)
                           .getVibrState1()
                           .getV_full_dim(nm);
    }

    // add the Boltzmann distribution to the initial state population
    double IExponent;
    if (temperature == 0)
      if (E_prime_prime == 0)
        IExponent = 0; // intensity unchanged
      else
        IExponent = 100; //(intensity < 10e-44 for non ground states
    else {
      IExponent = E_prime_prime / (temperature * KELVINS2EV);
      if (IExponent > 100)
        IExponent = 100; // keep the intensity >= 10e-44 == exp(-100)
    }
    dushinsky.getSpectrum().getSpectralPoint(pt).getIntensity() *=
        exp(-IExponent);

    dushinsky.getSpectrum().getSpectralPoint(pt).getEnergy() = energy;
    dushinsky.getSpectrum().getSpectralPoint(pt).getE_prime_prime() =
        E_prime_prime;

    // if intensity below the intensity threshold or energy above the threshold
    // -- do not print
    if ((dushinsky.getSpectrum().getSpectralPoint(pt).getIntensity() <
         job_config.get_intensity_thresh()) or
        (-(energy - E_prime_prime + elStates[targN].Energy()) >
         thresholds.target_eV()) or
        (E_prime_prime > thresholds.initial_eV())) {
      dushinsky.getSpectrum().getSpectralPoint(pt).setIfPrint(false);
      points_removed++;
    }
  }
  std::cout << "Done\n" << std::flush;

  if (dush_config.get_max_quanta_init() != 0) {
    if (points_removed > 0)
      std::cout << "  " << points_removed
                << " hot bands were removed from the spectrum\n";
    else
      std::cout << "All hot bands are above the intensity threshold\n";
    std::cout << "\n" << std::flush;
  }

  //--------------------------------------------------------------------------------
  // Print the updated spectrum:
  dushinsky.getSpectrum().Sort();
  std::cout << HorizontalLine << "\n";
  std::cout
      << "        Stick photoelectron spectrum (with Dushinsky rotations) \n";
  std::cout << HorizontalLine << "\n";
  if (elStates[targN].ifNMReorderedManually()) {
    std::cout
        << "\nWARNING! The normal modes of the target state were reordered!\n"
        << "         New order is used for the target state assignment.\n";
  }

  std::vector<int> nms_dushinsky = no_excite_subspace.get_active_subspace();
  if (nms_dushinsky.size() != n_molecular_normal_modes) {
    std::cout << "\nNOTE: only the following normal modes were excited "
                 "(\"excite subspace\"):\n  ";
    for (int nm = 0; nm < nms_dushinsky.size(); nm++)
      std::cout << nms_dushinsky[nm] << ' ';
    std::cout << "\n";
    if (elStates[targN].ifNMReorderedManually()) {
      std::cout << "\nWARNING! The normal modes of the target state were "
                   "reordered!\n"
                << "         New order is used for the \"excite subspace\"\n";
    }
  }
  std::cout << "\n";

  dushinsky.getSpectrum().PrintStickTable();

  // save the spectrum to the file
  std::string spectrumFName =
      InputFileName + std::string(".spectrum_dushinsky");
  dushinsky.getSpectrum().PrintStickTable(spectrumFName);
  std::cout << "\nStick spectrum was also saved in \"" << spectrumFName
            << "\" file \n";
  if (nms_dushinsky.size() != n_molecular_normal_modes)
    std::cout << " (Full list of the normal modes was used for assigning "
                 "transitions)\n";
  std::cout << "\n\n";
}
