#include "dushinsky.h"
#include "energy_thresholds.h"

/* This function contains exactly what the old constructor contained. There is
 * work needed here.
 *
 * This function evaluates the <0|0> integral and all required matrices for a
 * given set of normal modes. i.e. the layer L=0 and is ready for iterative
 * evaluation of next layers.
 *
 * "do_not_excite_subspace" in the old version was passed as an argument
 * ```nm_list``` which was a std::vector<int> that listed numbers of normal
 * modes (list of normal modes to include)
 *
 * molStates -- vector of molecular states loaded from the input file;
 * in_targN -- which target state to take from the molStates vector;
 *
 * IMPORTANT! All matrices and zero-zero integral are evaluated for the FULL
 * space of the normal modes of the molecule then space is shrinked to the
 * "active subspace" described by the ```DoNotExcite``` object. Excitations
 * will be added only to those normal modes.
 * */
void Dushinsky::old_constructor(std::vector<MolState> &molStates,
                                const int in_targN,
                                const DushinskyParameters &dush_parameters,
                                const JobParameters &job_parameters,
                                const DoNotExcite &no_excite_subspace) {
  Kp_max = 0;
  Kp_max_saved = 0;
  // fcf threshold (from the <job_parameters> tag)
  fcf_threshold = sqrt(job_parameters.get_intensity_thresh());

  // index of the inital electronic state
  int iniN = 0;
  targN = in_targN;

  // number of normal modes in the subset
  N = molStates[iniN].NNormModes();

  //----------------------------------------------------------------------
  // get frequencies, displacements and rotation matrix for the current
  // subspace:

  // Load frequencies:
  std::vector<double> selected_nm_freq_ini;
  selected_nm_freq_ini.clear();
  for (int nm = 0; nm < N; nm++)
    selected_nm_freq_ini.push_back(molStates[iniN].getNormMode(nm).getFreq());
  std::vector<double> selected_nm_freq_targ;
  selected_nm_freq_targ.clear();
  for (int nm = 0; nm < N; nm++)
    selected_nm_freq_targ.push_back(molStates[targN].getNormMode(nm).getFreq());

  // get w & w'; frequencies; p=='==target electronic state; [N];
  // frequency==w_i=nu_i[cm-1]*BOHR_2_CM*2*PI*SPEEDOFLIGHT_AU;
  arma::Col<double> w(N);
  for (int nm = 0; nm < N; nm++)
    w(nm) = selected_nm_freq_ini[nm] * BOHR_2_CM * 2 * PI * SPEEDOFLIGHT_AU;

  arma::Col<double> wp(N);
  for (int nm = 0; nm < N; nm++)
    wp(nm) = selected_nm_freq_targ[nm] * BOHR_2_CM * 2 * PI * SPEEDOFLIGHT_AU;

  // get normal modes: L and Lp (for the selected subspace):
  arma::Mat<double> L(CARTDIM * (molStates[iniN].NAtoms()), N,
                      arma::fill::zeros);
  for (int a = 0; a < molStates[iniN].NAtoms(); a++)
    for (int nm = 0; nm < N; nm++)
      for (int k = 0; k < CARTDIM; k++)
        L(a * CARTDIM + k, nm) =
            molStates[iniN].getNormMode(nm).getDisplacement()(a * CARTDIM + k);

  arma::Mat<double> Lp(CARTDIM * (molStates[targN].NAtoms()), N,
                       arma::fill::zeros);
  for (int a = 0; a < molStates[targN].NAtoms(); a++)
    for (int nm = 0; nm < N; nm++)
      for (int k = 0; k < CARTDIM; k++)
        Lp(a * CARTDIM + k, nm) =
            molStates[targN].getNormMode(nm).getDisplacement()(a * CARTDIM + k);

  // get S; [NxN]; S=(Lp^T)*L; n.m.rotation matrix; if S==I, than norm. modes
  // are parallel, det(S) sould be close to 1
  arma::Mat<double> S;
  S = Lp.t() * L;

  // get d displacements ; [N]; d= Lp^T*sqrt(T)*(x-xp); x & xp -- cartesian
  // geometries of two states; displacement of the target geometry in normal
  // coordinates  d[angstr*sqrt(amu)] *ANGSTROM2AU*sqrt(AMU_2_ELECTRONMASS) );
  // mass weighted -- so m=1amu everywhere;
  arma::Col<double> tmpX(CARTDIM * molStates[iniN].NAtoms());
  arma::Col<double> tmpXp(CARTDIM * molStates[iniN].NAtoms());
  for (int a = 0; a < molStates[iniN].NAtoms(); a++)
    for (int k = 0; k < CARTDIM; k++) {
      tmpX(a * CARTDIM + k) = molStates[iniN].getAtom(a).Coord(k);
      tmpXp(a * CARTDIM + k) = molStates[targN].getAtom(a).Coord(k);
    }
  //(x-xp)
  tmpX -= tmpXp;
  // Make sqrt(T)-matrix (diagonal matrix with sqrt(atomic masses) in cartesian
  // coordinates); the matrix is the same for all states;
  int mass_matrix_size = CARTDIM * (molStates[iniN].NAtoms());
  arma::Mat<double> SqrtT(mass_matrix_size, mass_matrix_size,
                          arma::fill::zeros);
  for (int a = 0; a < molStates[iniN].NAtoms(); a++) {
    double sqrt_of_mass = sqrt(molStates[iniN].getAtom(a).Mass());
    for (int loc_i = a * CARTDIM; loc_i < a * CARTDIM + 3; ++loc_i) {
      SqrtT(loc_i, loc_i) = sqrt_of_mass;
    }
  }
  // d = Lp^T * sqrt(T) * (x-xp);
  arma::Col<double> d;
  d = Lp.t() * SqrtT * tmpX;
  // switch to atomic units
  d *= ANGSTROM2AU * sqrt(AMU_2_ELECTRONMASS);

  //----------------------------------------------------------------------
  // now evaluate all the matrices:

  // temporary NxN matrices:
  arma::Mat<double> tmpM(N, N);

  // diag(1) NxN matrix:
  arma::Mat<double> I = arma::eye(N, N);

  // get Lm & Lmp; lamda & lamda'; [NxN] diag.; sqrt.freq.;
  // \lamda=diag(sqrt(w_i));
  arma::Mat<double> Lm(N, N, arma::fill::zeros);
  for (int nm = 0; nm < N; nm++)
    Lm(nm, nm) = sqrt(w(nm));

  arma::Mat<double> Lmp(N, N, arma::fill::zeros);
  for (int nm = 0; nm < N; nm++)
    Lmp(nm, nm) = sqrt(wp(nm));

  // get Dt; [N]; \delta=Lmp*d;
  arma::Col<double> Dt = Lmp * d;

  // get J;  [NxN]; J=Lmp*S*Lm^{-1};  ^{-1} -- inverse;
  arma::Mat<double> J(N, N);
  J = Lmp * S * Lm.i();

  // get Q;  [NxN] symm. pos.; Q = (1 + J^T * J)^{-1}; ^{T} -- transposed; ^{-1}
  // -- inverse;
  arma::Mat<double> Q(N, N);
  tmpM = I + J.t() * J;
  Q = tmpM.i();
  double detQ = arma::det(Q);

  // get P;  [NxN] symm.; P = J * Q * J^T
  arma::Mat<double> P(N, N);
  P = J * Q * J.t();

  // get R;  [NxN]; R = Q * J^T
  arma::Mat<double> R(N, N);
  R = Q * J.t();

  // get Det(S)
  double detS = arma::det(S);
  std::cout << "Determinant of the normal modes rotation matrix: |Det(S)| ="
            << std::fixed << std::setw(12) << std::setprecision(8) << fabs(detS)
            << "\n";
  if (fabs(detS) < 0.5) {
    std::cout << "\n"
              << "Error: |Det(S)| is too small (<0.5). Please see \"|Det(S)| "
                 "is less than one\n"
              << "   (or even zero)\" in the \"Common problems\" section of "
                 "the manual\n\n";
    exit(2);
  }

  //--------------------------------------------------------------------------------
  // zero_zero, K'=0: one element = <0|0>
  double zero_zero = 1;
  for (int nm = 0; nm < N; nm++)
    zero_zero *= w(nm) / wp(nm);
  // ORIGINAL: pow 0.25; MODIFIED TO: pow -0.25
  zero_zero =
      pow(2.0, N * 0.5) * pow(zero_zero, -0.25) * sqrt(detQ) / sqrt(fabs(detS));
  tmpM = -0.5 * Dt.t() * (I - P) * Dt; // It's a scalar now
  zero_zero *= exp(tmpM(0, 0));

  //--------------------------------------------------------------------------------
  // get the the frequently used vectors/matrices:

  // for target state excitations (no hot bands):

  // ompd=sqrt(2)*(1-P)*Dt; [N]
  ompd_full = sqrt(2) * (I - P) * Dt;

  // tpmo=(2P-1); [NxN]
  tpmo_full = 2 * P - I;

  // for the hot bands:

  // rd=sqrt(2)*R*Dt; [N]
  rd_full = sqrt(2) * R * Dt;

  // tqmo=(2Q-1); [NxN]
  tqmo_full = 2 * Q - I;

  // tr=2*R; [NxN]
  tr_full = 2 * R;

  //----------------------------------------------------------------------
  // shrink the space to the nm_list space
  const std::vector<int> active_nm_subspace =
      no_excite_subspace.get_active_subspace();

  if (N != active_nm_subspace.size()) {
    N = active_nm_subspace.size();

    // shrink "frequently used matrices"
    ompd = arma::Col<double>(N);
    rd = arma::Col<double>(N);
    tpmo = arma::Mat<double>(N, N);
    tqmo = arma::Mat<double>(N, N);
    tr = arma::Mat<double>(N, N);

    for (int nm1 = 0; nm1 < N; nm1++) {
      ompd(nm1) = ompd_full(active_nm_subspace[nm1]);
      rd(nm1) = rd_full(active_nm_subspace[nm1]);
      for (int nm2 = 0; nm2 < N; nm2++) {
        tpmo(nm1, nm2) =
            tpmo_full(active_nm_subspace[nm1], active_nm_subspace[nm2]);
        tqmo(nm1, nm2) =
            tqmo_full(active_nm_subspace[nm1], active_nm_subspace[nm2]);
        tr(nm1, nm2) =
            tr_full(active_nm_subspace[nm1], active_nm_subspace[nm2]);
      }
    }
  } else {
    ompd = ompd_full;
    tpmo = tpmo_full;
    rd = rd_full;
    tqmo = tqmo_full;
    tr = tr_full;
  }

  //----------------------------------------------------------------------
  // create the initial ground vibr. state <0000...000|
  state0.clear();
  state0.setElStateIndex(iniN);
  for (int nm = 0; nm < active_nm_subspace.size(); nm++)
    state0.addVibrQuanta(0, active_nm_subspace[nm]);

  // create a target ground vibr. state |0000...000>
  state.clear();
  state.setElStateIndex(targN);
  for (int nm = 0; nm < active_nm_subspace.size(); nm++)
    state.addVibrQuanta(0, active_nm_subspace[nm]);

  //----------------------------------------------------------------------
  // add the zero-th layer, K=0:
  std::vector<double> *layer_ptr = new std::vector<double>;
  (*layer_ptr).push_back(zero_zero);
  layers.push_back(layer_ptr);
  if (fabs(zero_zero) > fcf_threshold)
    // set energies later (for now 0.0):
    addSpectralPoint(zero_zero, state0, state);

  // combinations nChoosek(n,k) = C_n^k = n!/(k!*(n-k)!);
  // n=0..N+K-1; k=0..N-1=0..min(n,N-1);
  // N-total number of normal modes; K -- max number of quanta in target state

  int max_quanta_initial = dush_parameters.get_max_quanta_init();
  int max_quanta_target = dush_parameters.get_max_quanta_targ();

  int K = max_quanta_target;
  int size = (N) * (N + K);

  // Create an array of sqrt() 0..K+1
  K = std::max(max_quanta_target, max_quanta_initial);
  sqrtArray = new double[K + 2];
  for (int i = 0; i < K + 1; i++)
    sqrtArray[i] = sqrt(i);
}

/* Go over all layers and add points to the spectrum. */
void Dushinsky::evaluate_higher_levels(
    const DushinskyParameters &dush_parameters) {

  int max_quanta_targ = dush_parameters.get_max_quanta_targ();
  int Kp_max_to_save = dush_parameters.get_Kp_max_to_save();

  for (int Kp = 1; Kp <= max_quanta_targ; Kp++) {
    bool save_layer = Kp <= Kp_max_to_save;
    if (save_layer) {
      std::cout << "Layer K'=" << Kp
                << " is being evaluated (will be saved in memory)... "
                << std::flush;
    } else {
      std::cout << "Layer K'=" << Kp << " is being evaluated... " << std::flush;
    }

    int n_fresh_points = evalNextLayer(save_layer);

    std::cout << "Done\n";
    if (n_fresh_points > 0) {
      std::cout << n_fresh_points
                << " points above the intensity threhold were added to the "
                   "spectrum\n\n"
                << std::flush;
    } else {
      std::cout << "No points above the intensity threhold were found in "
                   "this layer\n\n"
                << std::flush;
    }
  }
}

/* Creator of the Dushinsky class that is responsible for calculating the FCFs
 * with inclusion of the Duschinsky rotation. */
Dushinsky::Dushinsky(std::vector<MolState> &molStates, const int in_targN,
                     const DushinskyParameters &dush_parameters,
                     const JobParameters &job_parameters,
                     const DoNotExcite &no_excite_subspace) {
  old_constructor(molStates, in_targN, dush_parameters, job_parameters,
                  no_excite_subspace);
  printLayersSizes(dush_parameters, no_excite_subspace);
  evaluate_higher_levels(dush_parameters);
}

Dushinsky::~Dushinsky() {
  for (int Kp = 0; Kp < layers.size(); Kp++)
    delete layers[Kp];
  layers.clear();

  // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
  // delete [] C;
  delete[] sqrtArray;
}

int Dushinsky::evalNextLayer(const bool if_save) {
  // debug check:
  if (if_save and (Kp_max) != (Kp_max_saved)) {
    std::cerr << "\nDebug Error!: Layer " << Kp_max + 1
              << " can not be saved. The number of already saved layers is "
              << Kp_max_saved << " but should be " << Kp_max << "\n\n";
    exit(2);
  }

  // points above the threshold
  int points_added = 0;

  // reset the target state
  for (int i = 0; i < N; i++) {
    state.getV()[i] = -1;
  }

  // create a new layer
  std::vector<double> *layer_ptr = new std::vector<double>;

  // check the memory avaliability (dirty):
  if (if_save) {
    unsigned long total_combs =
        nChoosek((Kp_max + 1 + state.getVibrQuantaSize() - 1),
                 (state.getVibrQuantaSize() - 1));
    double *buffer = (double *)malloc(total_combs * sizeof(double));
    if (buffer == NULL) {
      std::cout << "\n\nError: not enough memory available to store layer K'="
                << Kp_max_saved + 1 << "\n\n"
                << "Please reduce \"max_vibr_to_store\" value to limit the "
                   "memory use.\n"
                << " You also can reduce "
                   "\"max_vibr_excitations_in_target_el_state\" value,\n"
                << "and/or use \"single_excitation\" elements to add higher "
                   "excitations.\n";
      exit(2);
    }
    free(buffer);
  }

  unsigned long index_counter = 0;
  while (enumerateVibrStates(state.getVibrQuantaSize(), Kp_max + 1,
                             state.getV(), true)) {

    // ZZZ 4/11/2012 removed, and the combinations are calculated now on the fly
    // unsigned long index_rev=convVibrState2Index(state.getV(), N, C,
    // Kp_max+1);
    unsigned long index_rev = convVibrState2Index(state.getV(), N, Kp_max + 1);

    // check if the reverse index is ok
    if (index_rev != index_counter) {
      // if numbers are too large, factorials in nChoosek() function will be out
      // of "unsigned long" range ...
      std::cout << "\n Error!\n[Debug info: reverse index function "
                   "convVibrState2Index(state.getV()) for the state:\n";
      state.print();
      std::cout
          << "\n"
          << "returns index=" << index_rev
          << "; should be index=" << index_counter << "]\n\n"
          << "Layer #" << Kp_max
          << " is the maximum layer which can be handled for this system.\n"
          << "Please reduce \"max_vibr_excitations_in_target_el_state\" value "
             "to "
          << Kp_max << ".\n"
          << "You can also use \"single_excitation\" elements to add higher "
             "excitations manually.\n";
      exit(2);
    }

    double fcf = evalSingleFCF(state0, 0, state, Kp_max + 1);

    if (if_save)
      (*layer_ptr).push_back(fcf);

    if (fabs(fcf) > fcf_threshold) {
      points_added++;
      addSpectralPoint(fcf, state0, state);
    }

    index_counter++;
  }

  Kp_max++;
  // save the layer
  if (if_save) {
    layers.push_back(layer_ptr);
    Kp_max_saved++;
  }

  return points_added;
}

int Dushinsky::add_the_only_intial_state_transitions(
    const int Kp, VibronicState &the_only_initial_state) {
  // This is basically a slightly eddited version of Duschinsky::evalNextLayer()

  // points above the threshold
  int points_added = 0;

  // reset the target state
  for (int i = 0; i < N; i++) {
    state.getV()[i] = -1;
  }

  // create a new layer
  std::vector<double> *layer_ptr = new std::vector<double>;

  unsigned long index_counter = 0;
  while (enumerateVibrStates(state.getVibrQuantaSize(), Kp_max, state.getV(),
                             true)) {
    unsigned long index_rev = convVibrState2Index(state.getV(), N, Kp_max);

    // check if the reverse index is ok
    if (index_rev != index_counter) {
      // if numbers are too large, factorials in nChoosek() function will be out
      // of "unsigned long" range ...
      std::cout << "\n Error!\n[Debug info: reverse index function "
                   "convVibrState2Index(state.getV()) for the state:\n";
      state.print();
      std::cout
          << "\n"
          << "returns index=" << index_rev
          << "; should be index=" << index_counter << "]\n\n"
          << "Layer #" << Kp_max
          << " is the maximum layer which can be handled for this system.\n"
          << "Please reduce \"max_vibr_excitations_in_target_el_state\" value "
             "to "
          << Kp_max << ".\n"
          << "You can also use \"single_excitation\" elements to add higher "
             "excitations manually.\n";
      exit(2);
    }

    // evaluate FCF for each transition and add to the spectrum:
    int K = the_only_initial_state.getTotalQuantaCount();
    int Kp = state.getTotalQuantaCount();

    double fcf = evalSingleFCF_full_space(the_only_initial_state, K, state, Kp);

    if (fabs(fcf) > fcf_threshold) {
      points_added++;
      addSpectralPoint(fcf, the_only_initial_state, state);
    }

    index_counter++;
  }

  Kp_max++;

  return points_added;
}

double Dushinsky::evalSingleFCF(VibronicState &state_ini, int K,
                                VibronicState &state_targ, int Kp) {
  double fcf;

  if (K == 0) // if there are no excitations in the initial state left

    if (Kp <=
        Kp_max_saved) // if there are less than K'max excitations in the target
                      // state left (i.e. saved) ZZZ 4/11/2012 removed, and the
                      // combinations are calculated now on the fly
                      // fcf=(*layers[Kp])[convVibrState2Index(state_targ.getV(),
                      // N, C, Kp)];
      fcf = (*layers[Kp])[convVibrState2Index(state_targ.getV(), N, Kp)];

    else {
      // find first non zero quanta
      int ksi = 0;
      while (state_targ.getV()[ksi] == 0)
        ksi++;

      // get the first term
      state_targ.getV()[ksi]--;
      fcf = ompd[ksi] * evalSingleFCF(state_ini, K, state_targ, Kp - 1);

      // add the second term for K'>=2
      if (Kp > 1)
        for (int theta = ksi; theta < N; theta++)
          if (state_targ.getV()[theta] > 0) {
            double tmp_dbl =
                tpmo(ksi, theta) * sqrtArray[state_targ.getV()[theta]];
            state_targ.getV()[theta]--;
            tmp_dbl *= evalSingleFCF(state_ini, K, state_targ, Kp - 2);
            fcf += tmp_dbl;
            state_targ.getV()[theta]++;
          }
      fcf /= sqrtArray[state_targ.getV()[ksi] + 1];

      // get back to the original state ("to be calculated" state)
      state_targ.getV()[ksi]++;
    }
  else {
    // find first non zero quanta
    int ksi = 0;
    while (state_ini.getV()[ksi] == 0)
      ksi++;

    // get the first term
    state_ini.getV()[ksi]--;
    fcf = -rd(ksi) * evalSingleFCF(state_ini, K - 1, state_targ, Kp);

    // add the second term
    if (K > 1)
      for (int theta = ksi; theta < N; theta++)
        if (state_ini.getV()[theta] > 0) {
          double tmp_dbl =
              tqmo(ksi, theta) * sqrtArray[state_ini.getV()[theta]];
          state_ini.getV()[theta]--;
          tmp_dbl *= evalSingleFCF(state_ini, K - 2, state_targ, Kp);
          fcf += tmp_dbl;
          state_ini.getV()[theta]++;
        }

    // add the third term
    if (Kp > 0)
      for (int theta = 0; theta < N; theta++)
        if (state_targ.getV()[theta] > 0) {
          double tmp_dbl = tr(ksi, theta) * sqrtArray[state_targ.getV()[theta]];
          state_targ.getV()[theta]--;
          tmp_dbl *= evalSingleFCF(state_ini, K - 1, state_targ, Kp - 1);
          fcf += tmp_dbl;
          state_targ.getV()[theta]++;
        }

    // normolize
    fcf /= sqrtArray[state_ini.getV()[ksi] + 1];

    // get back to the original state ("to be calculated" state)
    state_ini.getV()[ksi]++;
  }

  return fcf;
}

double Dushinsky::evalSingleFCF_full_space(VibronicState &state_ini, int K,
                                           VibronicState &state_targ, int Kp) {
  double fcf;

  if (K == 0)    // if there are no excitations in the initial state left
    if (Kp == 0) // if there are no excitations in the target state left
      fcf = (*layers[0])[0]; // return <0|0> integral
    else {
      // find first non zero quanta
      int ksi = 0;
      while (state_targ.getV()[ksi] == 0)
        ksi++;

      // get the first term
      state_targ.getV()[ksi]--;
      fcf = ompd_full[ksi] *
            evalSingleFCF_full_space(state_ini, K, state_targ, Kp - 1);

      // add the second term for K'>=2
      if (Kp > 1)
        for (int theta = ksi; theta < N; theta++)
          if (state_targ.getV()[theta] > 0) {
            double tmp_dbl =
                tpmo_full(ksi, theta) * sqrt(state_targ.getV()[theta]);
            state_targ.getV()[theta]--;
            tmp_dbl *=
                evalSingleFCF_full_space(state_ini, K, state_targ, Kp - 2);
            fcf += tmp_dbl;
            state_targ.getV()[theta]++;
          }
      fcf /= sqrt(state_targ.getV()[ksi] + 1);

      // get back to the original state ("to be calculated" state)
      state_targ.getV()[ksi]++;
    }
  else {
    // find first non zero quanta
    int ksi = 0;
    while (state_ini.getV()[ksi] == 0)
      ksi++;

    // get the first term
    state_ini.getV()[ksi]--;
    fcf = -rd_full[ksi] *
          evalSingleFCF_full_space(state_ini, K - 1, state_targ, Kp);

    // add the second term
    if (K > 1)
      for (int theta = ksi; theta < N; theta++)
        if (state_ini.getV()[theta] > 0) {
          double tmp_dbl =
              tqmo_full(ksi, theta) * sqrt(state_ini.getV()[theta]);
          state_ini.getV()[theta]--;
          tmp_dbl *= evalSingleFCF_full_space(state_ini, K - 2, state_targ, Kp);
          fcf += tmp_dbl;
          state_ini.getV()[theta]++;
        }

    // add the third term
    if (Kp > 0)
      for (int theta = 0; theta < N; theta++)
        if (state_targ.getV()[theta] > 0) {
          double tmp_dbl = tr_full(ksi, theta) * sqrt(state_targ.getV()[theta]);
          state_targ.getV()[theta]--;
          tmp_dbl *=
              evalSingleFCF_full_space(state_ini, K - 1, state_targ, Kp - 1);
          fcf += tmp_dbl;
          state_targ.getV()[theta]++;
        }
    fcf /= sqrt(state_ini.getV()[ksi] + 1);

    // get back to the original state ("to be calculated" state)
    state_ini.getV()[ksi]++;
  }

  return fcf;
}

int Dushinsky::addHotBands(std::vector<MolState> &molStates,
                           std::vector<int> &nm_list, double fcf_threshold,
                           double temperature, int max_n_initial,
                           int max_n_target,
                           const EnergyThresholds &thresholds) {

  int points_added = 0;

  // vector of states below the energy thresholds and with total number of
  // excitations up to requested number
  std::vector<VibronicState> selected_states_ini, selected_states_targ;

  // reset the initial state
  for (int i = 0; i < N; i++)
    state0.getV()[i] = -1;
  // reset the target state
  for (int i = 0; i < N; i++)
    state.getV()[i] = -1;

  unsigned long total_combs_ini = 0, total_combs_targ = 0;
  for (int curr_max_ini = 0; curr_max_ini <= max_n_initial; curr_max_ini++)
    total_combs_ini +=
        nChoosek((curr_max_ini + nm_list.size() - 1), (nm_list.size() - 1));
  for (int curr_max_targ = 0; curr_max_targ <= max_n_target; curr_max_targ++)
    total_combs_targ +=
        nChoosek((curr_max_targ + nm_list.size() - 1), (nm_list.size() - 1));
  std::cout << "Maximum number of combination bands = "
            << total_combs_ini * total_combs_targ
            << "\n   = " << total_combs_ini
            << " (# of vibrational states in the initial electronic state)"
            << "\n   * " << total_combs_targ
            << " (# of vibrational states in the target electronic state)\n\n"
            << std::flush;

  // find INITIAL states with up to 'max_n_initial' vibrational quanta and with
  // energy below 'energy_threshold_initial':
  std::cout << "A set of initial vibrational states is being created...\n";
  if (thresholds.present()){
    std::cout << "  energy threshold = " << std::fixed
              << thresholds.initial_eV() << " eV ("
              << thresholds.initial_eV() / KELVINS2EV << " K)\n"
              << std::flush;
  }
  else {
    std::cout
        << "  energy threshold is not specified in the input (Please consider "
           "adding\n"
        << "  the \"energy_thresholds\" tag for a faster calculation.)\n\n";
  }

  // std::cout <<"NOTE: ezSpectrum may crash at this point if memory is
  // insufficient to store\n"
  //	    <<"      all vibrational states. If so, please reduce the initial
  // state's energy\n"
  //	    <<"      threshold or max_vibr_excitations_in_initial_el_state\n\n"
  //<< std::flush;

  // start with one quanta in the initial state (hot bands)
  for (int curr_max_ini = 1; curr_max_ini <= max_n_initial; curr_max_ini++) {
    while (enumerateVibrStates(state0.getVibrQuantaSize(), curr_max_ini,
                               state0.getV(), true)) {
      double energy = 0;

      for (int nm = 0; nm < molStates[state0.getElStateIndex()].NNormModes();
           nm++)
        energy +=
            molStates[state0.getElStateIndex()].getNormMode(nm).getFreq() *
            WAVENUMBERS2EV * state0.getV_full_dim(nm);

      if (energy < thresholds.initial_eV()) {
        // check memory available (dirty, but anyway one copying of the state is
        // requared...)
        // VibronicState * state_tmp = (VibronicState*) malloc (
        // sizeof(VibronicState)*2 ); if (state_tmp==NULL)
        //{
        //  std::cout << "\nError: not enough memory available to store all
        //  initial vibrational states\n\n"; exit(2);
        //}
        // free (state_tmp);

        selected_states_ini.push_back(state0);
      }
    }
    // reset the initial state's vibration "population"
    state0.getV()[0] = -1;
  }
  std::cout << "  " << selected_states_ini.size()
            << " vibrational states below the energy threshold\n\n"
            << std::flush;

  // find TARGET states with up to 'max_n_target' vibrational quanta and with
  // energy below 'energy_threshold_target':
  std::cout << "A set of target vibrational states is being created...\n";
  if (thresholds.present()) {
    std::cout << "  energy threshold = " << std::fixed
              << thresholds.target_eV() << " eV ("
              << thresholds.target_eV() / WAVENUMBERS2EV << " cm-1)\n"
              << std::flush;
  }
  else {
    std::cout
        << "  energy thresholds are not specified in the input (Please "
           "consider adding\n"
        << "  the \"energy_thresholds\" tag for a faster calculation.)\n\n";
  }

  // std::cout <<"NOTE: ezSpectrum may crash at this point if memory is
  // insufficient to store\n"
  //	    <<"      all vibrational states. If so, please reduce the target
  // state's energy\n"
  //	    <<"      threshold or max_vibr_excitations_in_target_el_state\n\n"
  //<< std::flush;

  for (int curr_max_targ = 0; curr_max_targ <= max_n_target; curr_max_targ++) {
    while (enumerateVibrStates(state.getVibrQuantaSize(), curr_max_targ,
                               state.getV(), true)) {
      double energy = 0;

      for (int nm = 0; nm < molStates[state.getElStateIndex()].NNormModes();
           nm++)
        energy += molStates[state.getElStateIndex()].getNormMode(nm).getFreq() *
                  WAVENUMBERS2EV * state.getV_full_dim(nm);

      if (energy < thresholds.target_eV()) {
        // check memory available (dirty, but anyway one copying of the state is
        // requared...)
        // VibronicState * state_tmp = (VibronicState*) malloc (
        // sizeof(VibronicState)*2 ); if (state_tmp==NULL)
        //	{
        //	  std::cout << "\nError: not enough memory available to store
        // all target vibrational states\n\n"; 	  exit(2);
        //	}
        // free (state_tmp);

        selected_states_targ.push_back(state);
      }
    }
    // reset the target state's vibration "population"
    state.getV()[0] = -1;
  }
  std::cout << "  " << selected_states_targ.size()
            << " vibrational states below the energy threshold\n\n";

  std::cout << "Total number of combination bands with thresholds applied: "
            << selected_states_ini.size() * selected_states_targ.size()
            << "\n\n"
            << std::flush;

  std::cout << "Hot bands are being calculated..." << std::flush;

  // evaluate all possible cross integrals
  for (int curr_ini = 0; curr_ini < selected_states_ini.size(); curr_ini++)
    for (int curr_targ = 0; curr_targ < selected_states_targ.size();
         curr_targ++) {
      int K = selected_states_ini[curr_ini].getTotalQuantaCount();
      int Kp = selected_states_targ[curr_targ].getTotalQuantaCount();

      double s_fcf = evalSingleFCF(selected_states_ini[curr_ini], K,
                                   selected_states_targ[curr_targ], Kp);

      if (fabs(s_fcf) > fcf_threshold) {
        points_added++;
        addSpectralPoint(s_fcf, selected_states_ini[curr_ini],
                         selected_states_targ[curr_targ]);
      }
    }

  std::cout << "Done\n\n" << std::flush;

  return points_added;
}

void Dushinsky::addSpectralPoint(const double fcf, VibronicState state_ini,
                                 VibronicState state_targ) {
  // set energies later
  spectrum.AddSpectralPoint(0.0, fcf * fcf, fcf, 0.0, state_ini, state_targ);
};

void Dushinsky::printLayersSizes(const DushinskyParameters &dush_parameters,
                                 const DoNotExcite &no_excite_subspace) {

  std::cout << "Number of normal modes active in this calculation: "
            << no_excite_subspace.get_active_subspace().size() << "\n\n";
  std::cout << "Size of layers with exactly K' excitations in the target "
               "state (in bytes):\n";

  // print estimated size of each layer up to K'_max:
  int max_quanta_ini = dush_parameters.get_max_quanta_init();
  int max_quanta_targ = dush_parameters.get_max_quanta_targ();
  int Kp_max_to_save = dush_parameters.get_Kp_max_to_save();

  const int uptoKp =
      (max_quanta_targ < Kp_max_to_save) ? max_quanta_targ : Kp_max_to_save;

  unsigned long elements_per_layer, size_per_layer = 0, size_per_layer_prev = 0;

  for (int Kp = 0; Kp <= uptoKp; Kp++) {
    elements_per_layer = nChoosek(Kp + N - 1, N - 1);
    size_per_layer_prev = size_per_layer;
    size_per_layer = elements_per_layer * sizeof(double);

    if (size_per_layer < size_per_layer_prev) {
      std::cout << "\nError! Layer #" << Kp - 1
                << " is the maximum layer which can be saved "
                << "for this molecule. Please use \"max_vibr_to_store\" or/and"
                << "\"max_vibr_excitations_in_target_el_state\" value to "
                << Kp - 1 << ".\n";
      exit(2);
    }

    std::cout << "   layer K'=" << Kp << ": ";
    std::cout.precision(2);
    if (size_per_layer < 1024)
      std::cout << size_per_layer;
    else if (size_per_layer < 1024 * 1024)
      std::cout << size_per_layer / 1024.0 << " K";
    else if (size_per_layer < 1024 * 1024 * 1024)
      std::cout << size_per_layer / (1024.0 * 1024.0) << " M";
    else
      std::cout << size_per_layer / (1024.0 * 1024.0 * 1024.0) << " G";
    std::cout << "\n";
  }

  std::cout
      << "Please be sure that you have enough memory to store all these "
         "layers,\n"
      << "otherwise use \"max_vibr_to_store\" tag to limit the memory use.\n"
      << "You also can reduce \"max_vibr_excitations_in_target_el_state\" "
         "value,\n"
      << "and/or use \"single_excitation\" elements to add higher "
         "excitations.\n";

  std::cout << "\n";
}
