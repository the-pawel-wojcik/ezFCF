#include "parallel_approximation.h"

Parallel::Parallel(std::vector <MolState>& molStates, std::vector<int>& nm_list, 
    double fcf_threshold, double temperature, 
    int max_n_initial, int max_n_target, 
    bool if_comb_bands, bool if_use_target_nm, bool if_print_fcfs, bool if_web_version, const char* nmoverlapFName,
    double energy_threshold_initial,  double energy_threshold_target)
{
  double intens_threshold=fcf_threshold*fcf_threshold;

  //initial state index
  int  iniN=0;

  // full space size
  int n_norm_modes = molStates[iniN].NNormModes();

  //
  std::vector <int> state_ini, state_targ, state_ini_subspace, state_targ_subspace;
  for (int i=0; i<n_norm_modes; i++)
  {
    state_ini.push_back(-1);
    state_targ.push_back(-1);
  }

  //"excite subspace" size
  for (int i=0; i<nm_list.size(); i++) 
  {
    state_ini_subspace.push_back(-1);
    state_targ_subspace.push_back(-1);
  }

  // stores set of the inital and target vibrational states (for each electronic state) -- after energy threshold applied
  std::vector < std::vector <int> > selected_states_ini, selected_states_targ;

  // matrix with FC factors (use the same variable for every target state)
  // max_n_initial stands for maximum allowed excitations in a single mode in the initial state
  // FCFs_tmp(n, m) gives the integral of an overlap of a state with 'n' excitation with 
  // a state of 'm' excitations, i.e, <n|m> 
  arma::Mat<double> FCFs_tmp(max_n_initial+1, max_n_target+1);
  std::vector <arma::Mat<double>> FCFs;

  // matrix with intensities of FC transitions = FCFs * population_of_initial_vibrational_levels(Temperature distrib.)
  arma::Mat<double> I_tmp (max_n_initial+1, max_n_target+1);
  std::vector <arma::Mat<double>> I;

  // energy position of each transition for a given normal mode; 1D; ofset by IP; when add for N dimensions, substract IP from each energy;
  arma::Mat<double> E_position_tmp (max_n_initial+1, max_n_target+1);
  std::vector <arma::Mat<double>> E_position;

  // reduced mass for FCF calcualtions -- mass weighted coordinates
  // TODO: This variable is not used in this scope. Is it still necessary? Paweł Apr'22
  double reducedMass=1.0;

  // ==========================================================================================
  // Calculate transformation matrix "cartesian->normal mode" coordinates (for each state):
  std::vector <arma::Mat<double>> Cart2Normal_Rs;

  for (int state=0; state<molStates.size(); state++)
  {
    arma::Mat<double> NormModes( CARTDIM*(molStates[state].NAtoms()), molStates[state].NNormModes(), arma::fill::zeros ); // NOT mass weighted normal modes i.e. (L~)=T^(0.5)*L
    //Get L (normal modes in cartesian coordinates mass unweighted (in Angstroms) ):
    for (int j=0; j < molStates[state].NAtoms(); j++) 
      for (int i = 0; i < molStates[state].NNormModes(); i++) 
        for (int k=0; k < CARTDIM; k++)
          NormModes(j*CARTDIM+k, i) = molStates[state].getNormMode(i).getDisplacement()(j*CARTDIM+k);

    //Make sqrt(T)-matrix (diagonal matrix with sqrt(atomic masses) in cartesian coordinates):
    arma::Mat<double> SqrtT( CARTDIM*(molStates[state].NAtoms()), CARTDIM*(molStates[state].NAtoms()), arma::fill::zeros);
    for (int i=0; i<molStates[state].NAtoms(); i++)
      SqrtT(i*CARTDIM, i*CARTDIM) 
        = SqrtT(i*CARTDIM+1, i*CARTDIM+1)
        = SqrtT(i*CARTDIM+2, i*CARTDIM+2)
        = sqrt(molStates[state].getAtom(i).Mass());

    // Cart->NormalModes transformation matrix R=L^T*sqrt(T)
    // q' = q+d = q+R*(x-x') (for the parallel normal mode approximation)
    // units of d are Angstr*sqrt(amu)
    NormModes = NormModes.t(); 
    NormModes *= SqrtT;

    Cart2Normal_Rs.push_back(NormModes);
  }

  // Initial geometry in cartesian coordinates:
  arma::Col<double> InitialCartCoord(CARTDIM*(molStates[iniN].NAtoms()));
  for (int i=0; i<molStates[iniN].NAtoms(); i++)
    for (int k=0; k<CARTDIM; k++)
      InitialCartCoord[i*CARTDIM+k] = molStates[iniN].getAtom(i).Coord(k);
  InitialCartCoord *= ANGSTROM2AU; 

  // initialize spectrlPoint: create vector of VibrQuantNumbers for the initital and target state 
  // (keep only non-zero qunta, i.e. from the "excite subspace")
  // Also add subspace mask (the rest of nmodes do not have excitations, and  would not be printers in the spectrum)
  SpectralPoint tmpPoint;
  for (int nm=0; nm<nm_list.size(); nm++)
  {
    tmpPoint.getVibrState1().addVibrQuanta(0, nm_list[nm]);
    tmpPoint.getVibrState2().addVibrQuanta(0, nm_list[nm]);
  }

  // For each target state, for each normal mode: calculate shift, FCFs, and print spectrum
  for (int targN=1; targN < molStates.size(); targN++)
  {
    std::cout << "===== Target state #" << targN << " =====\n\n"<< std::flush;

    // clear the FC factors and related data
    FCFs.clear();
    I.clear();
    E_position.clear();

    //----------------------------------------------------------------------
    // calculate dQ:

    // Target state geometry in cartesian coordinates:
    arma::Col<double> TargetCartCoord(CARTDIM*(molStates[iniN].NAtoms()));
    for (int i=0; i<molStates[targN].NAtoms(); i++)
      for (int k=0; k <CARTDIM; k++)
        TargetCartCoord[i*CARTDIM+k] = molStates[targN].getAtom(i).Coord(k);

    // input is in angstroms
    TargetCartCoord *= ANGSTROM2AU;

    // Shift of the target state relative to the initial (In cart. coord.)
    arma::Col<double> cartesian_shift = TargetCartCoord - InitialCartCoord;

    // Transform Cart->normal Mode coordinates

    // Geometry differences in normal coordinates
    arma::Col<double> NormModeShift_ini(molStates[iniN].NNormModes());
    arma::Col<double> NormModeShift_targ(molStates[iniN].NNormModes());

    NormModeShift_ini = Cart2Normal_Rs[iniN] * cartesian_shift;
    NormModeShift_targ = Cart2Normal_Rs[targN] * cartesian_shift;

    // Rescale it back & print:
    NormModeShift_ini *= AU2ANGSTROM;
    NormModeShift_targ *= AU2ANGSTROM;

    std::cout << "Difference (dQ) between the initial and the target state geometries.\n"
      << "Angstrom*sqrt(amu):\n\n"
      << "normal mode  dQ in initial  dQ in target   frequency   frequency   comments\n"
      << "  number      state coord.  state coord.    initial      target\n\n";
    for (int nm=0; nm<n_norm_modes; nm++)
    {
      std::cout << "     " << std::fixed << std::setw(3) << nm <<"      " 
        << std::setprecision(6) 
        << std::setw(9) <<  NormModeShift_ini[nm] << "       " 
        << std::setw(9) << NormModeShift_targ[nm] << "      " 
        << std::setprecision(2) 
        << std::setw(7) << molStates[iniN].getNormMode(nm).getFreq() << "    " 
        << std::setw(7) << molStates[targN].getNormMode(nm).getFreq() << "\n";
    }
    std::cout<<"\n\n";


    //----------------------------------------------------------------------
    // save the spectrumnm overlap to file (for normal mode reordering tool in the web interface)
    if (if_web_version)
    {
      std::vector <int> nondiagonal_list;
      arma::Mat<double> NMoverlap;
      bool if_overlap_diagonal = molStates[targN].getNormalModeOverlapWithOtherState(molStates[iniN], NMoverlap, nondiagonal_list);

      std::ofstream nmoverlapF;     
      nmoverlapF.open(nmoverlapFName, std::ios::out);
      nmoverlapF << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
        << "<nmoverlap\n  nm_initial=\"" <<n_norm_modes <<"\"\n  nm_target =\""<<n_norm_modes << "\">\n\n<nm_order_initial>\n";
      //inital state's nm order -- 0,1,2... -- always!
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<oi" << nm << ">"<< molStates[iniN].getNormModeIndex(nm) << "</oi" << nm << ">\n";
      nmoverlapF << "</nm_order_initial>\n\n<frequencies_initial>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<fi" << nm << ">"<< std::fixed << std::setprecision(1) << molStates[iniN].getNormMode(nm).getFreq()<< "</fi" << nm << ">\n";
      nmoverlapF << "</frequencies_initial>\n\n<displacements_initial>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<dqi" << nm << ">"<< std::fixed << std::setprecision(2) << NormModeShift_ini[nm]<< "</dqi" << nm << ">\n";
      nmoverlapF << "</displacements_initial>\n\n<nm_order_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<ot" << nm << ">"<< molStates[targN].getNormModeIndex(nm) << "</ot" << nm << ">\n";
      nmoverlapF << "</nm_order_target>\n\n<frequencies_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<ft" << nm << ">"<< std::fixed << std::setprecision(1) << molStates[targN].getNormMode(nm).getFreq()<< "</ft" << nm << ">\n";
      nmoverlapF << "</frequencies_target>\n\n<displacements_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<dqt" << nm << ">"<< std::fixed << std::setprecision(2) << NormModeShift_targ[nm]<< "</dqt" << nm << ">\n";
      nmoverlapF << "</displacements_target>\n\n<overlap_matrix>\n";
      for (int nmt=0; nmt<n_norm_modes; nmt++)
      {    
        nmoverlapF << "<row"<<nmt<<">\n";
        for (int nmi=0; nmi<n_norm_modes; nmi++)
          nmoverlapF << "<c"<<nmi<<">"<< NMoverlap(nmt, nmi)<<"</c"<<nmi<<">";
        nmoverlapF << "\n</row"<<nmt<<">\n";
      }
      nmoverlapF << "</overlap_matrix>\n\n</nmoverlap>";
      nmoverlapF.close();
    }
    //----------------------------------------------------------------------


    if (if_use_target_nm) 
      std::cout<<"TARGET state normal coordinates (displacements dQ) will be used.\n\n"; 
    else 
      std::cout<<"INITIAL state normal coordinates are used.\n\n"; 
    // ----------------------------------------------------------------------

    if (if_print_fcfs)
      std::cout << "Peak positions, Franck-Condon factors, and intensities along each normal mode:\n\n";

    // for each normal mode of the initial state
    // TODO: This appears as too 
    for (int nm=0; nm<n_norm_modes; nm++) 
    {      
      // Calculate the energy positions:
      for (int i=0; i<max_n_initial+1; i++)
        for (int j=0; j<max_n_target+1; j++)
          // energy position of each transition for a given normal mode; 1D; ofset by IP; when add for N dimensions, substract IP from each energy;
          E_position_tmp(i, j)
            =
            - molStates[targN].Energy()
            + (molStates[iniN].getNormMode(nm).getFreq() * i
                - molStates[targN].getNormMode(nm).getFreq() * j) * WAVENUMBERS2EV;
      E_position.push_back(E_position_tmp);

      // Calculate the Franck-Condon factors
      if (if_use_target_nm) // which normal modes (displacements) to use: of the initial or the target state?
        harmonic_FCf(FCFs_tmp, reducedMass, NormModeShift_targ[nm],
            molStates[iniN].getNormMode(nm).getFreq(), 
            molStates[targN].getNormMode(nm).getFreq());
      else
        harmonic_FCf(FCFs_tmp, reducedMass, NormModeShift_ini[nm],
            molStates[iniN].getNormMode(nm).getFreq(), 
            molStates[targN].getNormMode(nm).getFreq());

      FCFs.push_back(FCFs_tmp);

      // Calculate the line intensities (FCF^2 *  the temperature distribution in the initial state):
      for (int i=0; i<max_n_initial+1; i++)
      { 
        double IExponent=100*i; //(intensity < 10e-44 for non ground states if temperature=0)
        if (temperature>0)
        {
          IExponent = molStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * i 
            / (temperature * KELVINS2EV);
          if (IExponent > 100) 
            IExponent=100; // keep the intensity >= 10e-44 == exp(-100)
        }

        for (int j=0; j<max_n_target+1; j++)
          I_tmp(i, j) = FCFs_tmp(i, j) * FCFs_tmp(i, j) * exp(-IExponent);
      }
      I.push_back(I_tmp);

      if (if_print_fcfs)
      {
        std::cout << "NORMAL MODE #" << nm
          << " (" 
          << std::fixed << std::setw(8) << std::setprecision(2)
          << molStates[iniN].getNormMode(nm).getFreq() << "cm-1 --> " 
          << molStates[targN].getNormMode(nm).getFreq() << "cm-1, dQ="
          << std::setw(6) << std::setprecision(4);
        if (if_use_target_nm)
          std::cout << NormModeShift_targ[nm];
        else
          std::cout << NormModeShift_ini[nm];
        std::cout << ")\n";

        E_position[nm].print("Peak positions, eV");
        FCFs[nm].print("1D Harmonic Franck-Condon factors");
        I[nm].print("Intensities (FCFs^2)*(initial vibrational states termal population)");
      }

    } // end for each normal mode


    //================================================================================
    // evaluate the overal intensities as all possible products of 1D intensities (including combination bands)
    std::cout << "Maximum number of vibrational excitations: " << max_n_initial << " and "<< max_n_target 
      << "\n in the initial and each target state, respectively.\n\n";

    unsigned long total_combs_ini=0, total_combs_targ=0;
    for (int curr_max_ini=0; curr_max_ini<=max_n_initial; curr_max_ini++)
      total_combs_ini += Combination(  (curr_max_ini + nm_list.size() - 1), (nm_list.size() - 1)  );

    for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
      total_combs_targ += Combination( (curr_max_targ + nm_list.size() - 1), (nm_list.size() - 1) );

    std::cout << "Maximum number of combination bands = " << total_combs_ini * total_combs_targ  
      << "\n   = " << total_combs_ini << " (# of vibrational states in the initial electronic state)"
      << "\n   * " << total_combs_targ  << " (# of vibrational states in the target electronic state)\n\n" << std::flush;

    // This is just a reseting of these variables for a use with the next molecular state. 
    // TODO: It would be an easier read if definitions of these variabes were here. Paweł Apr '22
    for (int i=0; i<state_ini_subspace.size(); i++)
      state_ini_subspace[i]=-1;
    for (int i=0; i<state_targ_subspace.size(); i++)
      state_targ_subspace[i]=-1;

    // find INITIAL states with up to 'max_n_initial' vibrational quanta and with energy below 'energy_threshold_initial':
    std::cout << "A set of initial vibrational states is being created...\n";
    if (energy_threshold_initial < DBL_MAX)    
      std::cout << "  energy threshold = " << std::fixed << energy_threshold_initial <<" eV ("<< energy_threshold_initial/KELVINS2EV <<" K)\n" << std::flush;
    else
      std::cout << "  energy threshold is not specified in the input\n"<<std::flush;
    //      std::cout << "NOTE: ezSpectrum may crash at this point if memory is insufficient to store\n"
    //        	 <<"      all vibrational states. If so, please reduce the initial state's energy\n"
    //		 <<"      threshold or max_vibr_excitations_in_initial_el_state\n\n" << std::flush;

    // TODO: It's easier to read when a variable is defined where it is used. Paweł Apr '22
    selected_states_ini.clear();
    for (int curr_max_ini=0; curr_max_ini<=max_n_initial; curr_max_ini++)
    {
      // state_targ_subspace.size() is the number of normal modes in the 'excite subspace'
      while (enumerateVibrStates(state_targ_subspace.size(), curr_max_ini, state_ini_subspace, if_comb_bands))
      {
        // reset the index list state_ini
        for (int nm=0; nm < state_ini.size(); nm++)
          state_ini[nm]=0;

        // copy indexes from the subspace state_ini_subspace into the full space state_ini (the rest stays=0):
        int index=0;
        for (int nm=0; nm<nm_list.size(); nm++)
          state_ini[ nm_list[nm] ] = state_ini_subspace[nm];

        double energy = 0;
        for (int nm=0; nm < n_norm_modes; nm++) 
          energy += E_position[nm](state_ini[nm], 0) + molStates[targN].Energy(); 

        if (energy < energy_threshold_initial)
        {
          // TODO: Is the commented code still valuable? Paweł Apr '22
          // check memory available (dirty, but anyway one copying of the state is requared...)
          //int * buffer = (int*) malloc ( sizeof(int)*n_norm_modes*100 );
          //if (buffer==NULL)
          //  {
          //    std::cout << "\nError: not enough memory available to store all initial vibrational states\n\n";
          //    exit(2);
          //  }
          //free(buffer);

          selected_states_ini.push_back(state_ini);
        }

      }
      //reset the initial state's vibration "population"

      state_ini_subspace[0]=-1;
    }
    std::cout << "  " << selected_states_ini.size() << " vibrational states below the energy threshold\n\n"<<std::flush;


    // TODO: This is a duplicate of the above code. It should be a function. Paweł Apr '22
    // find TARGET states with up to 'max_n_target' vibrational quanta and with energy below 'energy_threshold_target':
    std::cout << "A set of target vibrational states is being created...\n";
    if (energy_threshold_target < DBL_MAX)
      std::cout << "  energy threshold = " << std::fixed <<energy_threshold_target <<" eV ("<< energy_threshold_target/WAVENUMBERS2EV <<" cm-1)\n"<<std::flush;
    else
      std::cout << "  energy threshold is not specified in the input\n"<<std::flush;
    //  std::cout << "NOTE: ezSpectrum may crash at this point if memory is insufficient to store\n"
    //   	 <<"      all vibrational states. If so, please reduce the target state's energy\n"
    //		 <<"      threshold or max_vibr_excitations_in_target_el_state\n\n" << std::flush;

    selected_states_targ.clear();
    for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
    {
      while (enumerateVibrStates(state_targ_subspace.size(), curr_max_targ, state_targ_subspace, if_comb_bands))
      {
        // reset the index list state_targ
        for (int nm=0; nm < state_targ.size(); nm++)
          state_targ[nm]=0;
        //copy indexes from the subspace state_targ_subspace into the full space state_targ (the rest stays=0):
        int index=0;
        for (int nm=0; nm<nm_list.size(); nm++)
          state_targ[ nm_list[nm] ] = state_targ_subspace[nm];


        double energy = 0;

        for (int nm=0; nm < n_norm_modes; nm++)
          // threshold -- energy above the ground state:
          energy += E_position[nm](0, state_targ[nm])+molStates[targN].Energy(); 

        if ( -energy < energy_threshold_target )
        {
          // check memory available (dirty, but anyway one copying of the state is requared...)
          //VibronicState * state_tmp = (VibronicState*) malloc ( sizeof(VibronicState)*2 );
          //if (state_tmp==NULL)
          //  {
          //    std::cout << "\nError: not enough memory available to store all target vibrational states\n\n";
          //    exit(2);
          //  }
          //free (state_tmp);

          selected_states_targ.push_back(state_targ);
        }


      }
      //reset the target state's vibration "population"
      state_targ_subspace[0]=-1;
    }

    std::cout << "  " << selected_states_targ.size() << " vibrational states below the energy threshold\n\n";

    std::cout << "Total number of combination bands with thresholds applied: " << selected_states_ini.size()*selected_states_targ.size() <<"\n\n"<<std::flush;

    // TODO: Here is a great place to add to the state the selected additional state that Marty asked for

    std::cout << "Intensities of combination bands are being calculated...\n" <<std::flush;

    for (int curr_ini=0; curr_ini<selected_states_ini.size(); curr_ini++)
      for (int curr_targ=0; curr_targ<selected_states_targ.size(); curr_targ++)
      {
        double intens = 1;
        double FCF = 1;
        double energy = -molStates[targN].Energy();
        double E_prime_prime = 0;

        for (int nm=0; nm < n_norm_modes; nm++)
        {
          intens *= I[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]);
          FCF *= FCFs[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]);
          energy += E_position[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]) + molStates[targN].Energy(); 
          E_prime_prime += E_position[nm](selected_states_ini[curr_ini][nm], 0) + molStates[targN].Energy(); 
          // [ cancell the IE in each energy, which is stupid but inexpensive; probably should be eliminated ]
        }

        // add the point to the spectrum if its intensity is above the threshold
        if (intens > intens_threshold)
        {
          tmpPoint.getIntensity() = intens;
          tmpPoint.getEnergy() = energy;
          tmpPoint.getE_prime_prime()= E_prime_prime;
          tmpPoint.getFCF()= FCF;
          tmpPoint.getVibrState1().reset();
          tmpPoint.getVibrState1().setElStateIndex(0);
          tmpPoint.getVibrState2().reset();
          tmpPoint.getVibrState2().setElStateIndex(targN);
          for (int nm=0; nm<nm_list.size(); nm++)
          {
            tmpPoint.getVibrState1().setVibrQuanta(nm, selected_states_ini[curr_ini][  nm_list[nm]  ]);
            tmpPoint.getVibrState2().setVibrQuanta(nm, selected_states_targ[curr_targ][  nm_list[nm]  ]);
          }
          spectrum.AddSpectralPoint( tmpPoint );
        }
      }
  } // end for each target state
}

// An exact copy of the previous constructor with some changes added to accomodate
// the Marty's request. 
// Prints transitions from a selected initial state.
// This constructor is activated by the keyword the_only_initial_state
Parallel::Parallel(std::vector <MolState>& molStates, std::vector<int> nm_active_space, 
    double fcf_threshold, double temperature, 
    std::vector<int> initial_state, int max_n_target, 
    bool if_comb_bands, bool if_use_target_nm, bool if_print_fcfs, 
    bool if_web_version, const char* nmoverlapFName,  // TODO: is the web version still active: Paweł Apr '22
    double energy_threshold_target)
{
  double intens_threshold=fcf_threshold*fcf_threshold;

  int N = nm_active_space.size();

  //initial state index
  int  iniN=0;

  // full space size
  int n_norm_modes = molStates[iniN].NNormModes();
  if (initial_state.size() != n_norm_modes) {
    std::cout << "Error in Parallel::Parallel: Invalid intial state:" << std::endl;
    std::cout 
      << " Excitations in all "  
      << n_norm_modes 
      << " must be defined."
      << std::endl;
    std::cout << " Only " 
      << initial_state.size() 
      << " input modes are given in the_only_initial_state."
      << std::endl;
    exit(2);
  }

  // Prohibit excitations of the modes that were disable in do_not_excite_subspace node
  auto end = nm_active_space.end();
  auto begin = nm_active_space.begin();
  for (int nm = 0; nm < N; nm++) {
    if (std::find(begin, end, nm) != end) {
      // mode is active so can be used
      continue;
    }
    else {
      // mode belongs to do_not_excite_subspace: must stay unexcited
      if (initial_state[nm] != 0) {
        std::cout << "Error in Parallel::Parallel." << std::endl;
        std::cout << " A mode from the do_not_excite_subspace #"
          << nm << " cannot be excited in the_only_initial_state."
          << std::endl;
        exit(1);
      }
    }
  }

  auto max_it = std::max_element(initial_state.begin(), initial_state.end());
  int max_n_initial = *max_it;

  std::vector <int> state_targ, state_targ_subspace;
  for (int i=0; i<n_norm_modes; i++)
    state_targ.push_back(-1);

  // "excite subspace" size
  for (int i=0; i<nm_active_space.size(); i++) 
    state_targ_subspace.push_back(-1);

  // stores set of the inital and target vibrational states (for each electronic state) -- after energy threshold applied
  std::vector < std::vector <int> > selected_states_ini, selected_states_targ;

  // matrix with FC factors (use the same variable for every target state)
  // max_n_initial stands for maximum allowed excitations in a single mode in the initial state
  // FCFs_tmp(n, m) gives the integral of an overlap of a state with 'n' excitation with 
  // a state of 'm' excitations, i.e, <n|m> 
  arma::Mat<double> FCFs_tmp(max_n_initial+1, max_n_target+1);
  // TODO: this might be more resource friendly if each FCFs in this vector had its own dims. Paweł Apr'22
  std::vector <arma::Mat<double>> FCFs;

  // matrix with intensities of FC transitions = FCFs*population_of_initial_vibrational_levels(Temperature distrib.)
  arma::Mat<double> I_tmp (max_n_initial+1, max_n_target+1);
  std::vector <arma::Mat<double>> I;

  // energy position of each transition for a given normal mode; 1D; ofset by IP; when add for N dimensions, substract IP from each energy;
  arma::Mat<double> E_position_tmp (max_n_initial+1, max_n_target+1);
  std::vector <arma::Mat<double>> E_position;

  // reduced mass for FCF calcualtions -- mass weighted coordinates
  // TODO: This is not used. Is it still necessary? Paweł Apr'22
  double reducedMass=1.0;

  // ==========================================================================================
  // Calculate transformation matrix "cartesian->normal mode" coordinates (for each state):
  std::vector <arma::Mat<double>> Cart2Normal_Rs;

  for (int state=0; state<molStates.size(); state++)
  {
    arma::Mat<double> NormModes( CARTDIM*(molStates[state].NAtoms()), molStates[state].NNormModes(), arma::fill::zeros ); // NOT mass weighted normal modes i.e. (L~)=T^(0.5)*L
    //Get L (normal modes in cartesian coordinates mass unweighted (in Angstroms) ):
    for (int j=0; j < molStates[state].NAtoms(); j++) 
      for (int i = 0; i < molStates[state].NNormModes(); i++) 
        for (int k=0; k < CARTDIM; k++)
          NormModes(j*CARTDIM+k, i) = molStates[state].getNormMode(i).getDisplacement()(j*CARTDIM+k);

    //Make sqrt(T)-matrix (diagonal matrix with sqrt(atomic masses) in cartesian coordinates):
    arma::Mat<double> SqrtT( CARTDIM*(molStates[state].NAtoms()), CARTDIM*(molStates[state].NAtoms()), arma::fill::zeros);
    for (int i=0; i<molStates[state].NAtoms(); i++)
      SqrtT(i*CARTDIM, i*CARTDIM) 
        = SqrtT(i*CARTDIM+1, i*CARTDIM+1)
        = SqrtT(i*CARTDIM+2, i*CARTDIM+2)
        = sqrt(molStates[state].getAtom(i).Mass());

    // Cart->NormalModes transformation matrix R=L^T*sqrt(T)
    // q' = q+d = q+R*(x-x') (for the parallel normal mode approximation)
    // units of d are Angstr*sqrt(amu)
    NormModes = NormModes.t(); 
    NormModes *= SqrtT;

    Cart2Normal_Rs.push_back(NormModes);
  }

  // Initial geometry in cartesian coordinates:
  arma::Col<double> InitialCartCoord(CARTDIM*(molStates[iniN].NAtoms()));
  for (int i=0; i<molStates[iniN].NAtoms(); i++)
    for (int k=0; k<CARTDIM; k++)
      InitialCartCoord[i*CARTDIM+k] = molStates[iniN].getAtom(i).Coord(k);
  InitialCartCoord *= ANGSTROM2AU; 

  // initialize spectrlPoint: create vector of VibrQuantNumbers for the initital and target state 
  // (keep only non-zero qunta, i.e. from the "excite subspace")
  // Also add subspace mask (the rest of nmodes do not have excitations, and  would not be printed in the spectrum)
  SpectralPoint tmpPoint;
  for (int nm=0; nm<nm_active_space.size(); nm++)
  {
    tmpPoint.getVibrState1().addVibrQuanta(0, nm_active_space[nm]);
    tmpPoint.getVibrState2().addVibrQuanta(0, nm_active_space[nm]);
  }

  // For each target state, for each normal mode: calculate shift, FCFs, and print spectrum
  for (int targN=1; targN < molStates.size(); targN++)
  {
    std::cout << "===== Target state #" << targN << " =====\n\n"<< std::flush;

    // clear the FC factors and related data
    FCFs.clear();
    I.clear();
    E_position.clear();

    //----------------------------------------------------------------------
    // calculate dQ:

    // Target state geometry in cartesian coordinates:
    arma::Col<double> TargetCartCoord(CARTDIM*(molStates[iniN].NAtoms()));
    for (int i=0; i<molStates[targN].NAtoms(); i++)
      for (int k=0; k <CARTDIM; k++)
        TargetCartCoord[i*CARTDIM+k] = molStates[targN].getAtom(i).Coord(k);

    // input is in angstroms
    TargetCartCoord *= ANGSTROM2AU;

    // Shift of the target state relative to the initial (In cart. coord.)
    arma::Col<double> cartesian_shift = TargetCartCoord - InitialCartCoord;

    // Transform Cart_>normal Mode coordinates

    // Geometry differences in normal coordinates
    arma::Col<double> NormModeShift_ini(molStates[iniN].NNormModes());
    arma::Col<double> NormModeShift_targ(molStates[iniN].NNormModes());

    NormModeShift_ini = Cart2Normal_Rs[iniN] * cartesian_shift;
    NormModeShift_targ = Cart2Normal_Rs[targN] * cartesian_shift;

    // Rescale it back & print:
    NormModeShift_ini *= AU2ANGSTROM;
    NormModeShift_targ *= AU2ANGSTROM;

    std::cout << "Difference (dQ) between the initial and the target state geometries.\n"
      << "Angstrom*sqrt(amu):\n\n"
      << "normal mode  dQ in initial  dQ in target   frequency   frequency   comments\n"
      << "  number      state coord.  state coord.    initial      target\n\n";
    for (int nm=0; nm<n_norm_modes; nm++)
    {
      std::cout << "     " << std::fixed << std::setw(3) << nm <<"      " 
        << std::setprecision(6) 
        << std::setw(9) <<  NormModeShift_ini[nm] << "       " 
        << std::setw(9) << NormModeShift_targ[nm] << "      " 
        << std::setprecision(2) 
        << std::setw(7) << molStates[iniN].getNormMode(nm).getFreq() << "    " 
        << std::setw(7) << molStates[targN].getNormMode(nm).getFreq() << "\n";
    }
    std::cout<<"\n\n";


    //----------------------------------------------------------------------
    // save the spectrumnm overlap to file (for normal mode reordering tool in the web interface)
    if (if_web_version)
    {
      std::vector <int> nondiagonal_list;
      arma::Mat<double> NMoverlap;
      bool if_overlap_diagonal=molStates[targN].getNormalModeOverlapWithOtherState(molStates[iniN], NMoverlap, nondiagonal_list);

      std::ofstream nmoverlapF;     
      nmoverlapF.open(nmoverlapFName, std::ios::out);
      nmoverlapF << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
        << "<nmoverlap\n  nm_initial=\"" <<n_norm_modes <<"\"\n  nm_target =\""<<n_norm_modes << "\">\n\n<nm_order_initial>\n";
      //inital state's nm order -- 0,1,2... -- always!
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<oi" << nm << ">"<< molStates[iniN].getNormModeIndex(nm) << "</oi" << nm << ">\n";
      nmoverlapF << "</nm_order_initial>\n\n<frequencies_initial>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<fi" << nm << ">"<< std::fixed << std::setprecision(1) << molStates[iniN].getNormMode(nm).getFreq()<< "</fi" << nm << ">\n";
      nmoverlapF << "</frequencies_initial>\n\n<displacements_initial>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<dqi" << nm << ">"<< std::fixed << std::setprecision(2) << NormModeShift_ini[nm]<< "</dqi" << nm << ">\n";
      nmoverlapF << "</displacements_initial>\n\n<nm_order_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<ot" << nm << ">"<< molStates[targN].getNormModeIndex(nm) << "</ot" << nm << ">\n";
      nmoverlapF << "</nm_order_target>\n\n<frequencies_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<ft" << nm << ">"<< std::fixed << std::setprecision(1) << molStates[targN].getNormMode(nm).getFreq()<< "</ft" << nm << ">\n";
      nmoverlapF << "</frequencies_target>\n\n<displacements_target>\n";
      for (int nm=0; nm<n_norm_modes; nm++)
        nmoverlapF << "<dqt" << nm << ">"<< std::fixed << std::setprecision(2) << NormModeShift_targ[nm]<< "</dqt" << nm << ">\n";
      nmoverlapF << "</displacements_target>\n\n<overlap_matrix>\n";
      for (int nmt=0; nmt<n_norm_modes; nmt++)
      {    
        nmoverlapF << "<row"<<nmt<<">\n";
        for (int nmi=0; nmi<n_norm_modes; nmi++)
          nmoverlapF << "<c"<<nmi<<">"<< NMoverlap(nmt, nmi)<<"</c"<<nmi<<">";
        nmoverlapF << "\n</row"<<nmt<<">\n";
      }
      nmoverlapF << "</overlap_matrix>\n\n</nmoverlap>";
      nmoverlapF.close();
    }
    //----------------------------------------------------------------------


    if (if_use_target_nm) 
      std::cout<<"TARGET state normal coordinates (displacements dQ) will be used.\n\n"; 
    else 
      std::cout<<"INITIAL state normal coordinates are used.\n\n"; 
    // ----------------------------------------------------------------------

    if (if_print_fcfs)
      std::cout << "Peak positions, Franck-Condon factors, and intensities along each normal mode:\n\n";

    // for each normal mode of the initial state
    for (int nm=0; nm < n_norm_modes; nm++) 
    {      
      // Calculate the energy positions:
      for (int i=0; i<max_n_initial+1; i++)
        for (int j=0; j<max_n_target+1; j++)
          // energy position of each transition for a given normal mode; 1D; ofset by IP; when add for N dimensions, substract IP from each energy;
          E_position_tmp(i, j)
            =
            - molStates[targN].Energy()
            + (molStates[iniN].getNormMode(nm).getFreq() * i
                - molStates[targN].getNormMode(nm).getFreq() * j) * WAVENUMBERS2EV;
      E_position.push_back(E_position_tmp);

      // Calculate the Franck-Condon factors
      if (if_use_target_nm) // which normal modes (displacements) to use: of the initial or the target state?
        harmonic_FCf(FCFs_tmp, reducedMass, NormModeShift_targ[nm],
            molStates[iniN].getNormMode(nm).getFreq(), 
            molStates[targN].getNormMode(nm).getFreq());
      else
        harmonic_FCf(FCFs_tmp, reducedMass, NormModeShift_ini[nm],
            molStates[iniN].getNormMode(nm).getFreq(), 
            molStates[targN].getNormMode(nm).getFreq());

      FCFs.push_back(FCFs_tmp);

      // Calculate the line intensities (FCF^2 *  the temperature distribution in the initial state):
      for (int i=0; i<max_n_initial+1; i++)
      { 
        double IExponent=100*i; //(intensity < 10e-44 for non ground states if temperature=0)
        if (temperature>0)
        {
          IExponent = molStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * i 
            / (temperature * KELVINS2EV);
          if (IExponent > 100) 
            IExponent=100; // keep the intensity >= 10e-44 == exp(-100)
        }
        // TODO: The intensity should neglected if one only considers fluorescence
        for (int j=0; j<max_n_target+1; j++)
          I_tmp(i,j) = FCFs_tmp(i,j) * FCFs_tmp(i,j) * exp(-IExponent);
      }

      I.push_back(I_tmp);

      if (if_print_fcfs)
      {
        std::cout << "NORMAL MODE #" << nm
          << " (" 
          << std::fixed << std::setw(8) << std::setprecision(2)
          << molStates[iniN].getNormMode(nm).getFreq() << "cm-1 --> " 
          << molStates[targN].getNormMode(nm).getFreq() << "cm-1, dQ="
          << std::setw(6) << std::setprecision(4);
        if (if_use_target_nm)
          std::cout << NormModeShift_targ[nm];
        else
          std::cout << NormModeShift_ini[nm];
        std::cout << ")\n";

        E_position[nm].print("Peak positions, eV");
        FCFs[nm].print("1D Harmonic Franck-Condon factors");
        I[nm].print("Intensities (FCFs^2)");
      }

    } // end for each normal mode


    //================================================================================
    // evaluate the overal intensities as all possible products of 1D intensities (including combination bands)
    std::cout 
      << "Maximum number of vibrational excitations in each target state: "
      << max_n_target 
      << ".\n\n";

    // HINT: Only one vibrational state in the initial vibrational state. Paweł Apr'22
    unsigned long total_combs_ini=1, total_combs_targ=0;
    for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
      total_combs_targ += Combination( (curr_max_targ + nm_active_space.size() - 1), (nm_active_space.size() - 1) );

    std::cout 
      << "Maximum number of combination bands = " 
      << total_combs_targ
      << "\n   = " 
      << " (# of vibrational states in the target electronic state)\n\n" 
      << std::flush;

    std::cout << "The SINGLE initial vibrational state is read from the input.\n";
    selected_states_ini.clear();
    selected_states_ini.push_back(initial_state);

    // find TARGET states with up to 'max_n_target' vibrational quanta and with energy below 'energy_threshold_target':
    std::cout << "A set of target vibrational states is being created...\n";
    // This is just a reseting of these variables for a use with the next molecular state. 
    for (int i=0; i<state_targ_subspace.size(); i++)
      state_targ_subspace[i]=-1;

    if (energy_threshold_target < DBL_MAX)
      std::cout << "  energy threshold = " << std::fixed <<energy_threshold_target <<" eV ("<< energy_threshold_target/WAVENUMBERS2EV <<" cm-1)\n"<<std::flush;
    else
      std::cout << "  energy threshold is not specified in the input\n"<<std::flush;
    //  std::cout << "NOTE: ezSpectrum may crash at this point if memory is insufficient to store\n"
    //   	 <<"      all vibrational states. If so, please reduce the target state's energy\n"
    //		 <<"      threshold or max_vibr_excitations_in_target_el_state\n\n" << std::flush;

    selected_states_targ.clear();
    for (int curr_max_targ=0; curr_max_targ<=max_n_target; curr_max_targ++)
    {
      while (enumerateVibrStates(state_targ_subspace.size(), curr_max_targ, state_targ_subspace, if_comb_bands))
      {
        // reset the index list state_targ
        for (int nm=0; nm < state_targ.size(); nm++)
          state_targ[nm]=0;
        //copy indexes from the subspace state_targ_subspace into the full space state_targ (the rest stays=0):
        int index=0;
        for (int nm=0; nm<nm_active_space.size(); nm++)
          state_targ[ nm_active_space[nm] ] = state_targ_subspace[nm];


        double energy = 0;

        for (int nm=0; nm < n_norm_modes; nm++)
          // threshold -- energy above the ground state:
          energy += E_position[nm](0, state_targ[nm])+molStates[targN].Energy(); 

        if ( -energy < energy_threshold_target )
        {
          // check memory available (dirty, but anyway one copying of the state is requared...)
          //VibronicState * state_tmp = (VibronicState*) malloc ( sizeof(VibronicState)*2 );
          //if (state_tmp==NULL)
          //  {
          //    std::cout << "\nError: not enough memory available to store all target vibrational states\n\n";
          //    exit(2);
          //  }
          //free (state_tmp);

          selected_states_targ.push_back(state_targ);
        }


      }
      //reset the target state's vibration "population"
      state_targ_subspace[0]=-1;
    }

    std::cout << "  " << selected_states_targ.size() << " vibrational states below the energy threshold\n\n";

    std::cout << "Total number of combination bands with thresholds applied: " << selected_states_ini.size()*selected_states_targ.size() <<"\n\n"<<std::flush;

    std::cout << "Intensities of combination bands are being calculated...\n" <<std::flush;

    for (int curr_ini=0; curr_ini<selected_states_ini.size(); curr_ini++)
      for (int curr_targ=0; curr_targ<selected_states_targ.size(); curr_targ++)
      {
        double intens = 1;
        double FCF = 1;
        double energy = -molStates[targN].Energy();
        double E_prime_prime = 0;

        for (int nm=0; nm < n_norm_modes; nm++)
        {
          intens *= I[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]);
          FCF *= FCFs[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]);
          energy += E_position[nm](selected_states_ini[curr_ini][nm], selected_states_targ[curr_targ][nm]) + molStates[targN].Energy(); 
          E_prime_prime += E_position[nm](selected_states_ini[curr_ini][nm], 0) + molStates[targN].Energy(); 
          // [ cancell the IE in each energy, which is stupid but inexpensive; probably should be eliminated ]
        }

        // add the point to the spectrum if its intensity is above the threshold
        if (intens > intens_threshold)
        {
          tmpPoint.getIntensity() = intens;
          tmpPoint.getEnergy() = energy;
          tmpPoint.getE_prime_prime()= E_prime_prime;
          tmpPoint.getFCF()= FCF;
          tmpPoint.getVibrState1().reset();
          tmpPoint.getVibrState1().setElStateIndex(0);
          tmpPoint.getVibrState2().reset();
          tmpPoint.getVibrState2().setElStateIndex(targN);
          for (int nm=0; nm<nm_active_space.size(); nm++)
          {
            tmpPoint.getVibrState1().setVibrQuanta(nm, selected_states_ini[curr_ini][  nm_active_space[nm]  ]);
            tmpPoint.getVibrState2().setVibrQuanta(nm, selected_states_targ[curr_targ][  nm_active_space[nm]  ]);
          }
          spectrum.AddSpectralPoint( tmpPoint );
        }
      }
  } // end for each target state
}
