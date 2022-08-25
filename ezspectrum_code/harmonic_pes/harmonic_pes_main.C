#include "harmonic_pes_main.h"

//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string& ex_str, int& qnt, int& nm );
//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state (i.e. vector of integers)
void fillVibrState(My_istringstream& vibr_str, VibronicState& v_state, const int nm_max);

void harmonic_pes_parallel(xml_node& node_input, std::vector <MolState>& elStates, const char *InputFileName);
void harmonic_pes_dushinksy(xml_node& node_input, std::vector <MolState>& elStates, const char *InputFileName);

bool harmonic_pes_main (const char *InputFileName, xml_node& node_input, xml_node& node_amu_table)
{
  //======= read "global" job variables  =====================================================
  xml_node node_jobparams(node_input,"job_parameters",0);
  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes=node_input.read_flag_value("print_normal_modes");

  //===========================================================================================
  //Read initial state and N target states; i.e. (N+1) electronic states total
  std::vector <MolState> elStates;

  xml_node node_istate(node_input,"initial_state",0);
  MolState elSt;
  std::cout << "\n====== Reading the initial state ======\n";
  elSt.Read(node_istate,node_amu_table);
  elStates.push_back(elSt);

  //std::cout << "\n====== Reading the initial state: Done ======\n";

  size_t n_target_states=node_input.find_subnode("target_state");

  for (int state_i=0; state_i<n_target_states; state_i++) {
    MolState elSt_t;
    xml_node node_t_state(node_input,"target_state",state_i);
    std::cout << "===== Reading the target state #" << state_i << " =====\n";
    elSt_t.Read(node_t_state,node_amu_table);
    elStates.push_back(elSt_t);
  }    
  std::cout << "Done reading states" << std::endl << std::endl;

  // Perform various checks and transformations
  if (elStates.size() <= 1) {
    std::cout << "\nError! No target states found in the input.\n\n";
    exit(2);
  }

  if (elStates[0].IfGradient()) {
    std::cout 
      << "\nError! Use of the vertical gradient method is allowed only in target states.\n\n";
    exit(2);
  }

  for (int state_i=0; state_i<elStates.size(); state_i++) {
    // ifSimilar checks:
    // - same number of atoms, 
    // - same order of the atomic names, 
    // - same "linearity"
    // Check for a consistent use of the vertical gradient method
    if (state_i>0) {
      if ( not(elStates[state_i].ifSimilar(elStates[0])) ) {
        std::cout << "Error: target state #" 
          << state_i 
          << " is different from the initial state\n\n";
        exit(2);
      }
      // ^ is a binary XOR: returns True if the two bools are different
      bool gradient_used_in_the_first_target = elStates[1].IfGradient();
      if (elStates[state_i].IfGradient() ^ gradient_used_in_the_first_target) {
        std::cout << "Error: target state #" 
          << state_i 
          << " breaks the consistent use of the VG method in target states.\n\n"
          << " Make sure to use vertical gradient in all target states.\n\n";
        exit(2);
      }
    }

    // apply automatic transformations to the last loaded state (if no manual were requested):
    if ( not(elStates[state_i].ifAlignedManually()) ) {

      // align each state: 
      // 1. center of mass in the coordinates origin
      // 2. moment of ineretia principal axes along the coordinate axes
      elStates[state_i].align();

      // align each target state with the initial one
      if (state_i>0)
        elStates[state_i].align(elStates[0]);

      // get clean zeros
      elStates[state_i].applyCoordinateThreshold(COORDINATE_THRESHOLD);
    }

    std::cout << "\nNew molecular geometry:\n";
    elStates[state_i].printGeometry();
    // centerOfMass=elStates[i].getCenterOfMass().applyThreshold(COORDINATE_THRESHOLD).print("Center of mass: ");
    elStates[state_i].getMomentOfInertiaTensor().print("\nMOI tensor:");

    if (if_print_normal_modes) {
      std::cout << "Normal modes after the geometry transformations:\n\n";
      elStates[state_i].printNormalModes();
    }
  } 

  std::cout << "Done with the transformations" << std::endl;
  std::cout << "------------------------------------------------------------------------------\n";

  // total number of normal modes (in the initial state)
  int n_norm_modes = elStates[0].NNormModes();

  // if parallel or dushinsky
  bool if_something_to_do=false;

  //======================================================================
  // Parallel approximation section

  if ( node_input.find_subnode("parallel_approximation")) {
    if_something_to_do=true;
    harmonic_pes_parallel(node_input,elStates,InputFileName);
  }

  //=========================================================================
  // Dushinski rotation (nonparallel approximation)
  //=========================================================================
  // the notation and equations are from [Berger et al. JPCA 102:7157(1998)]
  //

  if(node_input.find_subnode("dushinsky_rotations")) {
    if_something_to_do=true;
    harmonic_pes_dushinksy(node_input,elStates,InputFileName);
  }

  if (!if_something_to_do) {                 
    std::cout << "\nError! No \"parallel_approximation\" or \"dushinsky_rotations\" section\n       was found in the input. Nothing to do.\n\n";
    exit(2);
  }

  return true;
};


//! splits string of type "3v21" into two integers 3 and 21
void get_qnt_nm(std::string& ex_str, int& qnt, int& nm ) {
  if (ex_str.find("v")==std::string::npos) {
    std::cout << "\nFormat error in [" << ex_str << "] excitation: should contain symbol \'v\'\n\n";
    exit(1);
  }
  ex_str.replace( ex_str.find("v"), 1, " " );
  std::istringstream ex_strs(ex_str);
  ex_strs>>qnt>>nm;

  if (ex_strs.fail()) {
    ex_str.replace( ex_str.find(" "), 1, "v" );
    std::cout << "\nFormat error in [" << ex_str << "] excitation. Should be two integers separated by the symbol 'v' \n\n";
    exit(1);
  }
}


//! converts string of type "1v1,1v2,1v3,3v19" into a vibrational state (i.e. vector of integers)
void fillVibrState(My_istringstream& vibr_str, VibronicState& v_state, const int nm_max) {

  // quanta & normal mode number (for parsing strings like "3v19", where qnt=3 and nm=19)
  int qnt=0, nm=0;
  // string like 3v19"
  std::string ex_str;
  // string like "1v1,1v2,1v3,3v19"
  //std::cout << "fillVibrState: vibr_str= " << vibr_str.str() << std::endl;
  bool if_read=vibr_str.getNextWord(ex_str);
  /* std::cout << "fillVibrState: if_read=" << if_read << "  ex_str= " << ex_str << std::endl; */

  // reset vibrational state
  for (int i=0; i<v_state.getVibrQuantaSize(); i++)
    v_state.setVibrQuanta(i,0);

  // fill vibrational state (if == 0 -- nothing to do)
  if (ex_str!="0") {
    get_qnt_nm(ex_str, qnt, nm);
    /* std::cout << "The word " << ex_str << " was read as qnt = " << qnt << " nm = " << nm << std::endl; */

    /* std::cout << "nm_max = " << nm_max << std::endl; */

    if (nm>nm_max) {
      std::cout << "\nError! Normal mode " << nm 
        << " (in [" << qnt << 'v' << nm << "] excitation) is out of range.\n\n";
      exit(1);
    }

    v_state.setVibrQuanta(nm, qnt);
    while (not(vibr_str.fail())) {
      vibr_str.getNextWord(ex_str);
      get_qnt_nm(ex_str, qnt, nm);
      v_state.setVibrQuanta(nm,qnt);
    }

  }
}


void harmonic_pes_parallel(xml_node& node_input, std::vector <MolState>& elStates, const char *InputFileName) {

  xml_node node_parallel_approx(node_input,"parallel_approximation",0);
  xml_node node_jobparams(node_input,"job_parameters",0);

  // read global paramters
  double temperature=node_jobparams.read_double_value("temperature");
  // fcf threshold (from the <job_parameters> tag)
  double fcf_threshold=sqrt(node_jobparams.read_double_value("spectrum_intensity_threshold"));
  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes=node_input.read_flag_value("print_normal_modes");
  // check if the web version format of the output (do not print the input file & create a ".nmoverlap" file)
  bool if_web_version=node_input.read_flag_value("if_web_version");

  bool ifAnyNormalModesReordered=false;
  for (int state_i=0; state_i<elStates.size(); state_i++) {
    if ( elStates[state_i].ifNMReorderedManually() ) {
      ifAnyNormalModesReordered=true;
      if (state_i==0) {
        std::cout<<"\nError! Manual reordering of the normal modes is not allowed for the initial state\n\n";
        exit(2);
      }
    }
  }

  // total number of the normal modes (in the initial state)
  int n_norm_modes = elStates[0].NNormModes();

  std::cout << "\n=== Reading the parallel approximation job parameters ===\n"<< std::flush;

  // Maximum number of vibrational levels to take:
  int max_n_initial, max_n_target;  //maximum number of vibrational levels for initial and target state
  max_n_initial = node_parallel_approx.read_int_value("max_vibr_excitations_in_initial_el_state"); // i.e. =2 <=> three vibr. states GS and two excited states
  max_n_target  = node_parallel_approx.read_int_value("max_vibr_excitations_in_target_el_state");
  bool if_comb_bands = node_parallel_approx.read_bool_value("combination_bands");
  bool if_use_target_nm = node_parallel_approx.read_bool_value("use_normal_coordinates_of_target_states");

  if(temperature==0) {
    max_n_initial = 0 ;
    std::cout << "\nSince temperature=0, \"max_vibr_excitations_in_initial_el_state\" has been set to 0.\n"<< std::flush;
  }

  // check if print normal modes after transformations & overlap matrix
  //FIXIT: check if loc is correct
  // TODO: Does this still require any work? Pawel, Feb '22
  bool if_print_fcfs= node_parallel_approx.read_flag_value("print_franck_condon_matrices");

  // read energy thresholds (if provided)
  double energy_threshold_initial = DBL_MAX;//eV
  double energy_threshold_target = DBL_MAX; //eV
                                            // TODO: the next two ifs should be two functions;
  if(node_parallel_approx.find_subnode("energy_thresholds")) {

    xml_node node_energy_thresholds(node_parallel_approx,"energy_thresholds",0);

    if( node_energy_thresholds.find_subnode("initial_state")) {

      std::cout << "Reading energy thresholds. " <<  std::endl;
      xml_node node_istate(node_energy_thresholds,"initial_state",0);

      std::string units=node_istate.read_string_value("units");
      energy_threshold_initial=node_istate.read_node_double_value();
      //std::cout << "Thresh=" << energy_threshold_initial << " " << units << std::endl;
      if ( !covert_energy_to_eV(energy_threshold_initial,units) ) {
        std::cout << "\nError! Unknown units of the initial state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
        exit(1);
      }
    }

    if( node_energy_thresholds.find_subnode("target_state")) {

      xml_node node_tstate(node_energy_thresholds,"target_state",0);

      std::string units=node_tstate.read_string_value("units");
      energy_threshold_target=node_tstate.read_node_double_value();
      //std::cout << "Thresh=" << energy_threshold_target << " " << units << std::endl;
      if ( !covert_energy_to_eV(energy_threshold_target,units) ) {
        std::cout << "\nError! Unknown units of the target state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
        exit(1);
      }
    }
  }

  // read normal modes do_not_excite subspace (for parallel approximation only):
  std::set<int> do_not_excite_subspace;
  bool if_use_do_not_excite_subspace=false;
  int do_not_excite_subspace_size=0;
  int do_not_excite_subspace_max=0; // maximum value in lists -- error check later (dirty)

  if(node_parallel_approx.find_subnode("do_not_excite_subspace")) {

    xml_node node_do_not_excite_subspace(node_parallel_approx,"do_not_excite_subspace",0);

    if_use_do_not_excite_subspace=true;
    do_not_excite_subspace_size=node_do_not_excite_subspace.read_int_value("size");
    std::istringstream tmp_iStr(node_do_not_excite_subspace.read_string_value("normal_modes"));

    int tmpInt;

    for (int nm=0; nm<do_not_excite_subspace_size; nm++)  {
      tmp_iStr >> tmpInt;
      //input error check:
      if (tmp_iStr.fail()) {
        std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
          << "(non numeric symbol or less entries then specified by the \"size\" value)\n\n";
        exit(1);
      }
      if (tmpInt<0) {
        std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
          << "Entry ["<< tmpInt<<"] is negative.\n\n";
        exit(1);
      }
      // keep the maximum value of the list
      if (tmpInt>do_not_excite_subspace_max)
        do_not_excite_subspace_max=tmpInt;

      //check if tmpInt is already in the set:
      std::set<int>::const_iterator intSet_iter;
      intSet_iter = do_not_excite_subspace.find(tmpInt);
      if ( intSet_iter != do_not_excite_subspace.end() ) {
        std::cout << "\nFormat error in \"input\"->\"do_not_excite_subspace\"->\"normal_modes\"\n"
          << "Entry ["<< tmpInt<<"] is not unique.\n\n";
        exit(1);
      }
      do_not_excite_subspace.insert(tmpInt);
    }

    if (do_not_excite_subspace.size()!=0) {

      if(ifAnyNormalModesReordered)
        std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
          <<"         New order is used for the \"do_not_excite_subspace\".\n\n";

      std::cout << "The following normal modes will have no vibrational excitations:\n";

      for (std::set<int>::const_iterator intSet_iter=do_not_excite_subspace.begin(); intSet_iter!=do_not_excite_subspace.end(); intSet_iter++)
        std::cout << *intSet_iter << ' ';
      std::cout<<"\n";
    }
  }

  //check that numbers in do_not_excite_subspace are less than the number_of_normal_modes
  if ((do_not_excite_subspace_max >= n_norm_modes)and(if_use_do_not_excite_subspace==true)) {
    std::cout << "\nError! Maximum normal mode number in \"do_not_excite_subspace\" is ["<< do_not_excite_subspace_max<<"],\n"
      << "  which is greater than (number_of_normal_modes-1)="<< n_norm_modes-1 <<"\n\n";
    exit(2);
  }

  //create nms_parallel -- "excite subspace" (full_space-do_not_excite_subspace)
  std::vector<int> nms_parallel;
  for (int nm=0; nm<n_norm_modes; nm++) {
    std::set<int>::const_iterator intSet_iter;
    intSet_iter = do_not_excite_subspace.find(nm);
    if ( intSet_iter == do_not_excite_subspace.end( ) )
      nms_parallel.push_back(nm);
  }


  //================================================================================
  // print the overlap matrix with the initial state for each target states:

  for (int state=1; state<elStates.size(); state++) {
    std::cout << "\n===== Overlap matrix of the target state #" << state << " with the initial state =====\n";

    std::vector <int> normal_modes_list;
    arma::Mat<double> NMoverlap;  //normal modes overlap matrix (for each target state the same matrix is used)
    bool if_overlap_diagonal;

    // select nondiagonal submatrix of the overlap matrix:
    if_overlap_diagonal=elStates[state].getNormalModeOverlapWithOtherState(elStates[0], NMoverlap, normal_modes_list);
    // rows -- norm modes of the target state; colums norm modes of the initial state;

    // remove normal modes from normal_modes_list that are in the do_not_excite_subspace:
    std::set<int>::iterator iter_set;
    std::vector<int> new_normal_modes_list; 

    for (int nm=0; nm<normal_modes_list.size(); nm++) {
      // if nm is not in the do_not_excite set:
      iter_set = do_not_excite_subspace.find(normal_modes_list[nm]);
      if ( iter_set == do_not_excite_subspace.end( ) )
        // then copy it to the new list:
        new_normal_modes_list.push_back(normal_modes_list[nm]);
    }

    //Create an overlap submatrix:
    if ((if_overlap_diagonal) or (new_normal_modes_list.size()<=1)) {
      std::cout << "The normal modes overlap matrix with the initial state is diagonal\n";
      if (new_normal_modes_list.size()<=1)
        std::cout<<"  (do_not_excite_subspace is excluded)\n";
      std::cout<<"\n";
    }
    else {        
      std::cout << "WARNING! The normal modes overlap matrix with the initial state\n"
        << "         is non-diagonal! Consider reordering the normal modes.\n\n";
      // create a normal mode submatrix:
      arma::Mat<double> overlap_submatrix(new_normal_modes_list.size(), new_normal_modes_list.size(), arma::fill::zeros);
      for (int nm1=0; nm1<new_normal_modes_list.size(); nm1++)
        for (int nm2=0; nm2<new_normal_modes_list.size(); nm2++)
          overlap_submatrix(nm1, nm2) = NMoverlap(new_normal_modes_list[nm1],new_normal_modes_list[nm2]);

      //print the overlap_submatrix (with correct column/row labbels):
      std::cout << "  The non-diagonal part of the normal modes overlap matrix (do_not_excite_subspace is excluded):";
      std::cout << "\n     ";

      for (int j=0; j<new_normal_modes_list.size(); j++)
        std::cout << std::fixed << std::setprecision(0) << std::setw(8) << new_normal_modes_list[j];
      for (int i=0; i<new_normal_modes_list.size(); i++) {
        std::cout << "\n  "<< std::fixed << std::setprecision(0) << std::setw(3) << new_normal_modes_list[i];
        for (int j=0; j<new_normal_modes_list.size(); j++)
          if (fabs(overlap_submatrix(i, j)) >= 0.001)
            std::cout << std::fixed << std::setprecision(3) << std::setw(8) << overlap_submatrix(i, j);
          else
            std::cout << "      --";
      }
      std::cout <<"\n\n";
    }

    // print in a "fit 80 chars wide terminal" form
    if(if_print_normal_modes)
      NMoverlap.print("Normal modes overlap matrix with the initial state \n(if significantly non diagonal, please consider normal modes reordering)");
  }

  // for the web version: save the overlap matrix (with displacements) in an xml file
  std::stringstream nmoverlapFName; 
  nmoverlapFName << InputFileName << ".nmoverlap";


  std::cout << "------------------------------------------------------------------------------\n\n";
  std::cout << "Photoelectron spectrum in the parallel approximation will be evaluated\n\n"<< std::flush;

  //================================================================================
  //================================================================================
  //================================================================================
  //================================================================================

  // create a new parallel approximation object (evaluates and stores FCFs in the harmonic approximation)
  Parallel* parallel_ptr;

  if (node_parallel_approx.find_subnode("the_only_initial_state")) {
    xml_node node_the_only_initial_state(node_parallel_approx, "the_only_initial_state", 0);
    std::string text = node_the_only_initial_state.read_string_value("text");
    if (text.empty()) {
      std::cout << "Error in processing the_only_initial_state." 
        << std::endl
        << " The selected state must have at least one excitation." 
        << std::endl;
      exit(1);
    }
    std::queue<std::string> non_zero_modes;

    std::string non_zero_mode;
    size_t pos = text.find(",");
    while (true) {
      non_zero_mode = text.substr(0, pos);
      non_zero_modes.push(non_zero_mode);
      if (pos == std::string::npos) 
        break;
      text.erase(0, pos + 1);
      pos = text.find(",");
    } 

    // container for the only initial state
    // initialized as a state which has every mode with zero excitations 
    std::vector<int> the_only_initial_state(n_norm_modes, 0); 

    while (non_zero_modes.size()) {
      std::string excitation = non_zero_modes.front();
      non_zero_modes.pop();
      int no_of_quanta = -1;
      int mode_number = -1;
      get_qnt_nm(excitation, no_of_quanta, mode_number);

      if (mode_number >= n_norm_modes || mode_number < 0) {
        std::cout << "Error in processing the_only_initial_state." 
          << std::endl
          << " Please pick normal mode from the range 0 to " 
          << n_norm_modes - 1 
          << "." 
          << std::endl;
        exit(1);
      }
      if (no_of_quanta < 0) {
        std::cout << "Error in processing the_only_initial_state." 
          << std::endl
          << " Please pick a positive number of vibrational quanta in all excited modes."
          << std::endl;
        exit(1);
      }
      if (the_only_initial_state[mode_number] != 0) {
        std::cout << "Error in processing the_only_initial_state." 
          << std::endl
          << " The number of vibrational quanta in mode #"
          << mode_number
          << " is defined more than once."
          << std::endl
          << " Please define each normal mode excitations only once."
          << std::endl;
        exit(1);
      }
      the_only_initial_state[mode_number] = no_of_quanta;
    }

    parallel_ptr = new Parallel(elStates, nms_parallel, 
        fcf_threshold, temperature, 
        the_only_initial_state, max_n_target, 
        if_comb_bands, if_use_target_nm, if_print_fcfs, if_web_version, 
        nmoverlapFName.str().c_str(),
        energy_threshold_target);
  }
  else {
    parallel_ptr = new Parallel(elStates, nms_parallel, 
        fcf_threshold, temperature, 
        max_n_initial, max_n_target, 
        if_comb_bands, if_use_target_nm, if_print_fcfs, if_web_version,
        nmoverlapFName.str().c_str(),  
        energy_threshold_initial,  energy_threshold_target);
  }

  //================================================================================
  //================================================================================
  //================================================================================
  //================================================================================


  //--------------------------------------------------------------------------------
  // Print the updated spectrum:
  (*parallel_ptr).getSpectrum().Sort();
  std::cout << "------------------------------------------------------------------------------\n";
  std::cout << "           Stick photoelectron spectrum (parallel approximation)\n";
  std::cout << "------------------------------------------------------------------------------\n";
  if (ifAnyNormalModesReordered)
    std::cout <<"\nWARNING! The normal modes of one of the target states were reordered!\n"
      <<"         New order is used for the target state assignment.\n";
  if( nms_parallel.size()!=n_norm_modes )
  {
    std::cout << "\nNOTE: only the following normal modes were excited: (\"excite subspace\"):\n  ";
    for (int nm=0; nm< nms_parallel.size(); nm++)
      std::cout << nms_parallel[nm] << ' ';
    std::cout << "\n";
    if (ifAnyNormalModesReordered)
      std::cout <<"\nWARNING! The normal modes of one of the target states were reordered!\n"
        <<"         New order is used for the \"excite subspace\"\n";

  }
  std::cout << "\n";
  if ( (*parallel_ptr).getSpectrum().getNSpectralPoints()>0)
    (*parallel_ptr).getSpectrum().PrintStickTable();
  else                
    std::cout << "\n\n\n"
      <<"WARNING! The spectrum is empty.\n\n"
      <<"         Plese refer to \"My spectrum is empty!\" in the\n"
      <<"         \"Common problems\" section of the manual\n\n\n\n";

  std::cout << "------------------------------------------------------------------------------\n";

  // save this spectrum to the file
  std::stringstream spectrumFName; 
  spectrumFName << InputFileName << ".spectrum_parallel";
  (*parallel_ptr).getSpectrum().PrintStickTable(spectrumFName.str().c_str());
  std::cout << "\nStick spectrum was also saved in \"" << spectrumFName.str() << "\" file \n";
  if(if_use_do_not_excite_subspace)
    std::cout << " (Full list of the normal modes was used for assigning transitions)\n";

  std::cout << "------------------------------------------------------------------------------\n\n";

  delete parallel_ptr;
}


void harmonic_pes_dushinksy(xml_node& node_input, std::vector <MolState>& elStates, const char *InputFileName) {

  xml_node node_dushinsky_rotations(node_input,"dushinsky_rotations",0);
  xml_node node_jobparams(node_input,"job_parameters",0);


  // TODO: this piece of code is repeated in the harmonic_pes_parallel. Paweł June '22

  // read global paramters
  double temperature=node_jobparams.read_double_value("temperature");
  // fcf threshold (from the <job_parameters> tag)
  double fcf_threshold=sqrt(node_jobparams.read_double_value("spectrum_intensity_threshold"));
  // check if print normal modes after transformations & overlap matrix
  bool if_print_normal_modes=node_input.read_flag_value("print_normal_modes");
  // check if the web version format of the output (do not print the input file & create a ".nmoverlap" file)
  bool if_web_version=node_input.read_flag_value("if_web_version");

  bool ifAnyNormalModesReordered=false;

  for (int state_i=0; state_i<elStates.size(); state_i++) {	

    if ( elStates[state_i].ifNMReorderedManually() ) {

      ifAnyNormalModesReordered=true;
      if (state_i==0) {
        std::cout<<"\nError! Manual reordering of the normal modes is not allowed for the initial state\n\n";
        exit(2);
      }
    }
  }
  // total number of the normal modes (in the initial state)
  int n_norm_modes = elStates[0].NNormModes();
  //Done, can proceed to do Dushinksy calculations

  std::cout << "\n\n=== Reading the Dushinsky rotations job parameters ===\n\n"<< std::flush;
  //----------------------------------------------------------------------
  // load parameters 

  // HINT: target_state is exclusively used in Duschinsky node. Paweł June '22
  // indexes of the initial and target electronic states:
  int iniN=0;
  int targN=node_dushinsky_rotations.read_int_value("target_state"); 

  if ((targN>elStates.size()-1)or(targN<1)) {
    std::cout<<"\nFormat error: \"target_state\" value should be positive and not greater than "<<elStates.size()-1<< "\n\n";
    exit(1);
  }
  std::cout<<"Target state number "<<targN<<" from the input will be used\n\n";

  // TODO: This is a continuation of the copy of code from harmonic_pes_parallel. Paweł June '22
  // maximum number of quanta to store:     
  int max_quanta_ini = node_dushinsky_rotations.read_int_value("max_vibr_excitations_in_initial_el_state");
  int max_quanta_targ = node_dushinsky_rotations.read_int_value("max_vibr_excitations_in_target_el_state");

  if(temperature==0) {
    max_quanta_ini = 0 ;
    std::cout << "\nSince temperature=0, \"max_vibr_excitations_in_initial_el_state\" has been set to 0.\n"<< std::flush;
  }

  // HINT: Back to the Duschinsky-exclusive part of the code
  int Kp_max_to_save=32000;
  if(node_dushinsky_rotations.find_subnode("max_vibr_to_store")) {

    xml_node node_max_vibr_to_store(node_dushinsky_rotations,"max_vibr_to_store",0);
    Kp_max_to_save=node_max_vibr_to_store.read_int_value("target_el_state");
  }

  // TODO: do_not_excite_subspace is also processed in the harmonic_pes_parallel 
  // but the code below is slightly different . Paweł June '22
  // read "do not excite subspace"
  std::set<int> do_not_excite_subspace;

  if (node_dushinsky_rotations.find_subnode("do_not_excite_subspace")) {

    xml_node node_do_not_excite_subspace(node_dushinsky_rotations,"do_not_excite_subspace",0);
    int ss_size=node_do_not_excite_subspace.read_int_value("size");
    if ((ss_size<0) or (ss_size>n_norm_modes)) {
      std::cout << "\nError: subspace size is out of range\n\n";
      exit(1);
    }
    if (ss_size>0) {

      if (elStates[targN].ifNMReorderedManually())
        std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
          <<"         New order is used for the \"do_not_excite_subspace\".\n\n";

      std::istringstream tmp_iStr(node_do_not_excite_subspace.read_string_value("normal_modes"));
      int tmpInt;

      for (int nm=0; nm<ss_size; nm++)  {
        tmp_iStr >> tmpInt;
        //input error check:
        if (tmp_iStr.fail()) {
          std::cout << "\nError: non-numeric symbol or fewer entries than specified by the \"size\" value\n\n";
          exit(1);
        }
        if ((tmpInt<0)or(tmpInt>n_norm_modes-1)) {
          std::cout << "\nError: Entry ["<< tmpInt<<"] is out of range.\n\n";
          exit(1);
        }

        //check if tmpInt is already in the set:
        std::set<int>::const_iterator intSet_iter;
        intSet_iter = do_not_excite_subspace.find(tmpInt);
        if ( intSet_iter != do_not_excite_subspace.end( ) ) {
          std::cout << "\nEntry ["<< tmpInt<<"] is not unique.\n\n";
          exit(1);
        }
        do_not_excite_subspace.insert(tmpInt);
      }
    }
  }

  //----------------------------------------------------------------------
  // create the "excite subspace"; only normal modes from this subspace will be excited
  std::vector<int> nms_dushinsky;
  for (int nm=0; nm<n_norm_modes;nm++) {
    std::set<int>::const_iterator intSet_iter;
    intSet_iter = do_not_excite_subspace.find(nm);
    if ( intSet_iter == do_not_excite_subspace.end() )
      nms_dushinsky.push_back(nm);
  }

  // HINT: This if section is absent in the harmonic_pes_parallel but wouldn't be bad to
  // have it there as well. Paweł June '22
  if (nms_dushinsky.size()<n_norm_modes) {
    std::cout << "The following normal will be excited\n(for both states the order is same as in the input):\n ";
    for (int nm=0; nm<nms_dushinsky.size(); nm++)
      std::cout << nms_dushinsky[nm] << ' ';
    std::cout<<"\n\n";
  }
  else
    std::cout << "All normal modes modes will be excited\n\n";

  // HINT: in the harmonic_pes_parallel now comes a section which prints the normal modes overalp. Paweł June '22

  std::cout << "=== Photoelectron spectrum with the Dushinsky rotations will be evaluated ==== \n\n"<< std::flush;

  //================================================================================
  //================================================================================
  //================================================================================
  //================================================================================

  // create a new dushinsky object for a given set of normal modes
  // All matrices and zero-zero integral are evaluated for the full space;
  // then excitations are only added to the normal modes from "nms_dushinsky"
  Dushinsky* dushinsky_ptr = new Dushinsky(elStates, nms_dushinsky, fcf_threshold, targN, max_quanta_targ, max_quanta_ini);

  //================================================================================
  //================================================================================
  //================================================================================
  //================================================================================

  // print estimated size of each layer up to K'_max:
  std::cout << "Number of normal modes to excite: " << nms_dushinsky.size() << "\n\n";
  std::cout << "Size of layers with exactly K' excitations in the target state (in bytes):\n";
  (*dushinsky_ptr).printLayersSizes( (max_quanta_targ<Kp_max_to_save)? max_quanta_targ : Kp_max_to_save );
  std::cout << "\n";
  //--------------------------------------------------------------------------------

  // go over all layers and add points to the spectrum:
  for (int Kp=1; Kp<=max_quanta_targ; Kp++)
  {
    if ( Kp<=Kp_max_to_save ) {
      std:: cout << "Layer K'="<< Kp <<" is being evaluated (will be saved in memory)... " << std::flush;
    }
    else {
      std:: cout << "Layer K'="<< Kp <<" is being evaluated... " << std::flush;
    }

    int n_fresh_points=(*dushinsky_ptr).evalNextLayer( Kp<=Kp_max_to_save );

    std:: cout << "Done\n";
    if (n_fresh_points>0) {
      std:: cout << n_fresh_points <<" points above the intensity threhold were added to the spectrum\n\n" << std::flush;
    }
    else {
      std:: cout << "No points above the intensity threhold were found in this layer\n\n" << std::flush;
    }
  }

  //----------------------------------------------------------------------
  // add the hot bands if requested

  double energy_threshold_initial = DBL_MAX;//eV
  double energy_threshold_target = DBL_MAX; //eV
  if(max_quanta_ini!=0) {

    if (node_dushinsky_rotations.find_subnode("energy_thresholds")) {

      xml_node node_energy_thresholds(node_dushinsky_rotations,"energy_thresholds",0);
      // read the energy thresholds (if provided)

      if ( node_energy_thresholds.find_subnode("initial_state")) {
        xml_node node_istate(node_energy_thresholds,"initial_state",0);
        std::string units=node_istate.read_string_value("units");
        energy_threshold_initial=node_istate.read_node_double_value();
        if ( !covert_energy_to_eV(energy_threshold_initial,units) ) {
          std::cout << "\nError! Unknown units of the initial state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
          exit(1);
        }

      }

      if ( node_energy_thresholds.find_subnode("target_state")) {

        xml_node node_tstate(node_energy_thresholds,"target_state",0);
        std::string units=node_tstate.read_string_value("units");
        energy_threshold_target=node_tstate.read_node_double_value();
        if ( !covert_energy_to_eV(energy_threshold_target,units) ) {
          std::cout << "\nError! Unknown units of the target state threshold: \"" << units <<"\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
          exit(1);
        }
      }
    }

    std::cout << "T=" << temperature << " FCF thresh=" << fcf_threshold << std::endl;
    std::cout << "Max quanta ini=" << max_quanta_ini <<   "  Max quanta targ=" << max_quanta_targ << std::endl;
    std::cout <<  "Thresh[ini]=" << energy_threshold_initial <<  "  Thresh[targ]="  << energy_threshold_target << std::endl;

    int n_hot_bands=(*dushinsky_ptr).addHotBands(elStates, nms_dushinsky, fcf_threshold, temperature, 
        max_quanta_ini, max_quanta_targ, 
        energy_threshold_initial,energy_threshold_target);
    std:: cout << n_hot_bands <<" hot bands were added to the spectrum\n"
      <<"Note: the Boltzmann distribution will be applied later\n\n" << std::flush;
  }

  //----------------------------------------------------------------------
  // If the_only_initial_state node is present in the xml input.
  // Hide the already calculated spectrum. Add transitions from the_only_initial_state. 
  // No finesse. 
  // the_only_initial_state is in the "full space"; 
  // do_not_excite_subspace does not apply;

  size_t do_the_only_initial_state = node_dushinsky_rotations.find_subnode("the_only_initial_state");
  if (do_the_only_initial_state) {
    if (do_the_only_initial_state != 1) {
      std::cout << "Error! Only one the_only_initial_state node is allowed." << std::endl;
      exit(1);
    }

    // Don't show all the other points as this is expected to be the_only_initial_state
    for (int pt=0; pt<(*dushinsky_ptr).getSpectrum().getNSpectralPoints(); pt++) {
      (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).setIfPrint(false);
    }

    std::cout << "Clearing the spectrum is ready" << std::endl;

    xml_node node_the_only_initial_state(node_dushinsky_rotations, "the_only_initial_state", 0);
    std::string text = node_the_only_initial_state.read_string_value("text");
    My_istringstream ini_str(text);
    VibronicState the_only_initial_state;
    for (int nm = 0; nm < n_norm_modes; nm++) {
      the_only_initial_state.addVibrQuanta(0, nm);
    }
    fillVibrState(ini_str, the_only_initial_state, n_norm_modes);

    std::cout << "The initial state is ready." << std::endl;

    // go over all layers and add points to the spectrum:
    dushinsky_ptr->reset_Kp_max();
    for (int Kp=0; Kp<=max_quanta_targ; Kp++)
    {
      dushinsky_ptr->add_the_only_intial_state_transitions(Kp, the_only_initial_state);
    }
  }

  //----------------------------------------------------------------------
  // now load the list of single transitions to evaluate FCFs recursively
  // single_transition is in the "full space"; do_not_excite_subspace does not apply;
  SpectralPoint single_transition;
  for (int nm=0; nm<elStates[iniN].NNormModes(); nm++) {
    single_transition.getVibrState1().addVibrQuanta(0, nm);
    single_transition.getVibrState2().addVibrQuanta(0, nm);
  }
  single_transition.getVibrState1().setElStateIndex(iniN);
  single_transition.getVibrState2().setElStateIndex(targN);

  size_t n_single_ex = node_dushinsky_rotations.find_subnode("single_excitation");
  if ( n_single_ex ) {
    if (elStates[targN].ifNMReorderedManually())
    {
      std::cout <<"WARNING! The normal modes of the target state were reordered!\n"
        <<"         New order is used for the single transitions.\n\n";
    }

    std::cout << "The following single transitions were added to the spectrum:\n" << std::flush;

    for(size_t nsex=0; nsex< n_single_ex; nsex++) {

      xml_node node_single_ex(node_dushinsky_rotations,"single_excitation",nsex);

      My_istringstream ini_str(node_single_ex.read_string_value("ini"));
      //std::cout << "Single_ex [ini]=" << ini_str.str()  << std::endl;
      fillVibrState(ini_str, single_transition.getVibrState1(), n_norm_modes);
      //std::cout << "Vibronic state 1:" << std::endl;
      //single_transition.getVibrState1().print();

      My_istringstream targ_str(node_single_ex.read_string_value("targ"));
      //std::cout << "Single_ex [targ]" << targ_str.str()  << std::endl;
      fillVibrState(targ_str, single_transition.getVibrState2(), n_norm_modes);
      //std::cout << "Vibronic state 2:" << std::endl;
      //single_transition.getVibrState2().print();

      // evaluate FCF for each transition and add to the spectrum:
      int K =single_transition.getVibrState1().getTotalQuantaCount();
      int Kp=single_transition.getVibrState2().getTotalQuantaCount();
      //std::cout << "K=" << K << " Kp=" << Kp << std::endl;

      double s_fcf=(*dushinsky_ptr).evalSingleFCF_full_space(single_transition.getVibrState1(), K, single_transition.getVibrState2(), Kp);
      (*dushinsky_ptr).addSpectralPoint(s_fcf, single_transition.getVibrState1(), single_transition.getVibrState2()); 

      std::cout << "FCF=" << std::scientific << std::setprecision(6) << s_fcf << " ";
      single_transition.getVibrState1().print();
      std::cout << "->";
      single_transition.getVibrState2().print();
      std::cout << "\n" << std::flush;
    }
  }

  std:: cout << "\nUpdating the energies and applying the Boltzmann distribution..." << std::flush;
  //--------------------------------------------------------------------------------
  //update (fill) energies for every point in the spectrum and add the Boltzmann distribution:
  int points_removed=0;
  for (int pt=0; pt<(*dushinsky_ptr).getSpectrum().getNSpectralPoints(); pt++) {
    double energy = -elStates[targN].Energy();
    double E_prime_prime = 0; // no hot bands

    // run it over the full space, if nm not in the nms_dushinsky subspace, getV_full_dim() returns zero (no excitations):
    for (int nm=0; nm <elStates[iniN].NNormModes(); nm++) {
      energy += elStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState1().getV_full_dim(nm);
      energy -= elStates[targN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState2().getV_full_dim(nm);

      E_prime_prime += 
        elStates[iniN].getNormMode(nm).getFreq() * WAVENUMBERS2EV * (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getVibrState1().getV_full_dim(nm);
    }

    // add the Boltzmann distribution to the initial state population
    double IExponent;
    if (temperature==0)
      if (E_prime_prime==0)
        IExponent=0;   // intensity unchanged 
      else
        IExponent=100; //(intensity < 10e-44 for non ground states
    else
    {
      IExponent= E_prime_prime / (temperature * KELVINS2EV );
      if (IExponent > 100) 
        IExponent=100; // keep the intensity >= 10e-44 == exp(-100)
    }
    (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getIntensity() *= exp ( -IExponent); 

    (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getEnergy()=energy;
    (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getE_prime_prime()=E_prime_prime;

    // if intensity below the fcf_threshold^2 or energy above the threshold -- do not print
    if (
        ( (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).getIntensity() < fcf_threshold*fcf_threshold ) 
        or 
        ( -(energy-E_prime_prime+elStates[targN].Energy()) > energy_threshold_target ) 
        or 
        ( E_prime_prime > energy_threshold_initial )
       )
    {
      (*dushinsky_ptr).getSpectrum().getSpectralPoint(pt).setIfPrint(false);
      points_removed++;
    }
  }
  std:: cout << "Done\n" << std::flush;

  if (max_quanta_ini!=0) {
    if (points_removed>0)                 
      std::cout << "  " << points_removed << " hot bands were removed from the spectrum\n";
    else
      std::cout << "All hot bands are above the intensity threshold\n";
    std::cout << "\n" << std::flush;
  }

  //--------------------------------------------------------------------------------
  // Print the updated spectrum:
  (*dushinsky_ptr).getSpectrum().Sort();
  std::cout << "------------------------------------------------------------------------------\n";
  std::cout << "        Stick photoelectron spectrum (with Dushinsky rotations) \n";
  std::cout << "------------------------------------------------------------------------------\n";
  if (elStates[targN].ifNMReorderedManually()) {
    std::cout <<"\nWARNING! The normal modes of the target state were reordered!\n"
      <<"         New order is used for the target state assignment.\n";
  }
  if(nms_dushinsky.size()!=n_norm_modes)
  {
    std::cout << "\nNOTE: only the following normal modes were excited (\"excite subspace\"):\n  ";
    for (int nm=0; nm<nms_dushinsky.size();nm++)
      std::cout <<nms_dushinsky[nm]<< ' ';
    std::cout << "\n";
    if (elStates[targN].ifNMReorderedManually()) {
      std::cout <<"\nWARNING! The normal modes of the target state were reordered!\n"
        <<"         New order is used for the \"excite subspace\"\n";
    }
  }
  std::cout << "\n";
  (*dushinsky_ptr).getSpectrum().PrintStickTable();
  std::cout << "------------------------------------------------------------------------------\n";

  // save the spectrum to the file
  std::stringstream spectrumFName; 
  spectrumFName << InputFileName << ".spectrum_dushinsky";
  (*dushinsky_ptr).getSpectrum().PrintStickTable(spectrumFName.str().c_str());
  std::cout << "\nStick spectrum was also saved in \"" << spectrumFName.str() << "\" file \n";
  if(nms_dushinsky.size()!=n_norm_modes)
    std::cout << " (Full list of the normal modes was used for assigning transitions)\n";
  std::cout << "\n\n";

  delete dushinsky_ptr;

}
