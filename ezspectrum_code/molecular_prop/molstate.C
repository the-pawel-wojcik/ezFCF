#include "molstate.h"

/*! \file molstate.C
  \brief Molecular state: Stores Geometry, Normal Modes & Frequencies. 
  Also reads input data from the XML file (vm'06) 
  \ingroup DATA_CLASSES
  */

//------------------------------
// align each state:
// 1. center of mass in the coordinates origin
// 2. moment of ineretia principal axes along the coordinate axes
void MolState::align()
{
  // Shift center of mass to the origin:
  shiftCoordinates( getCenterOfMass() );

  arma::Mat<double> MOI_tensor(CARTDIM, CARTDIM);
  MOI_tensor = getMomentOfInertiaTensor();

  arma::Col<double> MOI_eigenValues(CARTDIM);
  arma::Mat<double> MOI_eigenVectors(CARTDIM, CARTDIM);

  // "std" indicates use of a standard method instead of 
  // the divide-and-conquer "dc" method. The "dc" method 
  // gives slightly different results than "std", but is 
  // considerably faster for large matrices:
  // http://arma.sourceforge.net/docs.html#eig_sym
  arma::eig_sym(MOI_eigenValues, MOI_eigenVectors, MOI_tensor, "std");
  // The previous ezSpectrum convention was:
  // a) eigenvalues in descending order
  // b) eigenvectors in ROWS
  // The armadillo convention is the exact opposite.
  
  // Move eigenvectors to rows
  MOI_eigenVectors = MOI_eigenVectors.t(); 
  // Flip the order of rows
  MOI_eigenVectors = arma::flipud(MOI_eigenVectors); 
  MOI_eigenValues = arma::reverse(MOI_eigenValues);

  // compute the determinant of a MOI matrix:
  double det_tmp =  arma::det(MOI_eigenVectors);
  // if determinant is = 1 than it is a proper rotation; 
  // if it is = -1 than it is rotation+reflection (switch from left handed to the right handed coord.system); swap x&y axes:
  if (det_tmp < 0)
  {
    MOI_eigenVectors.swap_rows(0, 1);
  }

  transformCoordinates(MOI_eigenVectors); // rotate coordinates using matrix of the MOI eigen vectors 
  std::cout << "Coordinate axes were aligned with MOI principal axes; center mass was shifted to the origin.\n";
}


//------------------------------
void MolState::align(const MolState& other)
{
  //align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  // by minimizing the "sum" of distances between the same atoms: SUM[i=1..Natoms](deltaRi)
  // DeltaRi: distance between i'th atoms of the initial and target states

  double min_diff=DBL_MAX; // minimum of the "sum"
  int min_x_rot, min_y_rot, min_z_rot; // rotations, that correspond to the minimum of the "sum"

  for (int z_rot=0;  z_rot<4; z_rot++)
  {
    for (int y_rot=0;  y_rot<4; y_rot++)
    {
      for (int x_rot=0;  x_rot<4; x_rot++)
      {
        double diff = getGeomDifference(other);
        if (diff <= min_diff) {
          min_diff=diff;
          min_x_rot=x_rot;
          min_y_rot=y_rot;
          min_z_rot=z_rot;
        }
        rotateX_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
      }
      rotateY_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry and starts over.
    }
    rotateZ_90deg(); // after 4 rotations by pi/2, it gets the "0" rotation geometry.
  }// at the end the molecule rotated by 4*pi/2 around x and y, i.e. it is non-rotated

  // rotate to the angles with the min norm of the deltaRi differences (min_x_rot, min_y_rot, min_z_rot):
  rotate(min_x_rot*PI/2.0, min_y_rot*PI/2.0, min_z_rot*PI/2.0);

  std::cout << "Also rotated by " << min_x_rot << "/2*pi, " << min_y_rot
            << "/2*pi, and " << min_z_rot << "/2*pi CCW around x, y, and z.\n";
  std::cout << "The norm of the geometry difference from the initial state is "
            << sqrt(min_diff) << "\n";
}

//------------------------------
bool MolState::ifSimilar(MolState& other)
{
  bool return_bool=true;
  std::string error_str="";
  if ( NAtoms()!=other.NAtoms() )
  {
    std::cout << "\nError: different number of atoms in the initial and the target states\n\n";
    return_bool=false;	
  }
  else if ( NNormModes()!=other.NNormModes() )
  {
    std::cout << "\nError: different number of normal modes in the initial and the target state\n\n";
    return_bool=false;	
  }

  return return_bool;
}


//------------------------------
arma::Col<double>& MolState::getCenterOfMass()
{
  centerOfMass = arma::Col<double> (CARTDIM, arma::fill::zeros);
  double totalMass = 0.0;
  for (int i=0; i<atoms.size(); i++)
  {
    for (int axis=0; axis<3; axis++)
      centerOfMass(axis) += atoms[i].getCoordMass(axis);
    totalMass+=atoms[i].getMass();
  }
  centerOfMass /= totalMass;
  return centerOfMass;
}

//------------------------------
arma::Mat<double> &MolState::getMomentOfInertiaTensor() {
  momentOfInertiaTensor =
      arma::Mat<double>(CARTDIM, CARTDIM, arma::fill::zeros);

  // I_ij = SUM_k[m_k * (r_k^2*delta_ij-r_ki*r_kj)]
  for (int i = 0; i < CARTDIM; i++) {
    for (int j = 0; j < CARTDIM; j++) {
      for (int atom = 0; atom < NAtoms(); atom++) {
        double atom_mass = atoms[atom].getMass();
        momentOfInertiaTensor(i, j) -=
            atoms[atom].getCoord(i) * atoms[atom].getCoord(j) * atom_mass;
        // add R^2, i.e. I_zz = m * (R^2 - r_z^2) = m * (r_x^2 + r_y^2)
        if (i == j)
          momentOfInertiaTensor(i, j) +=
              atoms[atom].getR() * atoms[atom].getR() * atom_mass;
      }
    }
  }

  // clean zeroes
  momentOfInertiaTensor.clean(MOI_THRESHOLD);

  return momentOfInertiaTensor;
}

//------------------------------
void MolState::shiftCoordinates(arma::Col<double>& vector)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].shiftCoordinates(vector);
  // there is no need to shift normal coordinates!!
}

//------------------------------
void MolState::transformCoordinates(const arma::Mat<double>& matrix_3x3)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].transformCoordinates(matrix_3x3);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].transformCoordinates(matrix_3x3);
}

//------------------------------
void MolState::rotateX_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateX_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateX_90deg();
}

//------------------------------
void MolState::rotateY_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateY_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateY_90deg();
}

//------------------------------
void MolState::rotateZ_90deg()
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].rotateZ_90deg();
  for (int i=0; i<normModes.size(); i++)
    normModes[i].rotateZ_90deg();
}

//------------------------------
void MolState::rotate(const double alpha_x, const double alpha_y, const double alpha_z)
{
  // three rotation matrices (instead of making one matrix) around x, y, z axes
  arma::Mat<double> Rx(CARTDIM, CARTDIM, arma::fill::zeros);
  Rx(0, 0) = 1;
  Rx(1, 1) = cos(alpha_x);
  Rx(2, 2) = cos(alpha_x);
  Rx(2, 1) = -sin(alpha_x);
  Rx(1, 2) = sin(alpha_x);

  arma::Mat<double> Ry(CARTDIM, CARTDIM, arma::fill::zeros); 
  Ry(1, 1) = 1;
  Ry(0, 0) = cos(alpha_y);
  Ry(2, 2) = cos(alpha_y);
  Ry(2, 0) = sin(alpha_y);
  Ry(0, 2) = -sin(alpha_y);

  arma::Mat<double> Rz(CARTDIM, CARTDIM, arma::fill::zeros);
  Rz(2, 2) = 1;
  Rz(0, 0) = cos(alpha_z);
  Rz(1, 1) = cos(alpha_z);
  Rz(1, 0) = -sin(alpha_z);
  Rz(0, 1) = sin(alpha_z);

  // Rotations are applied in the order: x, y, and z
  // which gives the overall rotation matrix
  // R = Rx Ry Rz
  arma::Mat<double> R = arma::eye(CARTDIM, CARTDIM);
  R *= Rx;
  R *= Ry;
  R *= Rz;

  // and now rotates using matix R:
  transformCoordinates(R);
}
//------------------------------
bool MolState::ifAlignedManually()
{
  return if_aligned_manually;
}
//------------------------------
bool MolState::ifNMReorderedManually() const
{
  return if_nm_reordered_manually;
}

//------------------------------
double MolState::getGeomDifference(const MolState &other) const {
  double diff = 0;
  for (int i = 0; i < atoms.size(); i++)
    for (int j = 0; j < CARTDIM; j++) {
      // diff += DeltaR ^ 2 [=DetaX^2+DeltaY^2+DeltaZ^2]
      double coord_diff = getAtom(i).getCoord(j) - other.getAtom(i).getCoord(j);
      diff += coord_diff * coord_diff;
    }
  return diff;
}

//------------------------------
void MolState::applyCoordinateThreshold(const double threshold)
{
  for (int i=0; i<atoms.size(); i++)
    atoms[i].applyCoordinateThreshold(threshold);
  for (int i=0; i<normModes.size(); i++)
    normModes[i].applyCoordinateThreshold(threshold);
}

bool MolState::getNormalModeOverlapWithOtherState(MolState& other, arma::Mat<double>& overlap, std::vector<int>& normal_modes_list)
{
  overlap = arma::Mat<double> (NNormModes(), NNormModes(), arma::fill::zeros);

  //normalization of initial and target normal modes (stored as un-mass-weighted)
  double norm_ini, norm_targ;
  // norm of the L
  for (int nm1=0; nm1<NNormModes(); nm1++)
    for (int nm2=0; nm2<NNormModes(); nm2++)
    {
      norm_ini = 0;
      norm_targ = 0;
      for (int a=0; a<NAtoms(); a++)
        for (int i=0; i<CARTDIM; i++ )
        {
          double dispnm1 = getNormMode(nm1).getDisplacement()(a*CARTDIM+i);
          double dispnm2 = other.getNormMode(nm2).getDisplacement()(a*CARTDIM+i);
          //x1*x1+y1*y1+..
          norm_ini+= dispnm1 * dispnm1;
          //x2*x2+y2*y2+..
          norm_targ+= dispnm2 * dispnm2;
          //x1*x2+y1*y2+...
          overlap(nm1, nm2)+= dispnm1 * dispnm2;
        }
      overlap(nm1, nm2)/= sqrt(norm_ini) * sqrt(norm_targ);
    }
  bool return_bool=true;


  // scan the diagonal; if the diagonal element is NOT the maximum element in the raw and the colum -- than ADD this nm to the list;
  double max;
  for (int nm=0; nm<NNormModes(); nm++)
  {
    max=fabs(overlap(nm, nm)); // maximum should be at the diagonal
    for (int nm_scan=0; nm_scan<NNormModes(); nm_scan++) // scan row&column
      if ( (fabs(overlap(nm, nm_scan))>max) or (fabs(overlap(nm_scan, nm))>max) ) // check nm-th row & column
      {
        return_bool=false;
        normal_modes_list.push_back(nm);
        normal_modes_list.push_back(nm_scan); // the column/row with larger element should be also included
      }
  }

  // Now Sort and Remove duplicates from the normal_modes_list:
  std::sort( normal_modes_list.begin(), normal_modes_list.end() );
  std::vector<int>::iterator new_end_pos;
  new_end_pos = std::unique( normal_modes_list.begin(), normal_modes_list.end() );
  normal_modes_list.erase( new_end_pos, normal_modes_list.end() );

  return return_bool;
}



void MolState::Print()
{
  std::cout << "Geometry=\n";
  printGeometry();    
  std::cout << "Normal modes=\n";
  for (int k=0; k<NNormModes(); k++)
  {
    std::cout << "   Frequency="; 
    std::cout << getNormMode(k).getFreq() << '\n';
    std::cout << "   Displacement=\n";
    for (int i=0; i<NAtoms(); i++)
    {
      for (int l=0; l<CARTDIM; l++)
        std::cout << getNormMode(k).getDisplacement()(i*CARTDIM+l) << ' ';
      std::cout << '\n';
    }
  }
  std::cout <<" end of the electronic state \n";   
}

void MolState::printGeometry()
{
  for (int i = 0; i < NAtoms(); i++) {
    std::cout << std::setw(4) << std::right << getAtom(i).Name();
    for (int k = 0; k < CARTDIM; k++)
      std::cout << std::setw(12) << std::right << std::fixed
                << std::setprecision(4) << std::showpoint << getAtom(i).Coord(k)
                << ' ';
    std::cout << '\n';
  }
}

arma::Col<double> MolState::get_geometry_as_col()
{
  arma::Col<double> geo_vec(this->NAtoms() * CARTDIM, arma::fill::zeros);

  for (int atm_idx = 0; atm_idx < NAtoms(); ++atm_idx) {
    for (int xyz_idx = 0; xyz_idx < CARTDIM; ++xyz_idx)
      geo_vec(atm_idx * 3 + xyz_idx) = getAtom(atm_idx).Coord(xyz_idx);
  }

  return geo_vec;
}

void MolState::printNormalModes()
{
  // print in qchem format (3 per line)
  int nModesPerLine=3;  // three number of vib. modes per Line
  int nLines = NNormModes() / nModesPerLine;
  if ( NNormModes() % nModesPerLine != 0)
    nLines++;

  for (int n = 0; n < nLines; n++)  // number of blocks with 3 norm.modes. ("lines")
  {
    int current_nModesPerLine = nModesPerLine;
    // for the last entree, nModesPerString may differ from 3
    if (nLines - 1 == n)
      if ( NNormModes() % nModesPerLine != 0 )
        current_nModesPerLine = NNormModes() % nModesPerLine;

    for (int a=0; a < NAtoms(); a++)   
    {
      for (int j=0; j < current_nModesPerLine; j++) 
      {
        for (int k=0; k < CARTDIM; k++)

          std::cout 
            << std::setw(7) << std::right << std::fixed << std::setprecision(3)
            << 
            getNormMode(n*nModesPerLine+j).getDisplacement()(a*CARTDIM+k) 
            * 
            sqrt( reduced_masses(n * nModesPerLine + j) / getAtom(a).Mass() )
            <<' ';
        std::cout <<  "  ";
      }
      std::cout << "\n";
    }
    std::cout << "\n";
  }
}

void MolState::printGradient()
{
  for (int i=0; i<NAtoms(); i++)
  {
    std::cout << std::setw(4) << std::right  << getAtom(i).Name();
    for (int k=0; k<CARTDIM; k++)
      std::cout << std::setw(12) << std::right << std::fixed << std::setprecision(4) 
        << std::showpoint << gradient(i * CARTDIM + k) << ' '; 
    std::cout << '\n';
  }
}

// Helper of 'Read_excitation_energy' and 'Read_vertical_energy' functions
void check_and_convert_energy(double &energy, const std::string &energy_text,
                              const std::string &units) {
  // print energy with its unit before conversion to eV
  if (units != "eV") {
    std::cout << energy_text << std::fixed << std::setprecision(6) << energy
              << " " << units << std::endl;
  }

  if (!covert_energy_to_eV(energy, units)) {
    std::cout << "\nError! Unknown units of the excitation energy: \"" << units
              << "\"\n  (should be equal to \"eV\", \"K\", or \"cm-1\")\n\n";
    exit(1);
  }

  // print energy in eV
  std::cout << energy_text << std::fixed << std::setprecision(6) << energy
            << " eV" << std::endl;
}

/* Parser of the "excitation_energy" subnode of the "initial_state" or
 * "target_state" node. Assigns value to the `energy` variable. Helper of the
 * `Read` function. If the node is not presents (e.g. `initial_state`), the
 * `energy` variable is set to zero. */
void MolState::Read_excitation_energy(xml_node &node_state) {

  // adiabatic ionization/excitation energy (bottom to bottom), previously
  // called IP
  energy = 0.0; // eV

  // If the node is not present there is nothing to read
  if (! node_state.find_subnode("excitation_energy"))
    return;

  // If the node is used together with the 'gradient' node -> it's an error
  if (IfGradientAvailable) {
    error("Vertical gradient method uses the 'vertical_excitation_energy' "
          "input node instead of the 'excitation_energy' node which is "
          "reserved for other types of calculations.");
  }

  xml_node node_eenergy(node_state, "excitation_energy", 0);
  std::string units = node_eenergy.read_string_value("units");
  energy = node_eenergy.read_node_double_value();
  std::string energy_text = "Excitation energy = ";

  check_and_convert_energy(energy, energy_text, units);
}

/* Parser of the "vertical_energy" subnode of the "target_state" node. Helper of
 * the `Read` function. Assigns values to the `vertical_energy` variable. If the
 * node is not presents the variable is set to zero. If the
 * "vertical_excitation_energy" node is present, this function returns error
 * unless the "gradient" node is also present in the state node. */
void MolState::Read_vertical_energy(xml_node &node_state) {

  // vertical excitation energy
  vertical_energy = 0.0; // eV
                         //
  // If the node is not present there is nothing to read
  if (! node_state.find_subnode("vertical_excitation_energy"))
    return;

  if (! IfGradientAvailable) {
    error("Use `vertical_excitation_energy` node only in calculations using "
          "the vertical gradient approximation.");
  }
  xml_node node_venergy(node_state, "vertical_excitation_energy", 0);
  std::string units = node_venergy.read_string_value("units");
  vertical_energy = node_venergy.read_node_double_value();
  std::string energy_text = "Vertical excitation energy = ";

  check_and_convert_energy(vertical_energy, energy_text, units);
}

/* Parser of the "geometry" subnode of the "initial_state" or "target_state"
 * node. Helper of the `Read` function. Assigns value to the `n_molecular_nm`,
 * and `atoms` variables. The `atoms` variable stores three cartesian components
 * of each atom (in Angstroms) and the name of the atom, the name is important
 * as it is used in conjunction with 'atomicMasses.xml' to tell the atom's mass.
 */
void MolState::Read_molecular_geometry(xml_node &node_state) {

  // If the node is not present there is nothing to read
  if (! node_state.find_subnode("geometry"))
    return;

  xml_node node_geometry(node_state, "geometry", 0);

  int nAtoms;
  nAtoms = node_geometry.read_int_value("number_of_atoms");

  n_molecular_nm = 3 * nAtoms - 6;
  bool ifLinear = node_geometry.read_bool_value("linear");
  if (ifLinear) {
    n_molecular_nm = 3 * nAtoms - 5;
  }

  My_istringstream geom_iStr(node_geometry.read_string_value("text"));
  atoms.clear();

  std::string units = node_geometry.read_string_value("units");
  double coeff = (units == "au") ? AU2ANGSTROM : 1.0;

  Atom tmp_atom;
  for (int i = 0; i < nAtoms; i++) {
    // Get atomic name:
    std::string tmp_atomName;
    geom_iStr.getNextWord(tmp_atomName);
    tmp_atom.Name() = tmp_atomName;

    // get coordinates & convert to Angstroms
    tmp_atom.Coord(0) = geom_iStr.getNextDouble() * coeff;
    tmp_atom.Coord(1) = geom_iStr.getNextDouble() * coeff;
    tmp_atom.Coord(2) = geom_iStr.getNextDouble() * coeff;

    if (geom_iStr.fail()) {
      std::cout << "\nMolState::Read_molecular_geometry(): Error. Incorrect "
                   "format of geometry:\n" +
                       geom_iStr.str() + "\n";
      exit(1);
    }
    atoms.push_back(tmp_atom);
  }
}

// Helper of `Read_normal_modes`
void MolState::Read_abinitio_atoms_masses(std::string &atoms_text) {
  My_istringstream atoms_iStr(atoms_text);

  ab_intio_atoms_masses.clear();
  Atom tmp_atom;
  for (int i = 0; i < NAtoms(); i++) {
    std::string tmp_atomName;
    // Get atomic name:
    atoms_iStr.getNextWord(tmp_atomName);
    tmp_atom.Name() = tmp_atomName;
    ab_intio_atoms_masses.push_back(tmp_atom);
  }
}

/* Parser of the "normal_modes" subnode of the "initial_state" or "target_state"
 * node. Normal modes are stored in the mass-un-weiged cartesian coordinates.
 * At this point the nms are read from the input xml file where sometimes they
 * are still in the mass-weighted form, for these cases the function
 * `un_mass_weight_nm` is later triggered by the variable
 * `ifInputNMmassweighted`.
 *
 * Helper of the `Read` function.
 * Assigns values to the `normModes`, `ifInputNMmassweighted`, and
 * `ab_intio_atoms_masses` variables.*/
void MolState::Read_normal_modes(xml_node &node_state) {

  // If the node is not present there is nothing to read
  if (! node_state.find_subnode("normal_modes"))
    return;

  xml_node node_nmodes(node_state, "normal_modes", 0);

  ifInputNMmassweighted = node_nmodes.read_bool_value("if_mass_weighted");

  std::string nm_atoms_list = node_nmodes.read_string_value("atoms");
  Read_abinitio_atoms_masses(nm_atoms_list);

  // the key node atribute
  My_istringstream nmodes_iStr(node_nmodes.read_string_value("text"));

  // number of normal modes per Line in the xml input
  int nModesPerLine = 3;
  int nLines;
  nLines = n_molecular_nm / nModesPerLine;
  if (n_molecular_nm % nModesPerLine != 0)
    nLines++;

  // initilize the container to empty values
  NormalMode nMode(NAtoms(), 0);
  normModes = std::vector<NormalMode> (n_molecular_nm, nMode);
  // number of blocks with 3 norm.modes. ("lines")
  for (int k = 0; k < nLines; k++) {
    int current_nModesPerLine = nModesPerLine;

    // in the last entry, nModesPerString may be less than 3
    if (nLines - 1 == k)
      if (n_molecular_nm % nModesPerLine != 0)
        current_nModesPerLine = n_molecular_nm % nModesPerLine;

    for (int i = 0; i < NAtoms(); i++)
      for (int j = 0; j < current_nModesPerLine; j++)
        for (int l = 0; l < CARTDIM; l++) {
          normModes[k * nModesPerLine + j].getDisplacement()(i * CARTDIM + l) =
              nmodes_iStr.getNextDouble();
          if (nmodes_iStr.fail()) {
            std::cout << "\nMolState::Read_normal_modes(): Error. Wrong format "
                         "of normal modes:\n" +
                             nmodes_iStr.str() + "\n";
            exit(1);
          }
        }
  }
}

/* Parser of the "frequencies" subnode of the "initial_state" or "target_state"
 * node. Helper of the `Read` function. Assigns frequency values to the
 * `NormalModes` in the `normModes` container. */
void MolState::Read_frequencies(xml_node &node_state) {

  // If the node is not present there is nothing to read
  if (! node_state.find_subnode("frequencies"))
    return;

  xml_node node_freq(node_state, "frequencies", 0);
  My_istringstream freq_iStr(node_freq.read_string_value("text"));

  for (int i = 0; i < n_molecular_nm; i++) {
    getNormMode(i).getFreq() = freq_iStr.getNextDouble();
    if (freq_iStr.fail()) {
      std::cout << "\nError: format error at frequency #" << i
                << " or less than " << n_molecular_nm << " frequencies found."
                << "\n";
      exit(1);
    }
    if (getNormMode(i).getFreq() <= 0) {
      std::cout << "\nError. The frequency [" << getNormMode(i).getFreq()
                << "] is negative\n";
      exit(1);
    }
  }
}

/* Parser of the "gradient" subnode of the "initial_state" or "target_state"
 * node. Helper of the `Read` function. Assigns values to the
 * 'IfGradientAvailable' and `gradient` (if available) variables. If not
 * available the `gradient` variable is left unassigned. */
void MolState::Read_vertical_gradient(xml_node &node_state) {

  IfGradientAvailable = bool(node_state.find_subnode("gradient"));

  // if the gradient node is not present in the input, this is the end
  if (! IfGradientAvailable) 
    return;

  xml_node node_gradient(node_state, "gradient", 0);

  // Read out the gradient in atomic units (Q-Chem default)
  std::string units = node_gradient.read_string_value("units");
  trim(units);
  if (units != "a.u.") {
    std::cout << "\nError! Gradient reported in unsupported units: \"" << units
              << "\"\n  (use \"a.u.\" or contact us to add a support of " 
	      << units << ")\n\n";
    exit(1);
  }

  std::string grad_text = node_gradient.read_string_value("text");
  std::istringstream grad_istr(trim(grad_text));

  // Read gradient withouth assuming number of it's components. This approach
  // allows to detect incorrect input no only if too few coordinates are given
  // but also if there were too many. This is checked in the
  // `test_vertical_gradient_dimension` function.
  std::queue<double> grad_queue;
  double cartesian_component;
  while (grad_istr.good()) {
    // Read the x, y, and z components of the gradient for every atom
    grad_istr >> cartesian_component;
    if (grad_istr.fail()) {
      std::cout << "\nError! Wrong gradient format:\n"
                << grad_istr.str() << std::endl;
      exit(1);
    }
    grad_queue.push(cartesian_component);
  }
  // copy the gradient into MolState storage
  size_t gradient_length = grad_queue.size();
  // gradient is stored as a column vector
  gradient = arma::Col<double>(gradient_length);
  for (int i = 0; i < gradient_length; ++i) {
    gradient(i) = grad_queue.front();
    grad_queue.pop();
  }
}

/* Parser of the "manual_coord_transform" subnode of the "initial_state" or
 * "target_state" node. Helper of the `Read` function. Assigns values to the
 * `manual_transofrmations` variable and sets 'if_aligned_manually' to false. */
void MolState::Read_manual_coord_transformations(xml_node &node_state) {
  if_aligned_manually = false;

  size_t manual_coord_transform =
      node_state.find_subnode("manual_coordinates_transformation");

  typedef arma::Col<double> vec;
  for (int i = 0; i < manual_coord_transform; ++i) {
    xml_node mct_node(node_state,"manual_coordinates_transformation", i);

    vec shift(CARTDIM, arma::fill::zeros);
    shift(0) -= mct_node.read_double_value("shift_along_x");
    shift(1) -= mct_node.read_double_value("shift_along_y");
    shift(2) -= mct_node.read_double_value("shift_along_z");

    vec rotation(CARTDIM, arma::fill::zeros);
    rotation(0) = mct_node.read_double_value("rotate_around_x");
    rotation(1) = mct_node.read_double_value("rotate_around_y");
    rotation(2) = mct_node.read_double_value("rotate_around_z");

    manual_transformations.push(std::pair<vec, vec>(shift, rotation));
  }
}

/* Parser of the "manual_normal_modes_reordering" subnode of the "initial_state"
 * or "target_state" node. Helper of the `Read` function. Assigns values to the 
 * `normModesOrder` and 'if_nm_reordered_manually` variables. */
void MolState::Read_normal_modes_reorder(xml_node &node_state) {

  if_nm_reordered_manually = false;
  size_t reorder_nmodes =
      node_state.find_subnode("manual_normal_modes_reordering");

  if (reorder_nmodes) {

    xml_node node_nm_order(node_state, "manual_normal_modes_reordering", 0);
    std::istringstream order_istr(node_nm_order.read_string_value("new_order"));

    normModesOrder.clear();

    std::cout << "New normal modes order was requested:\n"
              << order_istr.str() << "\n";
    int mode;
    for (int nm = 0; nm < NNormModes(); nm++) {
      order_istr >> mode;

      if (order_istr.fail()) {
        std::cout << "\nFormat error: non-numeric symbol or less entries then "
                     "the number of normal modes\n\n";
        exit(1);
      }
      if ((mode < 0) or (mode >= NNormModes())) {
        std::cout << "\nError: normal mode number [" << mode
                  << "] is out of range [0.." << NNormModes() - 1 << "].\n\n";
        exit(1);
      }
      normModesOrder.push_back(mode);
    }

    // check if there are duplicates in the list:
    std::set<int> s(normModesOrder.begin(), normModesOrder.end());
    if (s.size() != normModesOrder.size()) {
      std::cout << "\nFormat error: manual_normal_modes_reordering node "
                   "contains duplicates.\n";
      exit(0);
    }

    if_nm_reordered_manually = true;

  } else {
    for (int nm = 0; nm < NNormModes(); nm++)
      normModesOrder.push_back(nm);
  }
}

/* Parser of the "manual_atoms_reordering" subnode of the "initial_state" or
 * "target_state" node. Helper of the `Read` function. Assigns values to the
 * `atomOrder` and 'if_atoms_reordered_manually` variables. */
void MolState::Read_atoms_reorder(xml_node &node_state) {
  if_atoms_reordered_manually =
      bool(node_state.find_subnode("manual_atoms_reordering"));

  if (if_atoms_reordered_manually == false) {
    for (int i = 0; i < NAtoms(); ++i ) {
      atomsOrder.push_back(i);
    }
    // all job done leave the function
    return; 
  }

  xml_node reorder_atoms(node_state, "manual_atoms_reordering", 0);

  std::istringstream order_istr(reorder_atoms.read_string_value("new_order"));
  std::cout << "New atoms order:\n" << order_istr.str() << "\n";

  int atom;
  for (int nm = 0; nm < NAtoms(); nm++) {
    order_istr >> atom;
    if (order_istr.fail()) {
      std::cout << "\n"
                << "Input error in manual_atoms_reordering. Cannot process:\n  "
                << order_istr.str() << "\n";
      exit(1);
    }
    if ((atom < 0) or (atom >= NAtoms())) {
      std::cout << "\n"
                << "Input error in manual_atoms_reordering.\n"
                << "  atom number [" << atom << "] is out of range [0.."
                << NAtoms() - 1 << "].\n";
      exit(1);
    }
    atomsOrder.push_back(atom);
  }

  // check if there are duplicates in the list:
  std::set<int> s(atomsOrder.begin(), atomsOrder.end());
  if (s.size() != atomsOrder.size()) {
    std::cout << "\n"
              << "Input error in manual_atoms_reordering.\n"
              << "  Repeated entries in manual_atoms_reordering:\n  ";
    for (int atm : atomsOrder)
      std::cout << ' ' << atm << std::endl;
    exit(1);
  }
}

/* The function that read all keywords passed to the "initial_state" or
   "target_state" nodes. Arguments:
   -  a node pointing out to initial_state or target_state node
   -  a file where all masses are tabulated
   Collected input is saved into private variables of `MolState` processed
   as need in the `ApplyKeyWords` function.
   */
void MolState::Read(xml_node &node_state) {
  // Presence of the "gradient" node triggers the vertical gradient calculations
  Read_vertical_gradient(node_state);
  Read_excitation_energy(node_state);
  Read_vertical_energy(node_state);
  Read_molecular_geometry(node_state);
  Read_normal_modes(node_state);
  Read_frequencies(node_state);
  Read_manual_coord_transformations(node_state);
  Read_normal_modes_reorder(node_state);
  Read_atoms_reorder(node_state);
}

/* MolState object stores two vectors of atoms: one from the "geometry" node and
 * one from the "normal_modes" node. After sucessful reading of the xml input,
 * each vector stores the atom names. This function uses the atomic_masses.xml
 * file to add the information about the atom mass to the Atom containers. */
void MolState::convert_atomic_names_to_masses(xml_node &node_amu_table) {
  for (int i = 0; i < NAtoms(); i++)
    atoms[i].Mass() =
        node_amu_table.read_node_double_value(atoms[i].Name().c_str());

  for (int i = 0; i < NAtoms(); i++)
    ab_intio_atoms_masses[i].Mass() =
        node_amu_table.read_node_double_value(ab_intio_atoms_masses[i].Name().c_str());
}

/* Some ab-initio programs report mass-weighed normal modes while other report
 * the regular cartesian coordinates. MolState stores the normal modes in the
 * regual (mass-un-weighed) form. This function transforms the mass-weighed
 * normal modes into the mass-un-weighed nodes used in the reminder of the
 * program.
 *
 * Helper of the `ApplyKeyWords` function. Assigns value to the `reduced_masses`
 * variable. */
void MolState::un_mass_weight_nm() {
  // Mass-un-weight normal modes, using the ab_intio_atoms_masses
  for (int nm = 0; nm < NNormModes(); nm++)
    for (int a = 0; a < NAtoms(); a++) {
      double sqrt_mass = sqrt(ab_intio_atoms_masses[a].Mass());
      for (int i = 0; i < CARTDIM; i++)
        getNormMode(nm).getDisplacement()(a * CARTDIM + i) *= sqrt_mass;
    }

  // Normalize each normal mode (divide by sqrt(norm) which is the same as
  // division by sqrt(reduced mass))

  // Prepare reduced masses:
  reduced_masses = arma::Col<double>(NNormModes(), arma::fill::zeros);
  for (int nm = 0; nm < NNormModes(); nm++) {
    for (int a = 0; a < NAtoms(); a++)
      for (int i = 0; i < CARTDIM; i++) {
        double nm_displacement = getNormMode(nm).getDisplacement()(a * CARTDIM + i);
        reduced_masses(nm) += nm_displacement * nm_displacement;
      }
  }

  // Normalize:
  for (int nm = 0; nm < NNormModes(); nm++) {
    double sqrt_rm = sqrt(reduced_masses(nm));
    for (int a = 0; a < NAtoms(); a++)
      for (int i = 0; i < CARTDIM; i++)
        getNormMode(nm).getDisplacement()(a * CARTDIM + i) /= sqrt_rm;
  }
}

/* Helper of the `ApplyKeyWords` function.
 * Assigns variables:
 * 1) mass_matrix
 * 2) omega_matrix
 * 3) d_matrix
 *
 * `mass_matrix` is a diagonal matrix of atomic masses of dimension 3*(# atoms).
 * Masses are stored in AU (not amu!). `mass_matrix` uses the order of atoms
 * from `atoms` which in turn are read from the "geometry" node of an input xml
 * file. These masses can be defined by user to account for an isotope effect.
 *
 * The diagonal matrix `omega_matrix` stores harmonic frequencies also in AU
 * (not cm-1), it has the dimension matching the number of normal modes,
 * `n_molecular_nm`.
 *
 * `d_matrix` is a rectangular matrix that stores the normal modes as its
 * columns which implies the matrix dimensions: 3*(#atoms) x n_molecular_nm.
 * Normal modes matrix `d_matrix` is dimensionless. */
void MolState::create_matrices() {
  int threeN = CARTDIM * NAtoms();
  mass_matrix = arma::Mat<double>(threeN, threeN, arma::fill::zeros);

  for (int i = 0; i < NAtoms(); i++) {
    double mass = atoms[i].Mass() * AMU_2_ELECTRONMASS; // mass in AU
    for (int j = 0; j < CARTDIM; j++) {
      int index = CARTDIM * i + j;
      mass_matrix(index, index) = mass;
    }
  }

  // A diagonal matrix that stores the harmonic frequencies
  omega_matrix =
      arma::Mat<double>(n_molecular_nm, n_molecular_nm, arma::fill::zeros);
  for (int j = 0; j < n_molecular_nm; j++) {
    double freq = normModes[j].getFreq();
    freq *= WAVENUMBERS2EV * EV2HARTREE; // frequencies in AU
    omega_matrix(j, j) = freq;
  }

  // Matrix that stores normal modes as columns
  d_matrix = arma::Mat<double>(threeN, n_molecular_nm, arma::fill::zeros);
  for (int j = 0; j < n_molecular_nm; j++) {
    for (int i = 0; i < CARTDIM * NAtoms(); i++) {
      d_matrix(i, j) = normModes[j].getDisplacement()(i); // D is dimensionless
    }
  }
}

/* The `Read_vertical_gradient` function reads the gradient node input without
 * restricting the number of elements that are parsed. If it contains too few or
 * too many coordinates, this function captures it.*/
void MolState::test_vertical_gradient_dimension() {
  size_t gradient_dim = gradient.n_rows;
  if (gradient_dim != CARTDIM * NAtoms()) {
    std::cout << " Error. Gradient dimension (" << gradient_dim
              << ") doesn't match 3*(# atoms) = " << CARTDIM * NAtoms() << "\n";
    gradient.print("Gradient");
    exit(1);
  }
}

/* The `Read_vertical_gradient` function reads the gradient node input without
 * restricing the number of elements that are parsed. If it contains too few or
 * too many coordinates, this function captures it.*/
void MolState::test_for_small_frequencies_in_VGA() {

  auto min_nmod = std::min_element(normModes.begin(), normModes.end());
  double min_frequency = min_nmod->getFreq();

  if (min_frequency <= SMALL_FREQ_WARNING) {
    std::cout << " Warning!\n"
                 "State contains small frequencies (<"
              << std::fixed << std::setw(4) << std::setprecision(0)
              << SMALL_FREQ_WARNING
              << " wavenumbers)."
                 "\nVGA can overestimate displacement along the "
                 "corresponding modes.\n"
                 "Consider a separate treatment of these nodes.\n"
              << std::endl;
  }
}

/* Copies the bare minimum of data from the initial state for the VG
 * calculations. */
void MolState::copy_data_from_the_initial_state(const MolState &is){
  std::cout << "\n Using vertical gradient approximation (VGA).\n"
               "Properties copied from the initial state:\n"
               " - molecular geometry\n"
               " - normal modes\n"
               " - frequencies.\n\n";
  
  atoms = is.atoms;
  normModes = is.normModes;
  n_molecular_nm = is.n_molecular_nm;
  ab_intio_atoms_masses = is.ab_intio_atoms_masses;
  normModesOrder = is.normModesOrder;
  reduced_masses = is.reduced_masses;
  mass_matrix = is.mass_matrix;
  omega_matrix = is.omega_matrix;
  d_matrix = is.d_matrix;
}

/* Helper of the `vertical_gradient_method` function. Calculates the VG
 * geometry. 
 * `mmmh` is matrix of atomic masses to power minus half
 * `d` is a matrix with normal modes as its columns
 * `delta` is vector with displacement along normal modes
 * See comments in the `vertical_gradient_method` for more detail. */
void MolState::vg_calc_geom(const arma::Mat<double> &mmmh,
                            const arma::Mat<double> &d,
                            const arma::vec &delta) {
  // R _e ^{(2)} = R _e ^{(1)} - M ^{-1/2} D \Delta
  arma::Col<double> geometry_shift;
  geometry_shift = mmmh * d * delta;
  geometry_shift *= AU2ANGSTROM;

  for (int i = 0; i < NAtoms(); i++) {
    for (int j = 0; j < CARTDIM; j++) {
      atoms[i].Coord(j) -= geometry_shift(CARTDIM * i + j);
    }
  }
  std::cout << "Target-state VGA geometry:" << std::endl;
  printGeometry();
}

/* Helper of the `vertical_gradient_method` function. Calculates the correction
 * to the adiabatic excitation energy coming from the VG method.
 * `omm2` is the diagonal matrix of harmonic frequencies to power minus 2.
 * `delta` is a vector of geometry displacement in normal modes (see comments to
 * the `vertical_gradient_method`.) */
void MolState::vg_calc_energy(const arma::vec &delta,
                              const arma::Mat<double> &omm2) {
  double vg_e_shift = 0.0;
  for (int i = 0; i < NNormModes(); i++) {
    vg_e_shift -= 0.5 * delta(i) * delta(i) / omm2(i, i) / EV2HARTREE;
  }
  energy = vertical_energy + vg_e_shift;

  // HINT: Vertical gradient works within parallel approximation wo/ frequency
  // shifts. The adiabatic excitation energy is equal to the E^a _{00}, i.e.,
  // ZPE of the initial to the ZPE of the target states. This is the energy used
  // in the guts of the program where the intentities are calculated.
  std::cout << "\nAdiabatic excitation energy (within VGA) = " << energy << " eV "
            << std::endl;
}

/* Helper of the `ApplyKeyWords` function. Main part of the Vertical Gradient
 * implementation. Calculates the target state geometry and excitation energy
 * which can be forwarded to the parallel mode approximation code. */
void MolState::vertical_gradient_method() {
  /* The equilibrium geometry of the second (target) state) $R _e ^(2)$ is
   * derived from a vector $\Delta$ (defined below)
   * $$ R _e ^(2) = R _e ^{(1)} - M ^{-1/2} D \Delta $$
   * or in terms of a gradient
   * $$
   * R _e ^(2)
   * =
   * R _e ^{(1)}
   * -
   * M ^{-1/2} D \Omega ^{-2} D ^{-1} M ^{-1/2} g _{(2)} ^{(x)}.
   * $$
   *
   * The vector $\Delta$ is defined as
   * $$ \Delta = \Omega ^{-2} D ^{-1} M ^{-1/2} g _{(2)} ^{(x)} $$
   * it is a vector that allows for a transformation between the target
   * state normal coordinates $q ^{(2)}$ and the initial state normal
   * coordinates $q ^{(1)}$:
   * $$ q ^{(2)} = q ^{(1)} + \Delta $$
   * (in parallel approximation without frequency shifts there is no
   * Duschinsky matrix in the last equation).
   *
   * The two sets of normal modes are parallel, therfore $\Delta$ is also the
   * shift vector between two equilibrium geometries expressed in the normal
   * mode coordinates of the (1) state.
   *
   * $\Omega$ is a diagonal matrix (3N - 5/6 x 3N - 5/6) of harmonic
   * frequencies. $D$ stores normal modes as its columns. $D$ is a rectangular
   * matrix (3N x 3N - 5/6) that turns the the mass-weighted Hessian matrix
   * into its diagonal form:
   * $$ H = D \Omega ^2 D ^T $$
   * One dimension of $D$ is shorter because the normal coordinates
   * coresponding to the zero frequency modes do not contribute to the
   * vibrational spectrum as they describe translations and rotations. $M$ is
   * a diagonal matrix of dim 3N x 3N of atomic masses. $g _(2) ^{(x)}$ is the
   * gradient of the target state energy surface calculated at the geometry of
   * the initial state in cartesian (non-mass-weighted) coordinates
   * (${}^{(x)}$). The number ${}_(2)$ indicates that the value was calculated
   * at the second (target) state.
   *
   * $M$ is in the $-1/2$ power in both cases as the right most is a
   * transformation of the **gradient** from cartesian to mass-weighted
   * coordinates, while the leftmost is a transfromation of **coordinates**
   * from a mass-weighted vector to a mass-un-weighed vector.
   *
   * $\Delta$ as well as the cartesian displacement ($M ^{-1/2} D \Delta$) are
   * expressed in AU and are converted to Angstroms only just before
   * addition to the geometry. For this reason both matrices ($\Omega$ and
   * $M$) and the gradient vector have values stored in the atomic units.
   */

  std::cout << "State geometry is calculated within VGA."
            << std::endl;

  // arma::sqrt is an element-wise square root
  arma::Mat<double> Omega_matrix_minus2 = arma::inv(arma::square(omega_matrix));
  arma::Mat<double> mass_matrix_minus_half = arma::inv(arma::sqrt(mass_matrix));

  arma::Col<double> delta;
  delta =
      Omega_matrix_minus2 * d_matrix.t() * mass_matrix_minus_half * gradient;

  vg_calc_geom(mass_matrix_minus_half, d_matrix, delta);
  vg_calc_energy(delta, Omega_matrix_minus2);
}

/* Helper of the `ApplyKeyWords` function. Applies manual shifts and rotations
 * of the molecular geometry, which were specified in the
 * `manual_coordinates_transformation` nodes. `ground` is a `MolState` object
 * corresponig to the initial state (needed for printing geometry difference).
 */
void MolState::apply_manual_coord_transformation(const MolState &ground) {

  while (! manual_transformations.empty()) {
    arma::Col<double> shift = manual_transformations.front().first;
    arma::Col<double> rotation = manual_transformations.front().second;
    manual_transformations.pop();

    shiftCoordinates(shift);

    rotate(rotation(0) * PI, rotation(1) * PI, rotation(2) * PI);
    applyCoordinateThreshold(COORDINATE_THRESHOLD);

    std::cout << "Molecule was shifted by " << shift(0) << ", " << shift(1)
              << ", " << shift(2) << " in x, y, and z." << std::endl;

    std::cout << "Also rotated by " << rotation(0) << "*pi, " << rotation(1)
              << "*pi, and " << rotation(2) << "*pi around x, y, and z axes."
              << std::endl;

    double diff_from_ground = this->getGeomDifference(ground);
    std::cout << "Norm of the geometry difference to the initial state: "
              << sqrt(diff_from_ground) << std::endl;

    if_aligned_manually = true;
  }
}

/* Helper of the `ApplyKeyWords` function. Applies manual reodering of normal
 * modes specified in the `manual_normal_modes_reordering` node. */
void MolState::reorder_normal_modes() {
  // backup normal modes
  std::vector<NormalMode> oldNormModes(normModes);

  // copy normal modes using new order
  for (int nm = 0; nm < NNormModes(); nm++)
    getNormMode(nm) = oldNormModes[normModesOrder[nm]];

  std::cout << "Normal modes reordered." << std::endl;
}

/* Helper of the `ApplyKeyWords` function. Applies manual reodering of atoms
 * specified in the `manual_atoms_reordering` node.
 * The atoms reordergin has to be applied to both geometry and normal modes. */
void MolState::reorder_atoms() {

  // backup molecular geometry and normal modes
  std::vector<Atom> oldAtoms(atoms);
  std::vector<NormalMode> oldNormModes(normModes);

  // copy molecular geometry using the new order of atoms
  for (int a = 0; a < NAtoms(); a++)
    getAtom(a) = oldAtoms[atomsOrder[a]];

  // copy normal modes using the new order of atoms
  for (int nm = 0; nm < NNormModes(); nm++)
    for (int a = 0; a < NAtoms(); a++)
      for (int k = 0; k < CARTDIM; k++)
        getNormMode(nm).getDisplacement()(a * CARTDIM + k) =
            oldNormModes[nm].getDisplacement()(atomsOrder[a] * CARTDIM + k);

  std::cout << "Atoms reordered in both molecular geometry and normal modes."
            << std::endl;
}

/* The "initial_state" or "target_state" include have keywords that modify the
 * state. This function complese requested changes. Additionaly this function
 * applies transformations necessary for keeping `MolState` variables in line
 * with the ezFCF storage conventions. */
void MolState::ApplyKeyWords(xml_node &node_amu_table, MolState & initial_state) {

  if (!IfGradientAvailable) {
    convert_atomic_names_to_masses(node_amu_table);

    if (ifInputNMmassweighted) {
      un_mass_weight_nm();
    }
    // TODO: A sample job that presents how to use the
    // manual_coordinates_transformation is missing. A problem also for testing.
    if (!manual_transformations.empty()) {
      apply_manual_coord_transformation(initial_state);
    }
    if (if_nm_reordered_manually) {
      reorder_normal_modes();
    }
    if (if_atoms_reordered_manually) {
      reorder_atoms();
    }
    // Create diagonal matrices of atomic masses, harmonic frequencies and a
    // rectangular matrix of normal modes: these are used thorughout the program.
    create_matrices();
  }

  // Find geometry and excitation_energy using VG method
  if (IfGradientAvailable) {
    copy_data_from_the_initial_state(initial_state);
    test_vertical_gradient_dimension();
    test_for_small_frequencies_in_VGA();
    vertical_gradient_method();
  }
}
