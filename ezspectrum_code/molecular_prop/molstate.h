#ifndef _molstate_h_
#define _molstate_h_

/*! \file molstate.h
  \brief Molecular state: Stores Geometry, Normal Modes & Frequencies. 
  Also reads input data from the XML file (vm'06)
  20-12-30: Replace parser by aik_xml_parser 
  \ingroup MOLECULAR_PROP
  */

#include "genincludes.h"
#include "aik_xml_parser.h"
#include "atom.h"
#include "normalmode.h"
#include <vector>
#include <queue>
#include <set>
#include <algorithm>


class MolState
{
  //! N Atoms: these atoms come from parsing the "geometry" node
  std::vector<Atom> atoms;
  //! N Atoms: these atoms come from parsing the "normal_modes" node.
  //! HINT: `atoms` and `nm_atoms` can come in different order
  //! TODO: This can cause troubles down the way:
  //! -- it's difficult to tell when this is the case as the same atoms
  //! can be reshuffled
  //! -- If the order of cartesian coordinates in normal modes doesn't match 
  //!    the order in which they appear in the geometries the restuls would be 
  //!    garbage.
  std::vector<Atom> nm_atoms;
  //! n_molecular_nm Normal modes and frequencies (3N-5 or 3N-6)
  //! Stored in the mass unweighted format in Angstoms (i.e as in ACESII)
  std::vector<NormalMode> normModes;
  //! The xml input may store the normal modes using mass weighted convention
  //! The MolState::Read function reads them as they are 
  //! only in MolState::Transform the nodes get mass-un-weighted
  bool ifInputNMmassweighted;
  //! number of molecular normal modes: 3N-6 or 3N-5 for linear
  int n_molecular_nm;
  //! may be removed later
  bool ifLinear;
  //! excitation energy (formerly IP), the adiabatic energy gap to the initial state
  double energy;
  //! Gradient calculated in the caresian (non-mass-weighted) coordinates, a 3N vector
  arma::Col<double> gradient;
  //! calculate the state's properties using the vertical gradient method
  bool IfGradientAvailable;

  //! Three matrices below are initilized in create_matrices function
  //! mass matrix uses the order of atoms from MolState::atoms
  arma::Mat<double> mass_matrix;
  //! matrix of harmonic frequencies
  arma::Mat<double> omega_matrix;
  //! matrix that stores normal modes as its columns  
  arma::Mat<double> d_matrix;


  arma::Col<double> centerOfMass;
  arma::Mat<double> momentOfInertiaTensor;

  //! reduced masses
  arma::Col<double> reduced_masses;


  // == Variables for manual tweaks ==

  //! if geometry transformation was performed manually
  bool if_aligned_manually;
  //! stores manual coordinate transformations
  std::queue<std::pair<arma::Col<double>, arma::Col<double>>>
      manual_transformations;

  bool if_nm_reordered_manually;
  //! Normal modes order (relative to the input file's order)
  std::vector<int> normModesOrder;

  bool if_atoms_reordered_manually;
  //! Atoms order (relative to the input file's order)
  std::vector<int> atomsOrder;

  // ==  Helpers of the MolState::Read function ==
  void Read_excitation_energy(xml_node &);
  void Read_molecular_geometry(xml_node &);
  void Read_normal_modes(xml_node &);
  void Read_frequencies(xml_node &);
  void Read_vertical_gradient(xml_node &);
  void Read_manual_coord_transformations(xml_node &);
  void Read_normal_modes_reorder(xml_node &);
  void Read_atoms_reorder(xml_node &);

  // -- helpers to the MolState::Read_normal_modes function --
  void Read_normal_modes_atoms(std::string &);

  // ==  Helpers of the MolState::Transform function ==
  void convert_atomic_names_to_masses(xml_node &);
  void un_mass_weight_nm();
  void create_mass_matrix();
  void create_matrices();
  void calculate_vertical_gradient_geometry();
  void apply_manual_coord_transformation();
  void reorder_normal_modes();
  void reorder_atoms();

  bool ifLetterOrNumber(char Ch);

public:
  MolState ();
  MolState (const MolState& other);
  //  ~MolState();
  MolState& operator=(const MolState& other);
  // TODO: the rule of three

  //! Read the state data from file
  bool Read(xml_node& node_state, xml_node& node_amu_table);
  void Print();
  void printGeometry(); 
  void printNormalModes(); 
  void printGradient(); 
  // returns true if the overlap matrix is diagonal; makes a list of normal modes which form a non-diagonal minor.
  bool getNormalModeOverlapWithOtherState(MolState& other, arma::Mat<double>& overlap, std::vector<int>& normal_modes_list);

  //--- interface ---------------------------------------------------

  //! Returns Atom #i (this atoms can be different/isotops from the normal mode section)
  Atom& getAtom(int i){ return atoms[i]; } 
  //! Returns NormalMode #i 
  //! Stored in the mass unweighted format in Angstoms (i.e as in ACESII)
  NormalMode& getNormMode(int i) { return normModes[i]; } 
  //! Returns NormalMode order index
  int getNormModeIndex(int i) { return normModesOrder[i]; } 
  //! Linear?
  bool IfLinear () const { return ifLinear; }
  //! Excitation Energy (formerly IP)
  double Energy () const { return energy; }
  //! Returns number of atoms
  int NAtoms() const { return atoms.size(); }
  //! Returns number of molecular normal modes, i.e., 3N - 5/6
  int NNormModes() const { return normModes.size(); }
  //! returns True if vertical gradient was used
  bool IfGradient() const { return IfGradientAvailable; }

  //--- alignment ------------- ------------------------------------
  //! align each state: center of mass in the coordinates origin, moment of ineretia principal axes along the coordinate axes:
  void align();
  //! align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  void align(MolState& other);
  //! check that states are similar (same number of atoms and same "liniarity", i.e. same number of normal modes); 
  bool ifSimilar(MolState& other);

  //--- Geometry transformations ------------------------------------
  // Center eof mass vector:
  arma::Col<double>& getCenterOfMass();
  // Moment of inertia tensor:
  arma::Mat<double>& getMomentOfInertiaTensor();
  // Shifts coordinate's origin by "vector" 
  void shiftCoordinates(arma::Col<double>& vector);
  // matrix multiplication coordinates*matrix_3x3; 
  // "coordinates" matrix is of 3 columns: x,y,z. 
  // matrix_3x3 has eigen vectors of the transformation in rows.
  void transformCoordinates(const arma::Mat<double>& matrix_3x3);
  void applyCoordinateThreshold(const double threshold);
  double getGeomDifference(MolState& other);
  // three rotations around x, y, z by PI/2:
  void rotateX_90deg();
  void rotateY_90deg();
  void rotateZ_90deg();
  // arbitrary rotation around x ,y, z axes by angles alpha_x, alpha_y,alpha_z respectevly:
  void rotate(const double alpha_x,const double alpha_y,const double alpha_z);
  bool ifAlignedManually();
  bool ifNMReorderedManually();

};


#endif

