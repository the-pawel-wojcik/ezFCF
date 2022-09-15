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
#include <algorithm>


class MolState
{
  //! N Atoms 
  std::vector<Atom> atoms;
  //! M Normal modes and frequencies (3N-5 or 3N-6)
  //! Stored in the mass unweighted format in Angstoms (i.e as in ACESII)
  std::vector<NormalMode> normModes;
  //! Normal modes order (relative to the input file's order)
  std::vector<int> normModesOrder;
  //! number of molecular normal modes: 3N-6 or 3N-5 for linear
  int n_molecular_nm;
  //! Gradient calculated in the caresian (non-mass-weighted) coordinates a 3N vector
  //! TODO: handle the cases where gradient is available only in the mass-weighted coordinates
  arma::Col<double> gradient;
  //! may be removed later
  bool ifLinear;
  //! excitation energy (formerly IP)
  double energy;
  //! calculate the state's properties using the vertical gradient method
  bool IfGradientAvailable;

  //!move this to functions...
  bool ifLetterOrNumber(char Ch);

  arma::Col<double> centerOfMass;
  arma::Mat<double> momentOfInertiaTensor;

  //! reduced masses
  arma::Col<double> reduced_masses;

  //!if geometry transformation was performed manually
  bool if_aligned_manually;
  bool if_nm_reordered_manually;

  // ==  Helpers of the MolState::Read function ==
  
  void Read_excitation_energy(xml_node &node_state);
  void Read_molecular_geometry(xml_node &node_state);
  void Read_normal_modes(xml_node &node_state);

  // ==  Helpers of the MolState::Read function ==
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

