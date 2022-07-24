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
//FIXIT it: replace this class by kmatrix, which allows vector computations.
//Vector would be Kmatrix[N,1] or [1,N] and can be used in all operations
#include <vector>
#include "vector3d.h"
#include "kmatrix.h"
#include <algorithm>
 

class MolState
{
  //! N Atoms 
  std::vector<Atom> atoms;
  //! M Normal modes and frequencies (3N-5 or 3N-6) Stored in mass unweighted format in Angstoms (i.e as in ACESII)
  std::vector<NormalMode> normModes;
  //! Normal modes order (relative to the input file's order)
  std::vector<int> normModesOrder;
  //! Gradient calculated in the caresian (non-mass-weighted) coordinates, a 3N x 1 matrix TODO: hande the cases where gradient is available only in the mass-weighted coordinates
  KMatrix gradient;
  //! may be removed later
  bool ifLinear;
  //! IP
  double energy;
  
  //!move this to functions...
  bool ifLetterOrNumber(char Ch);

  Vector3D centerOfMass;
  KMatrix momentOfInertiaTensor;

  //! reduced masses
  KMatrix reduced_masses;
  
  //!if geometry transformation was performed manually
  bool if_aligned_manually;
  bool if_nm_reordered_manually;

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
  bool getNormalModeOverlapWithOtherState(MolState& other, KMatrix& overlap, std::vector<int>& normal_modes_list);

  //--- interface ---------------------------------------------------

  //! Returns Atom #i (this atoms can be different/isotops from the normal mode section)
  Atom& getAtom(int i){ return atoms[i]; } 
  //! Returns NormalMode #i (AIK: Are we using mass-weighted NModes or not? FIXIT: need to figure out and comment)
  NormalMode& getNormMode(int i) { return normModes[i]; } 
  //! Returns NormalMode order index
  int getNormModeIndex(int i) { return normModesOrder[i]; } 
  //! Linear?
  bool IfLinear () const { return ifLinear; }
  //! Energy (IP)
  double Energy () const { return energy; }
  //! Returns number of atoms
  int NAtoms() const { return atoms.size(); }
  //! Returns number of normal modes
  int NNormModes() const { return normModes.size(); }
  //! returns True if vertical gradient was used
  //FIXIT: this is bad, need to rewrite in a safer way
  bool IfGradient() const { return static_cast<bool>(gradient.Size()); }

  //--- alignment ------------- ------------------------------------
  //! align each state: center of mass in the coordinates origin, moment of ineretia principal axes along the coordinate axes:
  void align();
  //! align the state with the "other" state by rotating around every axes (x,y,z) by pi/2;
  void align(MolState& other);
  //! check that states are similar (same number of atoms and same "liniarity", i.e. same number of normal modes); 
  bool ifSimilar(MolState& other);

  //--- Geometry transformations ------------------------------------
  // Center eof mass vector:
  Vector3D& getCenterOfMass();
  // Moment of inertia tensor:
  KMatrix& getMomentOfInertiaTensor();
  // Shifts coordinate's origin by "vector" 
  void shiftCoordinates(Vector3D& vector);
  // matrix multiolication coordinates*matrix_3x3; "coordinates" matrix is of 3 columns: x,y,z. matrix_3x3 has eigen vectors of the transformation in rows.
  void transformCoordinates(const KMatrix& matrix_3x3);
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

