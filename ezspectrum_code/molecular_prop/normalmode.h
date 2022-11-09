#ifndef _normalmode_h_
#define _normalmode_h_

/*! \file normalmode.h
  \brief Normal mode: frequency and 3N displacements (N-number of atoms)
  \ingroup DATA_CLASSES
  */

#include "constants.h"
#include "genincludes.h"

// Stored in mass-unweighted form
class NormalMode {
  // 1D structure AtomNumber*CARTDIM+CoordinateNumber:
  arma::Col<double> displacement;
  double freq;
  int nAtoms;

public:
  NormalMode(int n, double fr)
      : displacement(n * CARTDIM, arma::fill::zeros), freq(fr), nAtoms(n){};
  NormalMode(const NormalMode &other);
  NormalMode& operator=(const NormalMode &other);
  inline bool operator<(const NormalMode &other) const {
    return this->getFreq() < other.getFreq();
  }

  //! Returns frequency
  double &getFreq() { return freq; }
  const double &getFreq() const { return freq; }
  // 3N displacements (x,y,z for N atoms)
  arma::Col<double> &getDisplacement() { return displacement; }

  void transformCoordinates(const arma::Mat<double> &matrix_3x3);
  void applyCoordinateThreshold(const double threshold);
  // three rotations around x, y, z by PI/2:
  void rotateX_90deg();
  void rotateY_90deg();
  void rotateZ_90deg();
};
#endif
