/*! \file normalmode.C
  \brief Normal mode: frequency and 3N displacements (N-number of atoms)
  \ingroup DATA_CLASSES
  */

#include "normalmode.h"

NormalMode::NormalMode(const NormalMode& other) : displacement(other.displacement),
  freq(other.freq),
  nAtoms(other.nAtoms) {} 

NormalMode& NormalMode::operator=(const NormalMode& other)
{
  if (this != &other)
  {
    freq = other.freq;
    displacement = other.displacement;
    nAtoms = other.nAtoms;
  }
  return *this;
}

// does matrix_3x3 * (x y z) ^T to every atoms' coordinates
// matrix_3x3 has eigen vectors of the transformation in rows.
void NormalMode::transformCoordinates(const arma::Mat<double>& matrix_3x3)
{
  arma::Col<double> old_displacement = displacement;
  displacement.fill(0.0);

  for (int a=0; a<nAtoms; a++)
    for (int i=0; i<CARTDIM; i++) // columns #, i.e. i=0 is x, i=1 is y, i=2 is z;
      for (int j=0; j<CARTDIM; j++) // rows #, i.e. eigen vectors, wich are in rows;
        displacement(a*CARTDIM + i) += old_displacement(a*CARTDIM + j) * matrix_3x3(i, j);
}

// set zero for every displacement element smaller than threshold
void NormalMode::applyCoordinateThreshold(const double threshold)
{
  for (int a=0; a<nAtoms; a++)
    for (int j=0; j<CARTDIM; j++ )
      if (fabs(displacement(a*CARTDIM+j))< threshold)
        displacement(a*CARTDIM+j)=0.0;
}

void NormalMode::rotateX_90deg() {
  double coord_tmp;
  // for each atom take (x,y,z) vector and transform it as in Atom class
  for (int a = 0; a < nAtoms; a++) {
    // x'=x; y'=z; z'=-y;
    coord_tmp=-displacement(a*CARTDIM+1);// =-y
    displacement(a*CARTDIM+1)=displacement(a*CARTDIM+2);  // y'=z
    displacement(a*CARTDIM+2)=coord_tmp; // z'=-y
  }
}

void NormalMode::rotateY_90deg()
{
  double coord_tmp;
  for (int a=0; a<nAtoms; a++)// for each atom take (x,y,z) vector and transform it as in Atom class
  {
    // x'=-z; y'=y; z'=x;
    coord_tmp=displacement(a*CARTDIM+0);// =x
    displacement(a*CARTDIM+0)=-displacement(a*CARTDIM+2);  // x'=-z
    displacement(a*CARTDIM+2)=coord_tmp; // z'=x
  }
}

void NormalMode::rotateZ_90deg()
{
  double coord_tmp;
  for (int a=0; a<nAtoms; a++)// for each atom take (x,y,z) vector and transform it as in Atom class
  {
    // x'=y; y'=-x; z'=z;
    coord_tmp=-displacement(a*CARTDIM+0);// =-x
    displacement(a*CARTDIM+0)=displacement(a*CARTDIM+1);  // x'=y
    displacement(a*CARTDIM+1)=coord_tmp; // y'=-x
  }
}


