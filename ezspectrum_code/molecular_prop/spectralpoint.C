#include "spectralpoint.h"

void SpectralPoint::print(std::ostream & os) const
{
  os << std::fixed << std::setprecision(4) << std::setw(7) << -initial2target_E_gap_eV;
  os << std::scientific << std::setprecision(6) << std::setw(18) << intensity << "  ";
  os << std::scientific << std::setprecision(6) << std::setw(11) << std::showpos << FCF << std::noshowpos ;
  os << std::fixed << std::setprecision(3) << std::setw(10) << Epp/KELVINS2EV << "  ";
  os << *this << std::endl;
}

// It's a a non-member friend function of SpectralPoint as it needs to access
// the obj.inital private variable
std::ostream &operator<<(std::ostream &os, const SpectralPoint &obj) {
  os << obj.intial << "->" << obj.target;
  return os;
}
