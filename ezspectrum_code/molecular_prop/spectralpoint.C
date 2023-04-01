#include "spectralpoint.h"

/* Print energy, intensity, FCF and transition */
void SpectralPoint::print() const
{
  std::cout << std::fixed << std::setprecision(4) << std::setw(7) << -E;
  std::cout << std::scientific << std::setprecision(6) << std::setw(18) << I << "  ";
  std::cout << std::scientific << std::setprecision(6) << std::setw(11) << std::showpos << FCF << std::noshowpos ;
  std::cout << std::fixed << std::setprecision(3) << std::setw(10) << Epp/KELVINS2EV << "  ";
  std::cout << *this << std::endl;
}

/* Print 0(1v13,2v12)->1(0) */
std::ostream &operator<<(std::ostream &os, const SpectralPoint &obj) {
  obj.intial.print(os);
  os << "->";
  obj.target.print(os);
  return os;
}
