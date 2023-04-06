#ifndef _spectrum_h_
#define _spectrum_h_

/*! \file spectrum.h
\brief Spectrum object
Implemented:
1) stores stick spectrum (vector of "SpectralPoint"s {I,E}) )
2) sorts it by energy
3) prints as a table

\ingroup MOLECULAR_PROP
*/
#include "genincludes.h"
#include "genutil.h"
#include "spectralpoint.h"

class Spectrum {
  std::vector<SpectralPoint> spectralPoints;

  struct SortByEnergy {
    bool operator()(const SpectralPoint &p1, const SpectralPoint &p2) const {
      return p1.get_energy() < p2.get_energy();
    }
  };

public:
  // Both getSpectralPoint and getSpectralPoints are bad practice -- might 
  // as well call this class a struct and keep spectralPoints a public 
  // member.
  SpectralPoint &getSpectralPoint(int i) { return spectralPoints[i]; };
  std::vector<SpectralPoint> &getSpectralPoints() { return spectralPoints; };

  int getNSpectralPoints() const { return spectralPoints.size(); };

  void AddSpectralPoint(const SpectralPoint &PointToAdd);
  void AddSpectralPoint(const double E, const double I, const double fcf,
                        const double Epp, VibronicState state_ini,
                        VibronicState state_targ,
                        const bool if_print_flag = true);

  void Sort() {
    sort(spectralPoints.begin(), spectralPoints.end(), SortByEnergy());
  };

  /* Print with comments to std::cout. */
  void PrintStickTable() const;

  /* Print without decorations to a file. */
  void PrintStickTable(const std::string spectrumFileName) const;

  void clear() { spectralPoints.clear(); }
};

#endif
