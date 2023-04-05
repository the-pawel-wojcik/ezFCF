#ifndef _spectral_point_h_
#define _spectral_point_h_

/*! \file spectralpoint.h
\brief Spectral point I(E)
stores pair Intensity-Energy and the transition (vibrational quanta of the
initial and target states)

\ingroup MOLECULAR_PROP
*/

#include "constants.h"
#include "vibronic_state.h"

class SpectralPoint {
  double intensity;
  double initial2target_E_gap_eV;

  // energy of the hot bands (vibr. excited states of the initial electronic
  // state)
  double Epp; // E''

  // Franck-Condon factor;
  double FCF;

  bool if_print;

  VibronicState intial, target;

public:
  SpectralPoint() : if_print(true){};

  SpectralPoint(const VibronicState &init_vibst,
                const VibronicState &targ_vibst)
      : intensity(0.0), initial2target_E_gap_eV(0.0), Epp(0.0), FCF(0.0),
        if_print(true), intial(init_vibst), target(targ_vibst){};

  // Peak intensity: FCF * Boltzmann
  double getIntensity() const { return intensity; }
  void set_intensity(const double _intensity) { intensity = _intensity; }

  // Energy gap (eV) between the initial and target vibronic states.
  double get_energy() const { return initial2target_E_gap_eV; }
  void set_energy(const double _energy) { initial2target_E_gap_eV = _energy; }

  // TODO: what is E_prime_prime
  double &getE_prime_prime() { return Epp; };

  // Franck-Condon factor, |<initial|target>|^2
  double getFCF() const { return FCF; };
  void set_FCF(const double fcf) { FCF = fcf; };

  bool getIfPrint() const { return if_print; };
  void setIfPrint(const bool flag) { if_print = flag; }

  // TODO: This is a bad interface -- if that's the only way it can be
  // implemented then initial and target should be public.
  VibronicState &getVibrState1() { return intial; };
  VibronicState &getVibrState2() { return target; };

  /* Print energy, intensity, FCF and transition */
  void print(std::ostream & = std::cout) const;

  // operator<< is a a non-member friend function of SpectralPoint as it needs
  // to access the obj.inital private variable
  /* Print 0(1v13,2v12)->1(0) */
  friend std::ostream &operator<<(std::ostream &os, const SpectralPoint &obj);
};

#endif
