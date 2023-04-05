#include "spectrum.h"
#include <fstream>

void Spectrum::PrintStickTable() const {
  if (getNSpectralPoints() > 0) {
    std::cout
        << "Energy,eV    Intensity     FC factor         E\",K   Transition"
        << "\n\n";
    for (const auto &spectral_point : spectralPoints) {
      if (spectral_point.getIfPrint()) {
        spectral_point.print(std::cout);
      }
    }
  } else {
    std::cout << "\n\n\n"
              << "WARNING! The spectrum is empty.\n\n"
              << "         Please refer to \"My spectrum is empty!\" in the\n"
              << "         \"Common problems\" section of the manual.\n\n\n\n";
  }
  std::cout << HorizontalLine << std::endl;
}

void Spectrum::PrintStickTable(const std::string spectrumFileName) const {
  std::ofstream spectrumF(spectrumFileName, std::ios::out);
  for (const auto &spectral_point : spectralPoints) {
    if (spectral_point.getIfPrint()) {
      spectral_point.print(spectrumF);
    }
  }
  spectrumF.close();
}

void Spectrum::AddSpectralPoint(const SpectralPoint &PointToAdd) {
  spectralPoints.push_back(PointToAdd);
};

void Spectrum::AddSpectralPoint(const double E, const double I,
                                const double fcf, const double Epp,
                                VibronicState state_ini,
                                VibronicState state_targ,
                                const bool if_print_flag) {
  SpectralPoint tmpPoint;
  for (int nm = 0; nm < state_ini.getVibrQuantaSize(); nm++)
    tmpPoint.getVibrState1().addVibrQuanta(state_ini.getV()[nm],
                                           state_ini.getEx()[nm]);
  for (int nm = 0; nm < state_targ.getVibrQuantaSize(); nm++)
    tmpPoint.getVibrState2().addVibrQuanta(state_targ.getV()[nm],
                                           state_targ.getEx()[nm]);
  tmpPoint.getVibrState1().setElStateIndex(state_ini.getElStateIndex());
  tmpPoint.getVibrState2().setElStateIndex(state_targ.getElStateIndex());

  tmpPoint.set_energy(E);
  tmpPoint.set_intensity(I);
  tmpPoint.set_FCF(fcf);
  tmpPoint.getE_prime_prime() = Epp;

  tmpPoint.setIfPrint(if_print_flag);

  AddSpectralPoint(tmpPoint);
}
