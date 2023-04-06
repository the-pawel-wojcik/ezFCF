#include "single_excitations.h"

SingleExcitations::SingleExcitations(xml_node &node_head,
                                     const MolState &targ_elst,
                                     const int n_molecular_nms,
                                     const int iniN) {
  size_t n_single_ex = node_head.find_subnode("single_excitation");
  if (n_single_ex == 0) {
    return;
  }

  targ_elst.warn_about_nm_reordering("single excitations");

  std::cout << "List of single excitations added to the spectrum:\n"
            << std::flush;

  for (int ex_idx = 0; ex_idx < n_single_ex; ++ex_idx) {
    xml_node node_single_ex(node_head, "single_excitation", ex_idx);

    std::string init_str = node_single_ex.read_string_value("ini");
    std::string targ_str = node_single_ex.read_string_value("targ");
    VibronicState init_vibronic_st(init_str, n_molecular_nms, iniN);
    VibronicState targ_vibronic_st(targ_str, n_molecular_nms, iniN);
    single_excitations.push_back(
        SpectralPoint(init_vibronic_st, targ_vibronic_st));
  }
}
