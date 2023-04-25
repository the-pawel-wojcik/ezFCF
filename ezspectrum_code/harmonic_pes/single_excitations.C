#include "single_excitations.h"

SingleExcitations::SingleExcitations(xml_node &node_head,
                                     const int n_molecular_nms, const int iniN,
                                     const int targN) {

  size_t n_single_ex = node_head.find_subnode("single_excitation");
  if (n_single_ex == 0) {
    return;
  }

  for (int ex_idx = 0; ex_idx < n_single_ex; ++ex_idx) {
    xml_node node_single_ex(node_head, "single_excitation", ex_idx);

    std::string init_str = node_single_ex.read_string_value("ini");
    VibronicState init_vibronic_st(init_str, n_molecular_nms, iniN);

    std::string targ_str = node_single_ex.read_string_value("targ");
    VibronicState targ_vibronic_st(targ_str, n_molecular_nms, targN);

    single_excitations.push_back(
        SpectralPoint(init_vibronic_st, targ_vibronic_st));
  }
}
