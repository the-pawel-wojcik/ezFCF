#ifndef _vibronic_state_h_
#define _vibronic_state_h_

/*! \file vibronic_state.h
\brief Vibronic state,
i.e. stores electronic state index
and excitation quantum number along each normal mode;

\ingroup MOLECULAR_PROP
*/

#include "aik_xml_parser.h"
#include "genincludes.h"

// TODO: Store the excitation values as a list or a vector of pairs
// the pairs will likely need its own class as a container.
//
// TODO: Fix the use of this class. Even though this class has the
// tempting idea of storing only the nodes that have non-zero excitations
// at the end it is still used as if it was a vector<int> where
// the index is the normal_mode_number and the value is the excitation
// with the excitations being set to either zero or -1 by default.
// This should be changed by inspecting the palces where the objects
// of this class are used.
class VibronicState {
  //! electronic state index: 0=initial state, 1,2... target states
  int elStateIndex;

  // TODO: vibrQuanta and excite_subspace are tied together.
  // std::map<int, int> would keep them safely together.

  //! Vibrational state quantum numbers, e.g. (0,2,1) for three normal modes
  std::vector<int> vibrQuanta;
  //! numbers of normal modes; e.g. if the "excite subspace" is {v3,v7,v12},
  //! then (0,2,1) excitaion will be (2v7, 1v12) in the full space
  std::vector<int> excite_subspace;

public:
  VibronicState() : elStateIndex(0){};

  // A ctor that accepts a string "3v4 5v1 1v14" as an input
  VibronicState(std::string &input, const int n_molecular_normal_modes,
                const int el_st_idx);

  //! remove all vibronic excitations and set the electronic state index to 0
  void reset() {
    elStateIndex = 0;
    vibrQuanta.assign(vibrQuanta.size(), 0);
  };

  //! empty all containers and set the electronic state index to 0
  void clear() {
    elStateIndex = 0;
    vibrQuanta.clear();
    excite_subspace.clear();
  };

  int getElStateIndex() const { return elStateIndex; };
  void setElStateIndex(const int index) { elStateIndex = index; };

  //! returns the number of normal modes with specified occupation
  int getNNormModes() const { return vibrQuanta.size(); };
  int getVibrQuantaSize() const { return vibrQuanta.size(); };

  // adds a quanta the normal mode number
  void addVibrQuanta(const int quantNumber, const int nm) {
    vibrQuanta.push_back(quantNumber);
    excite_subspace.push_back(nm);
  };
  int getVibrQuanta(const int which_mode) { return vibrQuanta[which_mode]; };
  void setVibrQuanta(const int which_mode, const int how_excited) {
    vibrQuanta[which_mode] = how_excited;
  };

  //! checks if nm is in the "excite subspace":
  //! if it is in there, return number of excitations
  //! otherwise returns zero
  int getV_full_dim(const int nm);

  //! returns the number of the normal mode nm (nm in the full space) in the
  //! "excite subspace"
  // int getEx_full_dim(const int nm);

  int getTotalQuantaCount();

  void incrIndex(const int index, const int add) { vibrQuanta[index] += add; };

  // avoid this:
  std::vector<int> &getV() { return vibrQuanta; };
  std::vector<int> &getEx() { return excite_subspace; };

  bool if_equal(VibronicState &other);

  // Note: operator<< offers the same functionality.
  //
  // print the state in format N(avk,bvl,cvm,...)
  // N is the electronic state index,
  // a, b, c, ... -- how how excited
  // k, l, m, ... -- normal mode number
  // 'v' is a character that separates how excited from mode number
  void print_to(std::ostream &os = std::cout) const;
};

// an alternative way of accessing print(ostream&)
std::ostream &operator<<(std::ostream &, const VibronicState &);

#endif
