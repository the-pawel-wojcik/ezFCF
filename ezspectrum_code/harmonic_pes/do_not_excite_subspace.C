#include "do_not_excite_subspace.h"

/* ```node``` points to an input node that can contain
 * "do_not_excite_subspace" as its subnode. If no
 * "do_not_excite_subspace" subnode is found, the object still works: use
 * the member function ```empty()``` to check the status of the object.
 *
 * ```n_mol_nms``` is the number of molecular normal modes, i.e.,
 * 3 (# atoms) - 5 or 6. It is needed for input tests and in the
 * ```get_active_subspace()``` function. */
DoNotExcite::DoNotExcite(const xml_node &node, const int n_mol_nms)
    : size(0), n_molecular_nms(n_mol_nms) {

  // TODO: Unify this test: see JobParameters constructor
  int number_of_inputs = node.find_subnode("do_not_excite_subspace");

  if (number_of_inputs == 0) {
    return;
  }

  xml_node node_input(node, "do_not_excite_subspace", 0);

  size = node_input.read_int_value("size");
  if ((size <= 0) or (size >= n_mol_nms)) {
    error("\"do_not_excite_subspace\"->\"size\" is out of range.");
  }
  std::string normal_modes_input = node_input.read_string_value("normal_modes");
  std::stringstream nmodes_stream(normal_modes_input);
  parse_normal_modes_stream(nmodes_stream);
  run_tests();
}

/* The "normal_modes" string is a list of integers. Each integer is a normal
 * number of normal mode that has to be excluded from later calculations.
 * The "normal_modes" string is put into the `nmodes_stream` as an
 * istringstream. */
void DoNotExcite::parse_normal_modes_stream(std::stringstream &nmodes_stream) {

  int mode_no;
  for (int nm = 0; nm < size; nm++) {
    nmodes_stream >> mode_no;

    // read check
    if (nmodes_stream.fail()) {
      std::stringstream msg;
      msg << "Format error in "
          << "\"do_not_excite_subspace\"->\"normal_modes\" \n"
          << " (should be a list of whitespaces-separated integers).";
      error(msg);
    }

    subspace.insert(mode_no);
  }
}

/* Prints summary and warrnings about the parsed "do_not_excite_subspace". */
void DoNotExcite::print_summary(const bool nm_were_reordered) const {
  if (empty()) {
    std::cout << "All normal modes modes will be excited\n\n";
    return;
  }

  if (nm_were_reordered)
    std::cout
        << "WARNING! Normal modes of the target state were reordered!\n"
        << "         New order is used for the \"do_not_excite_subspace\"."
        << "\n\n";

  std::cout << "The following normal modes will have no vibrational "
               "excitations:\n";

  for (int mode_no : subspace)
    std::cout << mode_no << ' ';
  std::cout << "\n\n";

  std::cout << "The following normal modes will be excited:\n";
            /* << " (for both states the order is same as in the input):\n"; */
  for (int mode_no : get_active_subspace())
    std::cout << mode_no << ' ';
  std::cout << "\n\n";
}

/* Make sure that the list of modes to be excluded makes sense. */
void DoNotExcite::run_tests() const {
  
  // Tests assume that the do_not_excite_subspace is NOT empty.
  if (empty()) {
    return;
  }

  check_the_largest(n_molecular_nms);
  check_the_smallest();
  check_for_duplicates();
}

/* Helper of `DoNotExcite::run_tests`. Make sure that the largest normal mode in
 * the "do_not_excite_subspace" is smaller than number of molecular normal
 * modes. */
void DoNotExcite::check_the_largest(const int n_mol_nms) const {
  auto largest_mode = std::max_element(subspace.begin(), subspace.end());
  int largest_exclude_nm_no = *largest_mode;

  if (largest_exclude_nm_no < n_mol_nms)
    // all good
    return;

  std::stringstream msg;
  msg << "A normal mode number " << largest_exclude_nm_no
      << " in the \"do_not_excite_subspace\" is too large.\n"
      << " This molecule has " << n_mol_nms << " normal modes.\n"
      << " Numbering of normal modes starts at 0.";
  error(msg);
}

/* Helper of `DoNotExcite::run_tests`. Make sure numbers of modes to be excluded
 * are not negative. */
void DoNotExcite::check_the_smallest() const {
  auto smallest_mode = std::min_element(subspace.begin(), subspace.end());
  int smallest_mode_number = *smallest_mode;

  if (smallest_mode_number >= 0)
    // all good
    return;

  std::stringstream msg;
  msg << "The normal mode number " << smallest_mode_number
      << " in the \"do_not_excite_subspace\" is negative.\n"
      << " Numbering of normal modes starts at 0.";
  error(msg);
}

/* Helper of `DoNotExcite::run_tests`. Assure that the user did not input
 * duplicates in the list of modes to be excluded. */
void DoNotExcite::check_for_duplicates() const {
  int subspace_size = subspace.size();

  if (subspace_size == size)
    return;

  std::stringstream msg;
  msg << "Duplicates dected in "
      << "\"do_not_excite_subspace\"->\"normal_modes\".\n"
      << " Input contains " << subspace_size << " unique elements:\n";
  for (auto i : subspace) {
    msg << " " << i;
  }
  msg << "\n";
  error(msg);
}

/* Create the "excite subspace": full_space - do_not_excite_subspace 
 * Only normal modes from this subspace will be excited */
std::vector<int> DoNotExcite::get_active_subspace() const {
  std::vector<int> active_nms;
  for (int nm = 0; nm < n_molecular_nms; nm++) {
    auto intSet_iter = subspace.find(nm);
    if (intSet_iter == subspace.end())
      active_nms.push_back(nm);
  }
  return active_nms;
}
