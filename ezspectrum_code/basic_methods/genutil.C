/*! \file genutil.C
\ingroup BASIC_FUNCTIONS
*/

#include "genutil.h"
#include "constants.h"
#include "genincludes.h"

// Return string with the current time:
std::string GetTime(){
  time_t rawtime;
  time ( &rawtime );
  return  ctime ( &rawtime );
}

time_t GetRawTime(){
  time_t rawtime;
  time ( &rawtime );
  return  rawtime;
}

void error(const std::string & msg)
{
  std::cout << "\nezFCF: Error!\n" << msg <<"\n\n";
  exit(EXIT_FAILURE);
}

/* Sometimes a string stream is more convenient */
void error(const std::stringstream & msg)
{
  error(msg.str());
}

// Prints a vector using a compact q-chem's format for gradient
void print_qchem_style_vector(arma::Col<double> v, std::string header)
{
  std::cout << header << std::endl;

  if (v.n_rows % CARTDIM != 0) {
    std::cout << "The vector length is not a mupltiple of 3.\n";
    return;
  }
  const int n_atoms = v.n_rows / CARTDIM;

  std::vector<std::string> num_to_xyz = {std::string("x"), std::string("y"),
                                         std::string("z")};

  // number of atoms per print block
  const int blk_size = 4;

  for (int block_idx = 0; block_idx < n_atoms/blk_size; ++block_idx) {
    // print header with atom numbers
    std::cout << std::string(4, ' ');
    for (int at_idx = block_idx * blk_size; at_idx < (block_idx + 1) * blk_size; ++at_idx) {
      std::cout << " " << std::setw(6) << at_idx << std::string(3, ' ');
    }
    std::cout << std::endl;

    // print the vector data
    for (int xyz = 0; xyz < 3; ++xyz) {
      // 4 characters. Coordinate labels
      std::cout << "  " << num_to_xyz[xyz] << " ";
      for (int at_no = 0; at_no < blk_size; ++at_no) {
        std::cout << " " << std::setw(9) << std::fixed << std::setprecision(5)
                  << v[(block_idx * blk_size + at_no * 3) + xyz];
      }
      // end of every line
      std::cout << std::endl;
    }
  }

  const int blk_leftover = n_atoms % blk_size;
  if (blk_leftover != 0) {
    // print the leftover
    // Copy of the first for loop content
    // TODO: this should be wrapped back up there

    const int n_prntd_atoms = n_atoms - n_atoms % blk_size;

    std::cout << std::string(4, ' ');
    for (int at_idx = n_prntd_atoms; at_idx < n_atoms; ++at_idx) {
      std::cout << " " << std::setw(6) << at_idx << std::string(3, ' ');
    }
    std::cout << std::endl;

    // print the vector data
    for (int xyz = 0; xyz < 3; ++xyz) {
      // 4 characters. Coordinate labels
      std::cout << "  " << num_to_xyz[xyz] << " ";
      for (int at_idx = n_prntd_atoms; at_idx < n_atoms; ++at_idx) {
        std::cout << " " << std::setw(9) << std::fixed << std::setprecision(5)
                  << v[at_idx * 3 + xyz];
      }
      // end of every line
      std::cout << std::endl;
    }
  }

  std::cout << std::endl;
}
