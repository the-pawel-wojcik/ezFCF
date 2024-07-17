/*! \file genutil.C
\ingroup BASIC_FUNCTIONS
*/

#include "genutil.h"
#include "constants.h"
#include "genincludes.h"

void get_atomic_masses_file(std::ifstream & xml_amu_file) {

  xml_amu_file.open(ATOMIC_MASSES_FILE);
  if (xml_amu_file.is_open()) {
    std::cout << "Using the local version of the atomic masses file:\n\t"
              << ATOMIC_MASSES_FILE << std::endl;
    return;
  }

  std::cout << "Atomic masses file " << ATOMIC_MASSES_FILE
    << " not found in the current directory.\n";

  // See if enviroment points to the data directory
  std::string data_path;
  char *val = getenv(ENVIRONMENT_VAR_NAME.c_str());
  if (val == NULL) {
    std::cout << "ezFCF's environmental variable \"" << ENVIRONMENT_VAR_NAME
      << "\" is not set.\n";
    data_path = GLOBAL_DATA_PATH;
  } else {
    std::cout << "ezFCF's root directory was set using the \""
      << ENVIRONMENT_VAR_NAME << "\" environmental variable.\n";
    data_path = std::string(val) + std::string("/share/ezFCF");
  }
  std::cout << "Looking for the atomic masses file in:\n\t" << data_path
    << "\n";

  const std::string atomic_masses_file =
    data_path + std::string("/") + ATOMIC_MASSES_FILE;
  xml_amu_file = std::ifstream(atomic_masses_file);
  if (xml_amu_file.is_open()) {
    std::cout << "Using the following atomic masses file:\n\t"
      << atomic_masses_file << std::endl;
    return;
  }

  std::cout
    << "Cannot open the atomic masses file\n\t" 
    << atomic_masses_file
    << "\nError: Make sure that the atomic masses file exists and is readable." 
    << std::endl;
  exit(1);
}


double Boltzmann_factor(double temperature, double energy) {
  if (temperature == 0 and energy == 0)
    return 1.0; // intensity unchanged

  if (temperature == 0 and energy > 0)
    return 0.0; // suppress hot bands at T = 0

  double exponent = energy / (temperature * KELVINS2EV);
  if (exponent > 100)
    return 0.0; // original code wanted the value to be > 10^-44 = exp(-100)

  return exp(-exponent);
}

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
