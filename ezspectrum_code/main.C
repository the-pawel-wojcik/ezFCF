#include "genincludes.h"
#include "genutil.h"
#include "harmonic_pes_main.h"
#include "aik_xml_parser.h"
#include <fstream>

int main(int argc, char *argv[]) {

  if (argc != 2) {
    std::cerr << argv[0]
              << ": one and only one argument required, <Job.xml>.\n";
    exit(1);
  }
  std::string arg = argv[1];

  std::cout << "Job \"" << argv[0] << ' ' << argv[1]
            << "\" has been started: " << GetTime() << '\n';

  std::ifstream xml_file(arg.c_str());
  if (!xml_file.is_open()) {
    std::cout << "Input file " << arg.c_str() << " is not found." << std::endl;
    exit(1);
  }
  xml_node node_input("input", xml_file);

  std::cout << "A copy of the \"" << argv[1] << "\" input:\n";
  std::cout << HorizontalLine << std::endl;
  node_input.print(std::cout);
  std::cout << HorizontalLine << std::endl << std::endl;

  std::ifstream xml_amu_file(ATOMIC_MASSES_FILE);
  if (!xml_amu_file.is_open()) {
    std::cout << "Atomic masses file " << ATOMIC_MASSES_FILE
              << " is not found in the current directory." 
              << "\nAttempting to use the global file:\n\t"
              << GLOBAL_ATOMIC_MASSES_FILE << std::endl;
    xml_amu_file = std::ifstream(GLOBAL_ATOMIC_MASSES_FILE);
    if (!xml_amu_file.is_open()) {
      std::cout << "Global version of the atomic masses file\n\n" 
                << ATOMIC_MASSES_FILE
                << "\n is also not found." << std::endl;
      exit(1);
    }
    else {
      std::cout << "Using the global version of the atomic masses file:\n\t"
        << GLOBAL_ATOMIC_MASSES_FILE << std::endl;
    }
  }
  else {
    std::cout << "Using the local version of the atomic masses file:\n\t"
      << ATOMIC_MASSES_FILE << std::endl;
  }
  xml_node node_amu_table("masses", xml_amu_file);

  bool done = false;
  std::string job = node_input.read_string_value("job");
  if (job == "harmonic_pes")
    done = harmonic_pes_main(arg, node_input, node_amu_table);

  if (!done)
    std::cout << "Method \"" << job << "\" is unknown, or it has failed. \n";

  std::cout << '\n'
            << "Job \"" << argv[0] << ' ' << argv[1]
            << "\" has finished: " << GetTime() << '\n';

  return EXIT_SUCCESS;
}
