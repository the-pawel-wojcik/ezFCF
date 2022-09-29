#include "genincludes.h"
#include "genutil.h"
#include "harmonic_pes_main.h"
#include "aik_xml_parser.h"

// TODO: fix new lines in make_xml.py

int main(int argc, char *argv[]) {

  if ((argc - 1) < 1) {
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

  bool if_web_version = node_input.read_flag_value("if_web_version");
  if (not(if_web_version)) {
    std::string line(80, '-');
    std::cout << "A copy of the \"" << argv[1] << "\" input:\n";
    std::cout << line << "\n";
    node_input.print(std::cout);
    std::cout << line << "\n\n";
  }

  std::string job = node_input.read_string_value("job");

  std::ifstream xml_amu_file(ATOMIC_MASSES_FILE);
  if (!xml_amu_file.is_open()) {
    std::cout << "Atomic masses file " << ATOMIC_MASSES_FILE
              << " is not found in the current directory." << std::endl;
    exit(1);
  }
  xml_node node_amu_table("masses", xml_amu_file);

  bool done = false;
  if (job == "harmonic_pes")
    done = harmonic_pes_main(arg.c_str(), node_input, node_amu_table);

  if (!done)
    std::cout << "Method \"" << job << "\" is unknown, or it has failed. \n";

  std::cout << '\n'
            << "Job \"" << argv[0] << ' ' << argv[1]
            << "\" has finished: " << GetTime() << '\n';

  return EXIT_SUCCESS;
}
