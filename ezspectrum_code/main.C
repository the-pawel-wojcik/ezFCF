#include "genincludes.h"
#include "genutil.h"
#include "harmonic_pes_main.h"
#include "aik_xml_parser.h"


int main(int argc, char* argv[])
{

  if ((argc-1)<1) {
    std::cerr << argv[0] << ": one and only one argument required, <Job.xml>.\n";
    exit(1);
  }
  std::string arg=argv[1];

  std::cout << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been started: " << GetTime() << '\n';


  std::ifstream xml_file(arg.c_str()); 
  xml_node node_input("input",xml_file);
  bool if_web_version=node_input.read_flag_value("if_web_version");
  if (not(if_web_version)) {
    std::cout << "A copy of the \"" << argv[1] << "\" input:\n";
    std::cout << "------------------------------------------------------------------------------\n";
    node_input.print(std::cout);
    std::cout << "------------------------------------------------------------------------------\n \n";
  }

  
  std::string job=node_input.read_string_value("job");
   
  bool done = false;
  if (job == "harmonic_pes")
    done=harmonic_pes_main(arg.c_str());
  
  if ( !done )
    std::cout << "Method \"" << job <<"\" is unknown, or it has been failed. \n";
  
  std::cout << '\n' << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been finished: " <<  GetTime() << '\n';

  return EXIT_SUCCESS;
}


