#include "job_parameters.h"

JobParameters::JobParameters(const xml_node &node_input)
    : temperature(300), peak_thresh(0.001) {

  // TODO: make sure that the node exists and there is only one
  // This is a common test: make it a generally available function
  xml_node node_jobparams(node_input, "job_parameters", 0);
  // read global paramters
  temperature = node_jobparams.read_double_value("temperature");
  // fcf threshold (from the <job_parameters> tag)
  peak_thresh =
      node_jobparams.read_double_value("spectrum_intensity_threshold");
}
