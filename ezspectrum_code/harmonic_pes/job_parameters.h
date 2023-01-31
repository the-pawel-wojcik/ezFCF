#ifndef _job_parameters_h_
#define _job_parameters_h_

#include "aik_xml_parser.h"

/* Class for the "job_parameters" input node. */
class JobParameters {
  private:
    double temperature;
    double peak_thresh;
    JobParameters() {}
  public:
    JobParameters(const xml_node &node_input);
    double get_temp() const { return temperature; }
    double get_intensity_thresh() const { return peak_thresh; }
};

#endif
