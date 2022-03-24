#pragma once

#include <vector>
#include "../../c++/src/DataFrame.hpp"

struct instance_data {

  int p_rules; /* number of rules to infer */
  double t; /* selection thresold for the rule atoms */
  DataFrame<uint> dataset; /* dataset of the problem */

  std::vector<uint> positives; /* positive examples idnexes */
  std::vector<uint> negatives; /* negative examples indexes */

};
