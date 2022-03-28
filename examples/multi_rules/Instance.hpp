#pragma once

#include <vector>
#include "../../c++/src/DataFrame.hpp"

namespace LocalSearch {

struct Instance {

  public:
    int p_rules; /* number of rules to infer */
    double t; /* selection thresold for the rule atoms */
    DataFrame<uint> dataset; /* dataset of the problem */

    std::vector<uint> positives; /* positive examples indexes */
    std::vector<uint> negatives; /* negative examples indexes */

};


}
