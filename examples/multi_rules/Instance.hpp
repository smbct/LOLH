#pragma once

#include <vector>
#include "../../c++/src/DataFrame.hpp"

namespace LocalSearch {

class Instance {

  public:

    /*!
     * \brief compute the score of an atom
     * \param atom the atom
     * \return the score of the atom
     */
    double computeAtomScore(std::pair<int,int> atom);

  public:
    
    int p_rules; /* number of rules to infer */
    double t; /* selection thresold for the rule atoms */
    DataFrame<uint> dataset; /* dataset of the problem */

    std::vector<uint> positives; /* positive examples indexes */
    std::vector<uint> negatives; /* negative examples indexes */

};


}
