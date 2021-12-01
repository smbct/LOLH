/*!
 * \file Parameters.hpp
 * \author S. Buchet
 * \brief definition of Parameters data structure
 */

#pragma once

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>

/*!
 * \class Parameters
 * \brief data structure to store the program parameters
 */
class Parameters {

  public: /* public methods */

    /*!
     * \brief default constructor */
    Parameters() {

      debug = false;

      verbose = false;

      compute_network = true;

      coexpression = true;

      input_matrix = "";

      output_file = "";

      threshold = -1.;

      max_edges = 0;

      /****************************/

      input_transitions = "";

      transition_rate = -1.;

      predecessor_neq = 0;

      transition_delay = 1;
    }


  public: /* public attributes */

    /* true: bypass the parameters and execute the debug function */
    bool debug;

    /* true: display addition information on the computation */
    bool verbose;

    /* true if iference of the coexpresion graph, false if dynamic graph */
    bool coexpression;

    /* true to compute a network, false to only assess the quality */
    bool compute_network;

    /* input matrix, single-cell gene expression csv file */
    std::string input_matrix;

    /* output file of the network */
    std::string output_file;

    /* atom score threshold for the rule inference */
    double threshold;


    /* parameters for the dynamic inference */

    /* input transition file, for the dynamic network inference */
    std::string input_transitions;

    /* percentage of transitions (between 0 and 1) that need to lead to the learned atom in order to be considered as a positive example */
    double transition_rate;

    /* if 1, transitions are selected only of the predecessor atom is not equal to the successor atom (force a difference), if 2: selection occurs only when the predecessor is equal, otherwise if 0, no constraint  */
    uint predecessor_neq;

    /* consider as transitions the extremity of a path of a certain length */
    uint transition_delay;

    /* maximum number of edges in the coexpression or dynamical graph */
    uint max_edges;

};

#endif
