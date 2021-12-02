/*!
 * \file NetworkInduction.hpp
 * \author S. Buchet
 * \brief definition of class NetworkInduction
 */

#pragma once

#ifndef NETWORK_INDUCTION_HPP
#define NETWORK_INDUCTION_HPP

#include <string>

#include "DataFrame.hpp"

/*!
* \class NetworkInduction
 * \brief functions to create a network of variables given their correlations (static class)
 */
class NetworkInduction {

  public:

    /*!
     * \brief compute a co-expression network for all genes (fast version)
     * \param dataset the dataset
     * \param selectionThreshold the min relative distance of selected genes
     * \param output_filename the name of the ouput file
     * \param negativeCells ???
     * \param maxEdges max number of edges, if 0, no limits
     */
    static void computeNetwork(DataFrame<uint>& dataset, double selectionThreshold, std::string output_filename, std::vector<bool>* negativeCells = nullptr, uint maxEdges = 0);

    /*!
     * \brief compute a regulatory network for all genes (fast version)
     * \param dataset the dataset
     * \param successors the neighbourhood graph
     * \param trRate the transition rate to select a positive state
     * \param predNeq indicates if only transitions with a different value for the target are selected
     * \param selectionThreshold the min relative distance of selected genes
     * \param output_filename the name of the ouput file
     * \param max_edges maximum number of edges for an atom. 0 is the defualt value, meaning no limitations
     */
    static void computeRegulationNetwork(DataFrame<uint>& dataset, NGraph successors, double trRate, uint predNeq, double selectionThreshold, std::string output_filename, uint max_edges = 0);

    /*!
     * \brief compute a co-expression network for all genes
     * \param dataset the dataset
     * \param bodyLength length of the supported solutions
     * \param nGenes number of genes to select
     * \param rate rate of the best rules to select from the preto front
     * \param output_filename the name of the ouput file
     */
    static void computeNetworkOld(DataFrame<uint>& dataset, uint bodyLength, uint nGenes, double rate, std::string output_filename);

};

#endif // NETWORK_INDUCTION_HPP
