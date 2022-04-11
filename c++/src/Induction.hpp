/*!
 * \file Induction.hpp
 * \author S. Buchet
 * \brief definition of class Induction
 */

#pragma once

#ifndef INDUCTION_HPP
#define INDUCTION_HPP

#include "DataFrame.hpp"
#include "Instance.hpp"

/*
 * \class Parameters
 * \brief Parameters of the learning alorithm
 */
class Parameters {

  public:

    Parameters() { }

    bool coExpression; /* learn a co-xpression or a regulation program */
    uint bodyLength; /* length of the body */
    uint clusterBodyLength; /* body length of a cluster instance */
    uint graphDelay; /* distance between cells of a transition in the graph */
    uint nRules; /* number of rules to learn */
    double posRate; /* percentage of positive successor to consider the transition as positive */
    uint matchingThreshold; /* min error to accept a sample on a rule */
    double falsePositiveRate; /* max false positive rate authorized -> other rules are discared */

};

/*
 * \class Induction
 * \brief Manager of LFIT program Induction from transitions dataset
 */
class Induction {

  public:

    /*!
     * \brief constructor
     */
    Induction();

    /*!
     * \brief induction of multiple rules from a classification instance
     * \param instance the classification instance
     * \param param parameters of the algorithm
     * \param bodies list of bodies computed
     * \param clusterInstance true if the instance corresponds to a cluster label
     */
    static void ruleInduction(const Instance& instance, Parameters& param, std::vector<Body>& bodies, bool clusterInstance);

    /*!
     * \brief create a logic program containing the rules for variables in [indMin,indMax]
     * \param datasetFileName the name of the dataset file
     * \param transitionsFileName the name the the transition file
     * \param labelsFileName name of the file containing the cell labels
     * \param parameters of the learning algorithm
     */
    static void createLFITProgram(std::string datasetFileName, std::string transitionsFileName, std::string labelsFileName, Parameters& param);

  private:


};

#endif /* INDUCTION_HPP */
