/*!
 * \file ClassificationQuality.hpp
 * \author S. Buchet
 * \brief definition of class ClassificationQuality
 */

#pragma once

#ifndef CLASSIFICATION_QUALITY_HPP
#define CLASSIFICATION_QUALITY_HPP

#include "DataFrame.hpp"

/*!
 * \class ClassificationQuality
 * \brief functions to measure the quality of a classification instance, based on the paret front
 */
class ClassificationQuality {

  public: /* public methods */

    /*!
     * \brief compute the pareto front quality of the co-expression instance
     * \param dataset dataset containing gene expressions
     * \param bodyLength length of the supported solutions
     * \param output_filename name of the output file
     */
    static void computeCoexprQuality(DataFrame<uint>& dataset, uint bodyLength, std::string output_filename);

    /*!
     * \brief compute the pareto front quality of the co-expression instance: compute only the relative distance of the farthest point
     * \param dataset dataset containing gene expressions
     * \param bodyLength length of the supported solutions
     * \param output_filename name of the output file
     */
    static void computeFastCoexprQuality(DataFrame<uint>& dataset, uint bodyLength, std::string output_filename);

    /*!
     * \brief compute the quality of the co-expression instance: quality as relative area computed directly from the atoms
     * \param dataset dataset containing gene expressions
     * \param output_filename name of the output file
     */
    static void computeAtomCoexprQuality(DataFrame<uint>& dataset, std::string output_filename);

    /*!
     * \brief compute the pareto front quality of the co-expression instance
     * \param dataset dataset containing gene expressions
     * \param successors the glist of successors from the neighborhood graph for each cell
     * \param bodyLength length of the supported solutions
     * \param parameter: trRate percentage of successors verifying the right value to consider as a transition
     * \param parameter: predNeq 1 if transitions are extracted only when the predecessor state is not equal to the target value, 2 if transitions extracted whene the value is equal, 0 otherwise
     * \param output_filename name of the output file
     */
    static void computeRegulQuality(DataFrame<uint>& dataset, NGraph successors, uint bodyLength, double trRate, uint predNeq, std::string output_filename);



  private: /* private methods */

    /*!
     * \brief default constructor
     */
    ClassificationQuality();

};

#endif
