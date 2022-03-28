#pragma once

#include "Instance.hpp"

/*!
 * \class Solver
 * \brief class managing the solutions to the problem
 */
class Solution {

  public:

    /*!
     * \brief constructor
     * \param instance the problem instance
     */
    Solution(LocalSearch::Instance& instance);

    /*!
    * \brief return the score of the solution
    * \return the score of the solution
    */
    double computeScore();

    /*!
     * \brief compute the score of one specific rule
     * \param ruleIndex the index of the rule
     * \return the score of the rule
     */
    double computeRuleScore(int ruleIndex);

    /*!
     * \brief random initialization of the values
     */
    void randomInit();

    /*!
     * \brief return the number of variables in the solution
     * \return the number of variables
     */
     int size();

    /*!
     * \brief get the value of a variable
     * \param variableIndex the index of the variable
     * \return the value of the variable
     */
    uint getValue(int variableIndex);

  private: /* private methods */



  private: /* private attributes */

    /* the problem instance */
    LocalSearch::Instance& _instance;

    /* the assignment of the solution */
    std::vector<uint> _var;

    /* the logical atoms associated to the instance */
    std::vector<std::pair<uint,uint>> _atoms;

    /* the positive and negative errors of the atoms for each rules */
    std::vector<std::vector<std::pair<uint,uint>>> _atomErrors;

};
