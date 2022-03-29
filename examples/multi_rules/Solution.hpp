#pragma once

#include "Instance.hpp"

#include <set>

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
    * \brief recompute the score of the solution
    */
    void recomputeScore();

    /*!
     * \brief get the score of the solution
     * \return the score of the solution
     */
    double score();

    /*!
     * \brief recompute the score of one specific rule
     * \param ruleIndex the index of the rule
     * \return the score of the rule
     */
    double recomputeRuleScore(int ruleIndex);

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

    /*!
     * \brief update the value of one variable, and recompute the score accordingly
     * \param varInd index of the variable to update
     * \param new value value of the variable
     */
    void updateVariable(int varInd, int value);


  private: /* private methods */



  private: /* private attributes */

    /* the problem instance */
    LocalSearch::Instance& _instance;

    /* the assignment of the solution */
    std::vector<uint> _var;

    /* the global score of the rules */
    double _score;

    /* the logical atoms associated to the instance */
    std::vector<std::pair<uint,uint>> _atoms;

    /* indexes of the atoms given the variable index and the logical value */
    std::vector<std::vector<int>> _atomIndexes;

    /* the positive and negative errors of the atoms for each rules */
    std::vector<std::vector<std::pair<uint,uint>>> _atomErrors;

    /* the score of the atoms for each rules */
    std::vector<std::vector<double>> _atomScores;

    /* sum of the scores of the atoms selected for each rule */
    std::vector<double> _sumAtomScores;

    /* set of indexes for the atoms selected for the body of the rules */
    std::vector<std::set<int>> _selectedAtoms;

    /* number of positive examples for each rule, deduced from the current assignment */
    std::vector<int> _nPositives;

    /* number of negative examples for each rule, deduced from the current assignment (if the positive example is not assigned to the rule, it is considered as negative) */
    std::vector<int> _nNegatives;

};
