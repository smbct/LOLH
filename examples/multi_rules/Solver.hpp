#pragma once

#include "Instance.hpp"
#include "Solution.hpp"

/*!
 * \class Solver
 * \brief main class for the local search alorithm
 */
class Solver {

  public:

    /*!
     * \brief constructor
     * \param instance the problem instance
     */
    Solver(LocalSearch::Instance& instance);

    /*!
     * \brief solve the instance
     */
    void solve();

    /*!
     * \brief perform best improvement search on the solution
     * \return true if the solution score has been improved
     */
    bool bestImprovement();

    /*!
     * \brief perform first improvement search on the solution
     * \return true if the solution score has been improved
     */
    bool firstImprovement();

    /*!
     * \brief save the solution into a file
     * \param filename name of the file
     */
    void saveSolution(std::string filename);

  private:

    /* the problem instance */
    LocalSearch::Instance& _instance;

    /* the solution of the problem */
    Solution _sol;

};
