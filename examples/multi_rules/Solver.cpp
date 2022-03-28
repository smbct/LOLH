#include "Solver.hpp"

#include <iostream>
#include <fstream>

/*----------------------------------------------------------------------------*/
Solver::Solver(LocalSearch::Instance& instance): _instance(instance), _sol(instance)
{

}

/*----------------------------------------------------------------------------*/
void Solver::solve() {

  /* random initialization */
  _sol.randomInit();

  /* local search */

  /* print the solution */
  for( int varIndex = 0; varIndex < _sol.size(); varIndex ++ ) {
    std::cout << _sol.getValue(varIndex) << " ";
  }
  std::cout << std::endl;

  std::cout << "cost: " << _sol.computeScore() << std::endl;

}

/*----------------------------------------------------------------------------*/
void Solver::saveSolution(std::string filename) {

  /* export the solution to a file */
  std::ofstream file(filename);
  file << _sol.size();
  for( int varIndex = 0; varIndex < _sol.size(); varIndex ++ ) {
    file << " " << _sol.getValue(varIndex);
  }
  file.close();

}
