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
  _sol.recomputeScore();

  /* local search */

  /* n random moves */
  for(int ind = 0; ind < 100; ind ++) {

    int varIndex = rand()%_sol.size();

    std::cout << "move " << ind << std::endl;
    _sol.updateVariable(varIndex, rand()%_instance.p_rules);
    std::cout << "move done" << std::endl << std::endl;


    std::cout << "cost: " << _sol.score() << std::endl;

    if(ind % 10 == 0) {
      _sol.recomputeScore();
      std::cout << std::endl << "cost re-computed: " << _sol.score() << std::endl << std::endl;
    }

  }

  std::cout << std::endl << "recompute solution" << std::endl << std::endl;
  _sol.recomputeScore();

  /* print the solution */
  // for( int varIndex = 0; varIndex < _sol.size(); varIndex ++ ) {
  //   std::cout << _sol.getValue(varIndex) << " ";
  // }
  // std::cout << std::endl;

  std::cout << "cost: " << _sol.score() << std::endl;

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
