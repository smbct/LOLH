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

  _sol.assignFromFile("B_assignment.txt");
  _sol.recomputeScore();
  std::cout << "cost: " << _sol.score() << std::endl;

  /* local search */

  /* n random moves */
  // for(int ind = 0; ind < 100; ind ++) {
  //
  //   int varIndex = rand()%_sol.size();
  //
  //   std::cout << "move " << ind << std::endl;
  //   _sol.updateVariable(varIndex, rand()%_instance.p_rules);
  //   std::cout << "move done" << std::endl << std::endl;
  //
  //   std::cout << "cost: " << _sol.score() << std::endl;
  //
  // }

  bool stop = true;
  bool improved = false;

  int nbIt = 0;
  int nbRestart = 0;
  int nbFullRestart = 0;

  Solution bestSol(_instance);

  while(!stop) {

    // improved = firstImprovement();

    // if(!improved) {
      improved = bestImprovement();
    // }

    if(!improved) {

      if(_sol.score() > bestSol.score()) {
        bestSol.copyFrom(_sol);
      }

      if(nbRestart < 4) {
        _sol.disruption(50);
        nbRestart ++;
      } else {
        if(nbFullRestart < 10) {
          nbRestart = 0;
          nbFullRestart ++;
          _sol.randomInit();
          _sol.recomputeScore();
        } else {
          stop = true;
        }
      }

    }

    // if(nbIt % 20 == 0) {
    //
    // }

    std::cout << "score: " << _sol.score() << std::endl;
    std::cout << " best so far: " << bestSol.score() << std::endl;
    nbIt += 1;
  }

  std::cout << "nb iterations: " << nbIt << std::endl;

  std::cout << std::endl << "recompute solution" << std::endl << std::endl;
  _sol.recomputeScore();

  /* print the solution */
  for( int varIndex = 0; varIndex < _sol.size(); varIndex ++ ) {
    std::cout << _sol.getValue(varIndex) << " ";
  }
  std::cout << std::endl;

  std::cout << "final score: " << _sol.score() << std::endl;

  if(_sol.score() < bestSol.score()) {
    _sol.copyFrom(bestSol);
  }

  saveSolution("B_sol.txt");
  std::cout << "solution saved!" << std::endl;

}

/*----------------------------------------------------------------------------*/
bool Solver::bestImprovement() {

  bool improved = false;

  double initialScore = _sol.score();

  int bestIndVar;
  int bestValue;
  double bestScore = -1;

  for(int indVar = 0; indVar < _sol.size(); indVar ++) {

    int prevVal = _sol.getValue(indVar);

    /* test all the changes */
    for(int value = 0; value < _instance.p_rules; value ++) {
      if(value != prevVal) {

        /* apply the move */
        _sol.updateVariable(indVar, value);

        if(bestScore < -1 || _sol.score() > bestScore) {
          bestIndVar = indVar;
          bestValue = value;
          bestScore = _sol.score();
        }

        /* undo the move */
        _sol.updateVariable(indVar, prevVal);
      }
    }

  }

  improved = (bestScore > initialScore);

  if(improved) {
    _sol.updateVariable(bestIndVar, bestValue);
  }

  return improved;

}

/*----------------------------------------------------------------------------*/
bool Solver::firstImprovement() {

  bool improved = false;

  double initialScore = _sol.score();

  bool stop = false;

  int indVar = 0;
  while(!stop && indVar < _sol.size()) {

    int prevVal = _sol.getValue(indVar);

    int value = 0;
    while(!stop && value < _instance.p_rules) {

      if(value != prevVal) {

        /* apply the move */
        _sol.updateVariable(indVar, value);

        if(_sol.score() > initialScore) {
          stop = true;
          improved = true;
        } else {
          /* undo the move */
          _sol.updateVariable(indVar, prevVal);
        }

      }

      value ++;
    }

    indVar ++;
  }

  return improved;

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
