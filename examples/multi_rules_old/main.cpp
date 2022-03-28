/*!
 * \file main.cpp
 * \author S. Buchet
 * \brief experiments for multi-rule learning
 */

#include <iostream>

#include <fstream>

#include <ghost/solver.hpp>

#include "model_builder.hpp"

using namespace std;

/*-----------------------------------------------------------------------------*/
int main() {


  cout << "hello GHOST" << endl;


  // Declaring the model builder
  LPOBuilder builder;

  // Defining the solver and calling it
  ghost::Solver<LPOBuilder> solver( builder );

  double cost;
  std::vector<int> solution;

  // Run the solver with a 500 microseconds budget
  // bool found = solver.solve( cost, solution, 500us );
  // bool found = solver.solve( cost, solution, 1800s );
  bool found = solver.solve( cost, solution, 30s );

  // After 500 microseconds, the solver will write in cost and solution the best solution it has found.
  found ? std::cout << "Solution found\n" : std::cout << "Solution not found\n";
  std::cout << "Cost: " << cost << "\nSolution:\n";

  // for( auto v : solution ) {
  //   std::cout << " " << v;
  // }

  std::cout << "\n";

  /* export the solution to a file */
  std::ofstream file("temp_sol.txt");
  file << solution.size();
  for( auto v : solution ) {
    file << " " << v;
  }
  file.close();

  return 0;

}
