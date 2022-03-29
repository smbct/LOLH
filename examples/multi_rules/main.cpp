/*!
 * \file main.cpp
 * \author S. Buchet
 * \brief experiments for multi-rule learning
 */

#include <iostream>

#include <fstream>

#include <cstdlib>

#include "Solver.hpp"
#include "Instance.hpp"


using namespace std;

/*-----------------------------------------------------------------------------*/
void readInstanceFile(string filename, LocalSearch::Instance& instance) {

  std::ifstream file(filename);

  uint n_positives, n_negatives;
  std::string label;
  file >> n_positives;
  std::cout << "n positives: " << n_positives << std::endl;
  for(uint i = 1; i <= n_positives; i ++) {
    file >> label;
    instance.positives.push_back(instance.dataset.getRowIndex(label));
  }

  file >> n_negatives;
  std::cout << "n negatives: " << n_negatives << std::endl;
  for(uint i = 1; i <= n_negatives; i ++) {
    file >> label;
    instance.negatives.push_back(instance.dataset.getRowIndex(label));
  }

  file.close();

}

/*-----------------------------------------------------------------------------*/
int main() {


  std::cout << "hello LOLH local search" << std::endl;

  /* random seed initialization */
  srand(42);

  LocalSearch::Instance instance;

  instance.p_rules = 3;
  // instance.t = 0.5;
  instance.t = 0.3;

  std::string matrix_filename = "../../dataset/Imagine/discrete_matrix.csv";
  instance.dataset.loadFromCSV(matrix_filename);
  instance.dataset.computeUniqueVal();
  std::cout << "done" << std::endl;

  /* read the instance (positive and negative examples) from a txt file */
  std::string instance_filename = "T_instance.txt";
  readInstanceFile(instance_filename, instance);

  Solver solver(instance);
  solver.solve();
  solver.saveSolution("T_temp_sol.txt");

  return 0;

}
