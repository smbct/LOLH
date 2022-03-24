#include "model_builder.hpp"

#include "mrule_objective.hpp"

#include <fstream>

#include <string>

#include <iostream>

using namespace ghost;

/*-----------------------------------------------------------------------------*/
LPOBuilder::LPOBuilder() : ghost::ModelBuilder() {

  _instance.p_rules = 3;
  _instance.t = 0.5;

  std::string matrix_filename = "../../dataset/Imagine/discrete_matrix.csv";
  _instance.dataset.loadFromCSV(matrix_filename);
  _instance.dataset.computeUniqueVal();
  std::cout << "done" << std::endl;


  /* read the instance (positive and negative examples) from a txt file */
  std::string instance_filename = "T_instance.txt";

  std::ifstream file(instance_filename);

  uint n_pos, n_neg;
  std::string label;
  file >> n_pos;
  std::cout << "n pos: " << n_pos << std::endl;
  for(uint i = 1; i <= n_pos; i ++) {
    file >> label;
    _instance.positives.push_back(_instance.dataset.getRowIndex(label));
  }

  file >> n_neg;
  std::cout << "n neg: " << n_neg << std::endl;
  for(uint i = 1; i <= n_neg; i ++) {
    file >> label;
    _instance.negatives.push_back(_instance.dataset.getRowIndex(label));
  }

  file.close();

}

/*-----------------------------------------------------------------------------*/
void LPOBuilder::declare_variables() {

  /* one variable per positive example */
  /* value = id of the rule associated to the positive example */
  for(int i = 0; i < _instance.positives.size(); i ++) {
    variables.emplace_back(0, _instance.p_rules);
  }



}

/*-----------------------------------------------------------------------------*/
void LPOBuilder::declare_constraints() {

  /* no constraint */

}

/*-----------------------------------------------------------------------------*/
void LPOBuilder::declare_objective() {

  objective = std::make_shared<MRuleObjective>(_instance, variables);

}
