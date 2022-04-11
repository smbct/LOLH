#include "mrule_objective.hpp"


/*-----------------------------------------------------------------------------*/
MRuleObjective::MRuleObjective(instance_data& instance, const std::vector< ghost::Variable >& var) : ghost::Maximize(var,  "Multi rule scores" ), _instance(instance)
{

  /* creation of the list of atoms */
  for(uint indVar = 0; indVar < _instance.dataset.nColumns(); indVar ++) {
    for(uint indValue = 0; indValue < _instance.dataset.nUnique(indVar); indValue ++) {
      _atoms.push_back(std::pair<uint,uint>(indVar, indValue));
    }
  }

  _prev_values.resize(var.size(), 0);

  // std::cout << "initial values" << std::endl;
  // for(uint ind = 0; ind < var.size(); ind ++) {
  //   std::cout << "var " << ind << ": "  << var[ind].get_value() << std::endl;
  // }

}

/*-----------------------------------------------------------------------------*/
double MRuleObjective::compute_rule_score(const std::vector< ghost::Variable*>& variables, uint rule_id) const {

  /* compute the positive examples selected for this rule */
  std::vector<uint> positives;
  for(ghost::Variable* var : variables) {
    if(var->get_value() == rule_id) {
      positives.push_back(_instance.positives[var->get_id()]);
    }
  }



  /* computation of the scores for each atoms */
  std::vector<uint> pos_scores(_atoms.size(), 0);
  std::vector<uint> neg_scores(_atoms.size(), 0);

  for(uint indAtom = 0; indAtom < _atoms.size(); indAtom ++) {
    auto atom = _atoms[indAtom];

    for(uint pos : positives) {
      uint val = _instance.dataset.getData(pos, atom.first);
      if(val != atom.second) {
        pos_scores[indAtom] ++;
      }
    }

    for(uint neg: _instance.negatives) {
      uint val = _instance.dataset.getData(neg, atom.first);
      if(val != atom.second) {
        neg_scores[indAtom] ++;
      }
    }

  }

  /* list of selected atoms for the rule */
  std::vector<uint> selectedAtoms;

  double rule_score;

  for(uint indAtom = 0; indAtom < _atoms.size(); indAtom ++) {

    double score = static_cast<double>(neg_scores[indAtom])/static_cast<double>(_instance.negatives.size());
    score -= static_cast<double>(pos_scores[indAtom])/static_cast<double>(positives.size());

    if(score >= _instance.t) {
      rule_score += score;
      selectedAtoms.push_back(indAtom);
    }

  }

  std::cout << selectedAtoms.size() << " atoms selected for rule " << rule_id << std::endl;

  rule_score /= static_cast<double>(selectedAtoms.size());

  return rule_score;
}

/*-----------------------------------------------------------------------------*/
double MRuleObjective::required_cost( const std::vector< ghost::Variable*>& variables ) const {

  /* compare the current and the previous variables assignment*/
  // std::cout << "variable assignment comparison" << std::endl;
  // for(uint ind = 0; ind < variables.size(); ind ++) {
  //   // if(variables[ind]->get_value() != _prev_values[ind]) {
  //     std::cout << "var " << ind << ": "  << _prev_values[ind] << " -> " << variables[ind]->get_value() << std::endl;
  //   // }
  // }

  double score = 0;
  for(uint indRule = 0; indRule < _instance.p_rules; indRule ++) {
    score += compute_rule_score(variables, indRule);
  }

  std::cout << "score: " << score << std::endl;



  return score;
}

/*-----------------------------------------------------------------------------*/
void MRuleObjective::conditional_update_data_structures(const std::vector<ghost::Variable*> &variables, int index, int new_value) {

  std::cout << "update v" << index << " -> " << new_value << std::endl;

  /* store the assignment*/
  // for(uint ind = 0; ind < variables.size(); ind ++) {
  //   _prev_values[ind] = variables[ind]->get_value();
  // }
}
