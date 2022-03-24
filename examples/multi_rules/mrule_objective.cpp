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

  // Notice the minus here.
  // GHOST's solver tries to minimize any objective function.
  // Thus, for maximization problems like this one, outputing '- returned_value' does the trick.

  double score = 0;
  for(uint indRule = 0; indRule < _instance.p_rules; indRule ++) {
    score += compute_rule_score(variables, indRule);
  }

  std::cout << "score: " << score << std::endl;

  return score;
}
