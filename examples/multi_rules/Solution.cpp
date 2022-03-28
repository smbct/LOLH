#include "Solution.hpp"

#include <iostream>

/*----------------------------------------------------------------------------*/
Solution::Solution(LocalSearch::Instance& instance): _instance(instance)
{

  /* initialize the variables */
  _var.resize(_instance.positives.size());

  /* creation of the list of atoms */
  for(uint indVar = 0; indVar < _instance.dataset.nColumns(); indVar ++) {
    for(uint indValue = 0; indValue < _instance.dataset.nUnique(indVar); indValue ++) {
      _atoms.push_back(std::pair<uint,uint>(indVar, indValue));
    }
  }

  /* creation of the list of atom errors for each rule */
  _atomErrors.resize(_instance.p_rules);
  for(int indRule = 0; indRule < _instance.p_rules; indRule ++) {
    _atomErrors[indRule].resize(_atoms.size());
  }

}

/*----------------------------------------------------------------------------*/
int Solution::size() {
  return _var.size();
}


uint Solution::getValue(int variableIndex) {
  return _var[variableIndex];
}

/*----------------------------------------------------------------------------*/
void Solution::randomInit() {

  for(int indVar = 0; indVar < _var.size(); indVar ++) {
    _var[indVar] = rand()%_instance.p_rules;
  }

}


/*----------------------------------------------------------------------------*/
double Solution::computeScore() {

  double score = 0;
  for(uint ruleIndex = 0; ruleIndex < _instance.p_rules; ruleIndex ++) {
    score += computeRuleScore(ruleIndex);
  }

  std::cout << "score: " << score << std::endl;

  return score;

}

/*----------------------------------------------------------------------------*/
double Solution::computeRuleScore(int ruleIndex) {

  /* compute the positive examples selected for this rule */
  std::vector<uint> positives;
  for(int indVar = 0; indVar < _var.size(); indVar ++) {
    if(_var[indVar] == ruleIndex) {
      positives.push_back(_instance.positives[indVar]);
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

  std::cout << selectedAtoms.size() << " atoms selected for rule " << ruleIndex << std::endl;

  rule_score /= static_cast<double>(selectedAtoms.size());

  return rule_score;

}
