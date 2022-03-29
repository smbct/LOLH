#include "Solution.hpp"

#include <iostream>

/*----------------------------------------------------------------------------*/
Solution::Solution(LocalSearch::Instance& instance): _instance(instance)
{

  /* initialize the variables */
  _var.resize(_instance.positives.size());

  _atomIndexes.resize(_instance.dataset.nColumns());

  int atomIndex = 0;

  /* creation of the list of atoms */
  for(uint colInd = 0; colInd < _instance.dataset.nColumns(); colInd ++) {

    uint nValues =_instance.dataset.nUnique(colInd);

    _atomIndexes[colInd].resize(nValues);

    for(uint value = 0; value < nValues; value ++) {
      _atoms.push_back(std::pair<uint,uint>(colInd, value));
      _atomIndexes[colInd][value] = atomIndex;
      atomIndex += 1;
    }
  }

  /* creation of the list of atom errors for each rule */
  _atomErrors.resize(_instance.p_rules);
  for(int ruleInd = 0; ruleInd < _instance.p_rules; ruleInd ++) {
    _atomErrors[ruleInd].resize(_atoms.size());
  }

  /* initialization of the number of positive examples of each rule */
  _nPositives.resize(_instance.p_rules);

  /* initialization of the number of negative examples of each rule */
  _nNegatives.resize(_instance.p_rules);

  /* indexes of the atoms selected for the rule bodies */
  _selectedAtoms.resize(_instance.p_rules);

  _atomScores.resize(_instance.p_rules);
  for(int ruleInd = 0; ruleInd < _instance.p_rules; ruleInd ++) {
    _atomScores[ruleInd].resize(_atoms.size());
  }

}

/*----------------------------------------------------------------------------*/
int Solution::size() {
  return _var.size();
}

/*----------------------------------------------------------------------------*/
uint Solution::getValue(int variableIndex) {
  return _var[variableIndex];
}

/*----------------------------------------------------------------------------*/
void Solution::randomInit() {

  for(int varInd = 0; varInd < _var.size(); varInd ++) {
    _var[varInd] = rand()%_instance.p_rules;
  }

}

/*----------------------------------------------------------------------------*/
void Solution::updateVariable(int varInd, int value) {

  if(_var[varInd] == value) {
    std::cout << "nothing to change" << std::endl;
    return;
  }

  /* two situations */
  /* 1) the example was a positive example and is not anymore */
  /* 2) the example was not a positive example and it is now one */

  /* update the score only for the atoms which do not verify the example values */

  /* 1) -> varInd is no longer a positive example for the rule of index "_var[varInd]" */
  int ruleIndex = _var[varInd];

  // _selectedAtoms[ruleIndex].clear();
  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];

    int dataVal = _instance.dataset.getData(_instance.positives[varInd], atom.first);


    double prevScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());
    prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    double newScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());

    if(atom.second != dataVal) {
      _atomErrors[ruleIndex][atomInd].first -= 1;
    }

    newScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]-1);

    /* check if the new score make the atom selected, and update the solution accordingly */
    if(prevScore < _instance.t && newScore >= _instance.t) {
      _selectedAtoms[ruleIndex].insert(atomInd);
    } else if(prevScore >= _instance.t && newScore < _instance.t) {
      _selectedAtoms[ruleIndex].erase(atomInd);
    }
    // if(newScore >= _instance.t) {
    //   _selectedAtoms[ruleIndex].insert(atomInd);
    // }

  }

  /* update the number of positive examples for this rule */
  _nPositives[ruleIndex] -= 1;




  /* 2) -> varInd is now positive for the rule of index "value" */
  ruleIndex = value;

  // _selectedAtoms[ruleIndex].clear();

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];
    int dataVal = _instance.dataset.getData(_instance.positives[varInd], atom.first);


    double prevScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());
    prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    double newScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());

    /* update the atom positive error */
    if(atom.second != dataVal) {
      _atomErrors[ruleIndex][atomInd].first += 1;
    }

    newScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]+1);

    /* check if the new score make the atom unselected, and update the solution accordingly */
    if(prevScore >= _instance.t && newScore < _instance.t) {
      _selectedAtoms[ruleIndex].erase(atomInd);
    } else if(prevScore < _instance.t && newScore >= _instance.t) {
      _selectedAtoms[ruleIndex].insert(atomInd);
    }

    // if(newScore >= _instance.t) {
    //   _selectedAtoms[ruleIndex].insert(atomInd);
    // }

  }


  /* update the number of positive examples for this rule */
  _nPositives[ruleIndex] += 1;

  ruleIndex = -1;


  


  /* perform the update */
  _var[varInd] = value;



  /* recompute the score */
  /* this is always required since the number of positive examples change for all atoms in the concerned rules */
  _score = 0;
  for(int ruleInd = 0; ruleInd < _instance.p_rules; ruleInd ++) {

    double ruleScore = 0.;

    for(int atomIndex : _selectedAtoms[ruleInd]) {
      double atomScore = static_cast<double>(_atomErrors[ruleInd][atomIndex].second)/static_cast<double>(_instance.negatives.size());
      atomScore -= static_cast<double>(_atomErrors[ruleInd][atomIndex].first)/static_cast<double>(_nPositives[ruleInd]);

      _atomScores[ruleInd][atomIndex] = atomScore;

      ruleScore += atomScore;
    }

    ruleScore /= static_cast<double>(_selectedAtoms[ruleInd].size());
    _score += ruleScore;
  }


}

/*----------------------------------------------------------------------------*/
double Solution::score() {
  return _score;
}

/*----------------------------------------------------------------------------*/
void Solution::recomputeScore() {

  double prevScore = _score;

  _score = 0;
  for(uint ruleIndex = 0; ruleIndex < _instance.p_rules; ruleIndex ++) {
    _score += recomputeRuleScore(ruleIndex);
  }

  std::cout << "recomputed score: " << _score << std::endl;

  if(prevScore != _score) {
    std::cout << std::endl << "!!!!!! error scores are not equal !!!!!! " << std::endl << std::endl;
  }

}

/*----------------------------------------------------------------------------*/
double Solution::recomputeRuleScore(int ruleIndex) {

  std::cout << "recomputation of rule " << ruleIndex << " score !" << std::endl;

  int nPosPrev = _nPositives[ruleIndex];

  _nPositives[ruleIndex] = 0;

  /* compute the positive examples selected for this rule */
  std::vector<uint> positives;
  for(int varInd = 0; varInd < _var.size(); varInd ++) {
    if(_var[varInd] == ruleIndex) {
      positives.push_back(_instance.positives[varInd]);
      _nPositives[ruleIndex] ++; /* one more positive exampe */
    }
  }

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];

    _atomErrors[ruleIndex][atomInd] = std::pair<uint,uint>(0,0);

    for(uint pos : positives) {
      uint val = _instance.dataset.getData(pos, atom.first);
      if(val != atom.second) {
        _atomErrors[ruleIndex][atomInd].first ++;
      }
    }

    for(uint neg: _instance.negatives) {
      uint val = _instance.dataset.getData(neg, atom.first);
      if(val != atom.second) {
        _atomErrors[ruleIndex][atomInd].second ++;
      }
    }

  }

  /* list of selected atoms for the rule */
  _selectedAtoms[ruleIndex].clear();

  double sumScore = 0.;

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    double atomScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());
    atomScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    _atomScores[ruleIndex][atomInd] = atomScore;

    if(atomScore >= _instance.t) {
      _selectedAtoms[ruleIndex].insert(atomInd);
      sumScore += atomScore;
    }

  }

  std::cout << _selectedAtoms[ruleIndex].size() << " atoms selected for rule " << ruleIndex << std::endl;

  double ruleScore = sumScore / static_cast<double>(_selectedAtoms[ruleIndex].size());

  std::cout << " debug rule " << ruleIndex << " score: " << ruleScore << std::endl;

  std::cout << "selected atoms: ";
  for(int atomInd : _selectedAtoms[ruleIndex]) {
    std::cout << atomInd << " ";
  }
  std::cout << std::endl;


  return ruleScore;

}
