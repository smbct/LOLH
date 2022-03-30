#include "Solution.hpp"

#include <iostream>
#include <fstream>

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

      double atomScore = _instance.computeAtomScore(std::pair<int,int>(colInd,value));
      if(atomScore >= 0) {
        _atoms.push_back(std::pair<uint,uint>(colInd, value));
        _atomIndexes[colInd][value] = atomIndex;
        atomIndex += 1;
      }

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

  _score = -1.;

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

  /* 1) -> varInd is no longer a positive example for the rule of index "_var[varInd]", it is now a negative example */
  int ruleIndex = _var[varInd];

  // _selectedAtoms[ruleIndex].clear();
  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];

    int dataVal = _instance.dataset.getData(_instance.positives[varInd], atom.first);

    double prevScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_nNegatives[ruleIndex]);
    prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    if(atom.second != dataVal) {
      _atomErrors[ruleIndex][atomInd].first -= 1;
      _atomErrors[ruleIndex][atomInd].second += 1;
    }

    double newScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_nNegatives[ruleIndex]+1);
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

  /* update the number of positive and negative examples for this rule */
  _nPositives[ruleIndex] -= 1;
  _nNegatives[ruleIndex] += 1;



  /* 2) -> varInd is now positive for the rule of index "value" */
  ruleIndex = value;

  // _selectedAtoms[ruleIndex].clear();

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];
    int dataVal = _instance.dataset.getData(_instance.positives[varInd], atom.first);

    double prevScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_nNegatives[ruleIndex]);
    prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    /* update the atom positive error */
    if(atom.second != dataVal) {
      _atomErrors[ruleIndex][atomInd].first += 1;
      _atomErrors[ruleIndex][atomInd].second -= 1;
    }

    double newScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_nNegatives[ruleIndex]-1);
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


  /* update the number of positive and negative examples for this rule */
  _nPositives[ruleIndex] += 1;
  _nNegatives[ruleIndex] -= 1;

  ruleIndex = -1;





  /* perform the update */
  _var[varInd] = value;



  /* recompute the score */
  /* this is always required since the number of positive examples change for all atoms in the concerned rules */
  _score = 0;
  for(int ruleInd = 0; ruleInd < _instance.p_rules; ruleInd ++) {

    double ruleScore = 0.;

    for(int atomIndex : _selectedAtoms[ruleInd]) {
      double atomScore = static_cast<double>(_atomErrors[ruleInd][atomIndex].second)/static_cast<double>(_nNegatives[ruleInd]);
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

}

/*----------------------------------------------------------------------------*/
double Solution::recomputeRuleScore(int ruleIndex) {

  _nPositives[ruleIndex] = 0;
  _nNegatives[ruleIndex] = _instance.negatives.size();

  /* compute the positive and negative examples selected for this rule */
  std::vector<uint> positives, negatives;
  for(int varInd = 0; varInd < _var.size(); varInd ++) {
    if(_var[varInd] == ruleIndex) {
      positives.push_back(_instance.positives[varInd]);
      _nPositives[ruleIndex] ++; /* one more positive exampe */
    } else {
      negatives.push_back(_instance.positives[varInd]);
      _nNegatives[ruleIndex] ++; /* one more negative exampe */
    }
  }

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    auto atom = _atoms[atomInd];

    _atomErrors[ruleIndex][atomInd] = std::pair<uint,uint>(0,0);

    for(uint positive : positives) {
      uint val = _instance.dataset.getData(positive, atom.first);
      if(val != atom.second) {
        _atomErrors[ruleIndex][atomInd].first ++;
      }
    }

    /* take into account the positive examples not selected for this rule as negative examples */
    for(uint negative : negatives) {
      uint val = _instance.dataset.getData(negative, atom.first);
      if(val != atom.second) {
        _atomErrors[ruleIndex][atomInd].second ++;
      }
    }

    for(uint negative : _instance.negatives) {
      uint val = _instance.dataset.getData(negative, atom.first);
      if(val != atom.second) {
        _atomErrors[ruleIndex][atomInd].second ++;
      }
    }



  }

  /* list of selected atoms for the rule */
  _selectedAtoms[ruleIndex].clear();

  double sumScore = 0.;

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    double atomScore = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_nNegatives[ruleIndex]);
    atomScore -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    _atomScores[ruleIndex][atomInd] = atomScore;

    if(atomScore >= _instance.t) {
      _selectedAtoms[ruleIndex].insert(atomInd);
      sumScore += atomScore;
    }

  }

  double ruleScore = sumScore / static_cast<double>(_selectedAtoms[ruleIndex].size());

  return ruleScore;

}

/*----------------------------------------------------------------------------*/
void Solution::disruption(int nb) {

  for(int i = 1; i <= nb; i ++) {
    int varInd = rand()%_var.size();
    int value = rand()%_instance.p_rules;
    if(_var[varInd] == value) {
      value ++;
      value = value % _instance.p_rules;
    }
    updateVariable(varInd, value);
  }

}

/*----------------------------------------------------------------------------*/
void Solution::copyFrom(Solution& other) {
  _var = other._var;
  _score = other._score;
}

/*----------------------------------------------------------------------------*/
void Solution::assignFromFile(std::string filename) {

  std::ifstream file(filename);

  int nb = 0;
  file >> nb;
  for(int varInd = 0; varInd < _var.size(); varInd ++) {
    file >> nb;
    _var[varInd] = nb;
  }

  file.close();
}
