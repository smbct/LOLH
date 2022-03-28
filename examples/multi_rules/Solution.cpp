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

  /* initialize the vector of the sum of atom scores */
  _sumAtomScores.resize(_instance.p_rules);

  /* initialization of the umber of atoms selected in the body of each rule */
  _nAtomsBody.resize(_instance.p_rules);

  /* initialization of the number of positive examples of each rule */
  _nPositives.resize(_instance.p_rules);

  /* indexes of the atoms selected for the rule bodies */
  _selectedAtoms.resize(_instance.p_rules);
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

  /* two situations */
  /* 1) the example was a positive example and is not anymore */
  /* 2) the example was not a positive example and it is now one */

  /* update the score only for the atoms which do not verify the example values */

  bool change = false;

  /* 1) -> varInd is no longer a positive example for the rule of index "_var[varInd]" */
  int ruleIndex = _var[varInd];

  for(int colInd = 0; colInd < _atomIndexes.size(); colInd ++) {
    int dataVal = _instance.dataset.getData(varInd, colInd);

    /* iterate over the atoms not verifying dataVal */
    uint nValues =_instance.dataset.nUnique(colInd);


    for(uint atomValue = 0; atomValue < nValues; atomValue ++) {

      int atomIndex = _atomIndexes[colInd][atomValue];

      double prevScore, newScore;

      prevScore = static_cast<double>(_atomErrors[ruleIndex][atomIndex].second)/static_cast<double>(_instance.negatives.size());
      prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first)/static_cast<double>(_nPositives[ruleIndex]);

      newScore = static_cast<double>(_atomErrors[ruleIndex][atomIndex].second)/static_cast<double>(_instance.negatives.size());

      if(atomValue != dataVal) {
        newScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first-1)/static_cast<double>(_nPositives[ruleIndex]-1);
      } else {
        newScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first)/static_cast<double>(_nPositives[ruleIndex]-1);
      }

      /* check if the new score make the atom selected, and update the solution accordingly */
      if(prevScore < _instance.t && newScore > _instance.t) {
        change = true;
        _nAtomsBody[ruleIndex] += 1;
        _sumAtomScores[ruleIndex] += newScore-prevScore;
        _selectedAtoms[ruleIndex].insert(_atomIndexes[colInd][atomValue]);
      }

      /* update the atom positive error */
      if(atomValue != dataVal) {
        _atomErrors[ruleIndex][atomIndex].first -= 1;
      }

    }

  }





  /* 2) -> varInd is now positive for the rule of index "value" */
  ruleIndex = value;

  for(int colInd = 0; colInd < _atomIndexes.size(); colInd ++) {

    int dataVal = _instance.dataset.getData(varInd, colInd);

    /* iterate over the atoms not verifying dataVal */
    uint nValues =_instance.dataset.nUnique(colInd);

    for(uint atomValue = 0; atomValue < nValues; atomValue ++) {

      int atomIndex = _atomIndexes[colInd][atomValue];

      double prevScore, newScore;

      prevScore = static_cast<double>(_atomErrors[ruleIndex][atomIndex].second)/static_cast<double>(_instance.negatives.size());
      prevScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first)/static_cast<double>(_nPositives[ruleIndex]);

      newScore = static_cast<double>(_atomErrors[ruleIndex][atomIndex].second)/static_cast<double>(_instance.negatives.size());
      if(atomValue == dataVal) {
        newScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first+1)/static_cast<double>(_nPositives[ruleIndex]+1);
      } else {
        newScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first)/static_cast<double>(_nPositives[ruleIndex]+1);
      }

      /* check if the new score make the atom selected, and update the solution accordingly */
      if(prevScore > _instance.t && newScore < _instance.t) {
        change = true;
        _nAtomsBody[ruleIndex] -= 1;
        _sumAtomScores[ruleIndex] -= newScore-prevScore;
        _selectedAtoms[ruleIndex].erase(_atomIndexes[colInd][atomValue]);
      }

      /* update the atom positive error */
      if(atomValue == dataVal) {
        _atomErrors[ruleIndex][atomIndex].first += 1;
      }

    }

  }

  if(change) {
    _score = 0;
    for(int ruleInd = 0; ruleInd < _instance.p_rules; ruleInd ++) {

      // _score += _sumAtomScores[ruleInd]/static_cast<double>(_nAtomsBody[ruleInd]);

      double ruleScore = 0;
      for(int atomIndex : _selectedAtoms[ruleIndex]) {

        auto atom = _atoms[atomIndex];

        double atomScore = static_cast<double>(_atomErrors[ruleIndex][atomIndex].second)/static_cast<double>(_instance.negatives.size());
        atomScore -= static_cast<double>(_atomErrors[ruleIndex][atomIndex].first)/static_cast<double>(_nPositives[ruleIndex]);

        ruleScore += atomScore;

      }

      ruleScore /= static_cast<double>(_selectedAtoms.size());

      _score += ruleScore;

    }
  }

}

/*----------------------------------------------------------------------------*/
double Solution::score() {
  return _score;
}

/*----------------------------------------------------------------------------*/
void Solution::recomputeScore() {

  _score = 0;
  for(uint ruleIndex = 0; ruleIndex < _instance.p_rules; ruleIndex ++) {
    _score += computeRuleScore(ruleIndex);
  }

}

/*----------------------------------------------------------------------------*/
double Solution::computeRuleScore(int ruleIndex) {

  _nPositives[ruleIndex] = 0;

  /* compute the positive examples selected for this rule */
  std::vector<uint> positives;
  for(int indVar = 0; indVar < _var.size(); indVar ++) {
    if(_var[indVar] == ruleIndex) {
      positives.push_back(_instance.positives[indVar]);
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

  _sumAtomScores[ruleIndex] = 0;
  _nAtomsBody[ruleIndex] = 0;

  _selectedAtoms[ruleIndex].clear();

  for(uint atomInd = 0; atomInd < _atoms.size(); atomInd ++) {

    double score = static_cast<double>(_atomErrors[ruleIndex][atomInd].second)/static_cast<double>(_instance.negatives.size());
    score -= static_cast<double>(_atomErrors[ruleIndex][atomInd].first)/static_cast<double>(_nPositives[ruleIndex]);

    if(score >= _instance.t) {
      _nAtomsBody[ruleIndex] ++;
      _sumAtomScores[ruleIndex] += score;
      _selectedAtoms[ruleIndex].insert(atomInd);
    }

  }

  std::cout << _selectedAtoms[ruleIndex].size() << " atoms selected for rule " << ruleIndex << std::endl;

  double ruleScore = _sumAtomScores[ruleIndex] / static_cast<double>(_nAtomsBody[ruleIndex]);

  return ruleScore;

}
