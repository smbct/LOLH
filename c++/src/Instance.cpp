/*!
 * \file Instance.cpp
 * \author S. Buchet
 * \brief implementation of Instance class
 */


#include <fstream>

#include <sstream>

#include <algorithm>

#include "Instance.hpp"

#include "Combinatorics.hpp"

#include "Utils.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
Instance::Instance(const DataFrame<uint>& dataset, const vector<bool>& classVector):
_dataset(dataset)
{
  Combinatorics::generateRange(_dataset.nColumns(), _features);
  initialize(classVector);
}

/*----------------------------------------------------------------------------*/
Instance::Instance(const DataFrame<uint>& dataset, const vector<bool>& classVector, uint nval):
_dataset(dataset)
{
  Combinatorics::generateRange(_dataset.nColumns(), _features);
  initialize(classVector, nval);
}

/*----------------------------------------------------------------------------*/
Instance::Instance(const DataFrame<uint>& dataset, const std::vector<uint>& features, const std::vector<bool>& classVector, uint nval):
_dataset(dataset),
_features(features)
{
  initialize(classVector, nval);
}

/*----------------------------------------------------------------------------*/
Instance::Instance(const DataFrame<uint>& dataset, const std::vector<uint>& features, const std::vector<bool>& classVector):
_dataset(dataset),
_features(features)
{
  initialize(classVector);
}


/*----------------------------------------------------------------------------*/
Instance::Instance(const Instance& instance, const vector<uint>& excluded, uint nVal):
_dataset(instance._dataset)
{
  _features = instance._features;

  initializeSubInstance(instance, excluded, nVal);
}

/*----------------------------------------------------------------------------*/
void Instance::initialize(const vector<bool>& classVector) {

  for(uint index = 0; index < classVector.size(); index ++) {
    if(classVector[index]) {
      _positives.push_back(index);
    } else {
      _negatives.push_back(index);
    }
  }

  if(_positives.empty()) {
    cout << "WARNING: the instance should contain positive samples" << endl;
  }
  if(_negatives.empty()) {
    cout << "WARNING: the instance should contain negative samples" << endl;
  }

  /* number of variables */
  _nvar = _dataset.nColumns();

  /* init number of values for each variable */
  _nval.resize(_dataset.nColumns(), 0);

  /* indexing from the data */
  _posScore.resize(_dataset.nColumns());
  _negScore.resize(_dataset.nColumns());

  _atomIndexes.resize(_dataset.nColumns());


  /* lsit of counts for each value of each variable to compute the positive and negative scores */
  vector<vector<uint>> posCounts, negCounts;
  /* counting occurences in positive and negative samples */
  _dataset.uniqueCount(classVector, posCounts, negCounts);


  // for(uint indCol = 0; indCol < _dataset.nColumns(); indCol ++) {
  //
  //   /* initialization of the values and scores */
  //   _nval[indCol] = _dataset.nUnique(indCol);
  //
  //   _posScore[indCol].resize(_nval[indCol], 0);
  //   _negScore[indCol].resize(_nval[indCol], 0);
  //
  //   _atomIndexes[indCol].resize(_nval[indCol]);
  //
  //   for(uint indVal = 0; indVal < _nval[indCol]; indVal ++) {
  //
  //     /* create the list of atoms */
  //     _atoms.push_back(pair<uint,uint>(indCol, indVal));
  //     _atomIndexes[indCol][indVal] = static_cast<uint>(_atoms.size())-1;
  //
  //     /* initialization of the score for this atom */
  //     _posScore[indCol][indVal] = nPositives()-posCounts[indCol][indVal];
  //     _negScore[indCol][indVal] = nNegatives()-negCounts[indCol][indVal];
  //
  //   }
  //
  // }

  /* only use prediction features as atoms */
  for(auto& indCol : _features) {

    /* initialization of the values and scores */
    _nval[indCol] = _dataset.nUnique(indCol);

    _posScore[indCol].resize(_nval[indCol], 0);
    _negScore[indCol].resize(_nval[indCol], 0);

    _atomIndexes[indCol].resize(_nval[indCol]);

    for(uint indVal = 0; indVal < _nval[indCol]; indVal ++) {

      /* create the list of atoms */
      _atoms.push_back(pair<uint,uint>(indCol, indVal));
      _atomIndexes[indCol][indVal] = static_cast<uint>(_atoms.size())-1;

      /* initialization of the score for this atom */
      _posScore[indCol][indVal] = nPositives()-posCounts[indCol][indVal];
      _negScore[indCol][indVal] = nNegatives()-negCounts[indCol][indVal];

    }

  }

}

/*----------------------------------------------------------------------------*/
void Instance::initialize(const vector<bool>& classVector, uint nval) {

  /* _features.empty() means that all features are used */

  for(uint index = 0; index < classVector.size(); index ++) {
    if(classVector[index]) {
      _positives.push_back(index);
    } else {
      _negatives.push_back(index);
    }
  }

  if(_positives.empty()) {
    cout << "WARNING: the instance should contain positive samples" << endl;
  }
  if(_negatives.empty()) {
    cout << "WARNING: the instance should contain negative samples" << endl;
  }

  /* number of variables: if _features is empty, all variables are used */
  _nvar = _dataset.nColumns();
  _nval.resize(_dataset.nColumns(), nval);

  _atoms.resize(_nvar*nval);

  /* positive and negative scores for the atoms: size = dataset size but depending on the features everything is not used */
  _posScore.resize(_dataset.nColumns(), vector<uint>(nval, 0));
  _negScore.resize(_dataset.nColumns(), vector<uint>(nval, 0));

  /* same here: some atoms do not exist -> index=0 */
  _atomIndexes.resize(_dataset.nColumns(), vector<uint>(nval, 0));

  if(_features.size() == _nvar) {

    /* counting occurences in positive and negative samples */
    _dataset.uniqueCount(classVector, _posScore, _negScore);

    /* create all atoms of the instance */
    uint indAtom = 0;
    for(uint indCol = 0; indCol < _dataset.nColumns(); indCol ++) {
      for(uint indVal = 0; indVal < nval; indVal ++) {
        /* create the list of atoms */
        _atoms[indAtom] = pair<uint,uint>(indCol, indVal);
        _atomIndexes[indCol][indVal] = indAtom;
        indAtom += 1;

        /* initialization of the score for this atom */
        _posScore[indCol][indVal] = nPositives()-_posScore[indCol][indVal];
        _negScore[indCol][indVal] = nNegatives()-_negScore[indCol][indVal];
      }
    }

  } else { /* only a subset of the features are used */

    _dataset.uniqueCount(_features, classVector, _posScore, _negScore);

    /* creation of the atoms: only in the features list */
    uint indAtom = 0;
    for(uint& indCol: _features) {
      for(uint indVal = 0; indVal < nval; indVal ++) {
        /* create the list of atoms */
        _atoms[indAtom] = pair<uint,uint>(indCol, indVal);
        _atomIndexes[indCol][indVal] = indAtom;
        indAtom += 1;

        /* initialization of the score for this atom */
        _posScore[indCol][indVal] = nPositives()-_posScore[indCol][indVal];
        _negScore[indCol][indVal] = nNegatives()-_negScore[indCol][indVal];
      }
    }

  }

}

/*----------------------------------------------------------------------------*/
void Instance::initializeSubInstance(const Instance& instance, const vector<uint>& excluded, uint nValues) {

  // cout << "sub instance creation: " << endl;

  _nvar = instance._nvar;
  _nval = instance._nval;

  /* list of atoms are copied (same prediction features) */
  _atoms = instance._atoms;
  _atomIndexes = instance._atomIndexes;

  _features = instance._features;

  /* remove excluded positive samples */
  for(uint ind = 0; ind < instance._positives.size(); ind ++) {
    auto it = find(excluded.begin(), excluded.end(), ind);
    if(it == excluded.end()) {
      _positives.push_back(instance._positives[ind]);
    }
  }


  _negatives = instance._negatives;

  _posScore = instance._posScore;
  _negScore = instance._negScore;

  // vector<vector<uint>> counts(_nvar, vector<uint>(nValues,0));
  vector<vector<uint>> counts(_dataset.nColumns(), vector<uint>(nValues,0));

  /* list of sample indexes to exclude */
  vector<uint> excludedInd(excluded.size());
  for(uint ind = 0; ind < excluded.size(); ind ++) {
    excludedInd[ind] = instance._positives[excluded[ind]];
  }

  /* compute atom occurences in excluded positive samples */
  _dataset.uniqueCount(excludedInd, counts);

  if(_features.size() == _nvar) { /* all features of the dataset are used */

    /* compute the new scores considering the occurences to remove */
    for(uint indCol = 0; indCol < _nvar; indCol ++) {
      for(uint indVal = 0; indVal < nValues; indVal ++) {
        /* number of samples where value != indVal */
        uint nDiff = static_cast<uint>(excludedInd.size())-counts[indCol][indVal];
        _posScore[indCol][indVal] -= nDiff;
      }
    }

  } else { /* only a subset of the features are used in the dataset */

    /* compute the new scores considering the occurences to remove */
    for(uint& indCol: _features) {
      for(uint indVal = 0; indVal < nValues; indVal ++) {
        /* number of samples where value != indVal */
        uint nDiff = static_cast<uint>(excludedInd.size())-counts[indCol][indVal];
        _posScore[indCol][indVal] -= nDiff;

      }
    }

  }


  // cout << "n excluded: " << excluded.size() << endl;
  // cout << "original n positives: " << instance._positives.size() << endl;
  // cout << "n positives: " << _positives.size() << endl;


}

/*----------------------------------------------------------------------------*/
Instance::~Instance() {

}

/*----------------------------------------------------------------------------*/
uint Instance::getAtomIndex(Atom atom) const {
  return _atomIndexes[atom.first][atom.second];
}

/*----------------------------------------------------------------------------*/
Atom Instance::getAtom(uint index) const {
  return _atoms[index];
}

/*----------------------------------------------------------------------------*/
pair<std::string, uint> Instance::getAtomLabel(uint index) const {
  return pair<string, uint>(_dataset.getColLabel(_atoms[index].first), _atoms[index].second);
}


/*----------------------------------------------------------------------------*/
Score Instance::getAtomScore(uint index) const {
  return getAtomScore(_atoms[index]);
}

/*----------------------------------------------------------------------------*/
Score Instance::getAtomScore(Atom atom) const {
  return Score(_posScore[atom.first][atom.second], _negScore[atom.first][atom.second]);
}

/*----------------------------------------------------------------------------*/
uint Instance::getAtomPosScore(uint index) const {
  auto atom = _atoms[index];
  return getAtomPosScore(atom);
}

/*----------------------------------------------------------------------------*/
uint Instance::getAtomPosScore(Atom& atom) const {
  return _posScore[atom.first][atom.second];
}

/*----------------------------------------------------------------------------*/
uint Instance::getAtomNegScore(uint index) const {
  auto atom = _atoms[index];
  return getAtomNegScore(atom);
}

/*----------------------------------------------------------------------------*/
uint Instance::getAtomNegScore(Atom& atom) const {
  return _negScore[atom.first][atom.second];
}

/*----------------------------------------------------------------------------*/
uint Instance::matchingScore(const Body& body, uint indSample) const {
  uint score = 0;
  auto& row = _dataset.getRow(indSample);
  for(auto& ind : body) {
    if(row[_atoms[ind].first] != _atoms[ind].second) {
      score += 1;
    }
  }
  return score;
}

/*----------------------------------------------------------------------------*/
uint Instance::positiveMatchingScore(const Body& body, uint indPositive) const {
  return matchingScore(body, _positives[indPositive]);
}

/*----------------------------------------------------------------------------*/
uint Instance::negativeMatchingScore(const Body& body, uint indNegative) const {
  return matchingScore(body, _negatives[indNegative]);
}

/*----------------------------------------------------------------------------*/
const Sample& Instance::getSample(uint index) const {
  return _dataset.getRow(index);
}

/*----------------------------------------------------------------------------*/
const Sample& Instance::getPositive(uint ind) const {
  return _dataset.getRow(_positives[ind]);
}

/*----------------------------------------------------------------------------*/
const Sample& Instance::getNegative(uint ind) const {
  return _dataset.getRow(_negatives[ind]);
}

/*----------------------------------------------------------------------------*/
uint Instance::getPosClassIndex(uint indPos) const {
  return _positives[indPos];
}

/*----------------------------------------------------------------------------*/
uint Instance::getNegClassIndex(uint indNeg) const {
  return _negatives[indNeg];
}

/*----------------------------------------------------------------------------*/
uint Instance::nSamples() const {
  return _dataset.nRows();
}


/*----------------------------------------------------------------------------*/
uint Instance::nPositives() const {
  return static_cast<uint>(_positives.size());
}

/*----------------------------------------------------------------------------*/
uint Instance::nNegatives() const {
  return static_cast<uint>(_negatives.size());
}

/*----------------------------------------------------------------------------*/
uint Instance::nVariables() const {
  return _nvar;
}

/*----------------------------------------------------------------------------*/
uint Instance::nValues(uint varInd) const {
  return _nval.at(varInd);
}

/*----------------------------------------------------------------------------*/
uint Instance::nAtoms() const {
  return static_cast<uint>(_atoms.size());
}

/*----------------------------------------------------------------------------*/
uint Instance::nFeatures() const {
  return static_cast<uint>(_features.size());
}

/*----------------------------------------------------------------------------*/
const vector<uint>& Instance::getFeatures() const {
  return _features;
}

/*----------------------------------------------------------------------------*/
const DataFrame<uint>& Instance::getDataset() const {
  return _dataset;
}


/*----------------------------------------------------------------------------*/
string Instance::toString() const {

  string res = "";

  res += "n variables: " + to_string(_nvar) + "\n";
  res += "n atoms: " + to_string(_atoms.size()) + "\n";
  res += "n positive samples: " + to_string(_positives.size()) + "\n";
  res += "n negative samples: " + to_string(_negatives.size()) + "\n";

  return res;
}

/*----------------------------------------------------------------------------*/
uint Instance::initRandomClass(double proportion, vector<bool>& classVector) {

  uint nPos = static_cast<uint>(proportion*static_cast<double>(classVector.size()));

  /* make sure there are positive and negative samples */
  if(nPos == classVector.size()) {
    nPos --;
  } else if(nPos == 0) {
    nPos ++;
  }

  vector<uint> indexes;
  Combinatorics::generateRange(static_cast<uint>(classVector.size()), indexes);
  Utils::getInstance().shuffle(indexes);

  for(uint ind = 0; ind < nPos; ind ++) {
    classVector[indexes[ind]] = true;
  }

  return nPos;
}

/*----------------------------------------------------------------------------*/
uint Instance::initClusterClass(DataFrame<uint>& dataset, DataFrame<string>& labels, string cellLabel, vector<bool>& classVector) {

  uint nPos = 0;

  for(uint indRow = 0; indRow < labels.nRows(); indRow ++) {
    string label = labels.getRowLabel(indRow);

    if(labels.getData(indRow, 0) == cellLabel && dataset.containsRowIndex(label)) {
      classVector[dataset.getRowIndex(label)] = true;
      nPos ++;
    }

  }

  return nPos;
}

/*----------------------------------------------------------------------------*/
uint Instance::initCoexprClass(DataFrame<uint>& dataset, uint geneIndex, uint value, vector<uint>& features, vector<bool>& classVector) {

  uint nPos = 0;

  /* class vector: true iif gene value is the target value */
  for(uint indRow = 0; indRow < classVector.size(); indRow ++) {
    auto& sample = dataset.getRow(indRow);
    if(sample[geneIndex] == value) {
      classVector[indRow] = true;
      nPos ++;
    }
  }

  /* list of features: all features except targetFeature */
  features.resize(dataset.nColumns()-1);
  uint indFeature = 0;
  for(uint indCol = 0; indCol < dataset.nColumns(); indCol ++) {
    if(indCol != geneIndex) {
      features[indFeature] = indCol;
      indFeature ++;
    }
  }

  return nPos;
}

/*----------------------------------------------------------------------------*/
uint Instance::initRegulationClass(DataFrame<uint>& dataset, NGraph& graph, uint geneIndex, uint value, double proportion, uint predNeq, vector<bool>& classVector) {

  uint nPos = 0;

  for(uint indCell = 0; indCell < dataset.nRows(); indCell ++) {

    if(!graph[indCell].empty() && (predNeq == 0 || (predNeq == 1 && dataset.getRow(indCell)[geneIndex] != value) || (predNeq == 2 && dataset.getRow(indCell)[geneIndex] == value) ) ) { /* look at all successors */

      uint occurences = 0;
      for(auto& indSuc : graph[indCell]) {
        if(dataset.getRow(indSuc)[geneIndex] == value) {
          occurences ++;
        }
      }



      /* if the target atom is detected in enough successors, the sample is considered as positive */
      if( (proportion <= -0.5 and occurences > 0) || (static_cast<double>(occurences) / static_cast<double>(graph[indCell].size()) >= proportion) ) {
        classVector[indCell] = true;
        nPos ++;
      }

      // if(occurences > 0) {
      //   classVector[indCell] = true;
      //   nPos ++;
      // }

    }
  }

  return nPos;

}

/*----------------------------------------------------------------------------*/
ostream& operator<<(ostream& str, const Instance& instance) {
  str << instance.toString();
  return str;
}


/*----------------------------------------------------------------------------*/
ostream& operator<<(ostream& str, const pair<uint,uint>& tuple) {
  str << "(" << tuple.first << "," << tuple.second << ")";
  return str;
}

/*----------------------------------------------------------------------------*/
ostream& operator<<(ostream& str, const Body& body) {
  str << "[";
  for(uint ind = 0; ind < body.size(); ind ++) {
    str << body[ind];
    if(ind < body.size()-1) {
      str << ", ";
    }
  }
  str << "]";
  return str;
}

/*----------------------------------------------------------------------------*/
string bodyToString(const Body& body) {
  ostringstream str;
  str << body;
  return str.str();
}
