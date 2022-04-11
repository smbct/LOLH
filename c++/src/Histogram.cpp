/*!
 * \file Histogram.cpp
 * \author S. Buchet
 * \brief implementation of class Histogram
 */

#include "Histogram.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
Histogram::Histogram(const Instance& instance, const Body& solution, Type type):
_instance(instance),
_solution(solution),
_type(type)
{

  /* build the histograms */

  if(type != Type::POSITIVE) {
    buildNegative();
  }

  if(type != Type::NEGATIVE) {
    buildPositive();
  }

}

/*----------------------------------------------------------------------------*/
void Histogram::buildPositive() {
  _positives.resize(_solution.size()+1);
  for(uint indPos = 0; indPos < _instance.nPositives(); indPos ++) {
    _positives[_instance.positiveMatchingScore(_solution, indPos)].push_back(indPos);
  }
}

/*----------------------------------------------------------------------------*/
void Histogram::buildNegative() {
  _negatives.resize(_solution.size()+1);
  for(uint indNeg = 0; indNeg < _instance.nNegatives(); indNeg ++) {
    _negatives[_instance.negativeMatchingScore(_solution, indNeg)].push_back(indNeg);
  }
}

/*----------------------------------------------------------------------------*/
Histogram::Type Histogram::type() const {
  return _type;
}

/*----------------------------------------------------------------------------*/
uint Histogram::positiveCovered(uint threshold) {
  uint covering = 0;
  for(uint error = 0; error <= threshold; error ++) {
    covering += static_cast<uint>(_positives[error].size());
  }
  return covering;
}

/*----------------------------------------------------------------------------*/
uint Histogram::negativeCovered(uint threshold) {
  uint covering = 0;
  for(uint error = 0; error <= threshold; error ++) {
    covering += static_cast<uint>(_negatives[error].size());
  }
  return covering;
}

/*----------------------------------------------------------------------------*/
void Histogram::getPositiveCovered(uint threshold, vector<uint>& indexes) {
  uint ncov = 0;
  for(uint error = 0; error <= threshold; error ++) {
    ncov += static_cast<uint>(_positives[error].size());
  }
  indexes.resize(ncov);
  uint ind = 0;
  for(uint error = 0; error <= threshold; error ++) {
    for(uint ind2 = 0; ind2 < _positives[error].size(); ind2 ++) {
      indexes[ind] = _positives[error][ind2];
      ind += 1;
    }
  }
}

/*----------------------------------------------------------------------------*/
const vector<uint>& Histogram::getPosSampleScore(uint score) {
  return _positives[score];
}

/*----------------------------------------------------------------------------*/
const vector<uint>& Histogram::getNegSampleScore(uint score) {
  return _negatives[score];
}

/*----------------------------------------------------------------------------*/
ostream& Histogram::display1(ostream& str) const {
  /* length of the maximum bar */
  uint barLength = 75;

  /* display the negative histogram */
  if(_type != Histogram::Type::NEGATIVE) {
    str << "POSITIVE HISTOGRAM" << endl;
    str << "n positive samples: " << _instance.nPositives() << endl;
    for(uint error = 0; error < _positives.size(); error ++) {
      if(error < 10) {
        str << "0";
      }
      str << error << ": ";
      uint nSymb = static_cast<uint>(_positives[error].size()*barLength)/_instance.nPositives();
      for(uint ind = 1; ind <= nSymb; ind ++) {
        str << "#";
      }
      str << endl;
    }
  }

  str << endl;

  /* display the positive histogram */
  if(_type != Histogram::Type::POSITIVE) {
    str << "NEGATIVE HISTOGRAM" << endl;
    str << "n negative samples: " << _instance.nPositives() << endl;
    for(uint error = 0; error < _negatives.size(); error ++) {
      if(error < 10) {
        str << "0";
      }
      str << error << ": ";
      uint nSymb = static_cast<uint>(_negatives[error].size()*barLength)/_instance.nNegatives();
      for(uint ind = 1; ind <= nSymb; ind ++) {
        str << "#";
      }
      str << endl;
    }

  }

  return str;
}

/*----------------------------------------------------------------------------*/
ostream& Histogram::display2(ostream& str) const {

  uint barLength = 50;

  vector<uint> posHeights, negHeights;

  uint maxPos = 0, maxNeg = 0;

  if(_type != Type::NEGATIVE) {
    posHeights.resize(_positives.size());
    for(uint error = 0; error < posHeights.size(); error ++) {
      posHeights[error] = static_cast<uint>(_positives[error].size()*barLength)/_instance.nPositives();
      if(posHeights[error] > maxPos) {
        maxPos = posHeights[error];
      }
    }
  }

  if(_type != Type::POSITIVE) {
    negHeights.resize(_negatives.size());
    for(uint error = 0; error < negHeights.size(); error ++) {
      negHeights[error] = static_cast<uint>(_negatives[error].size()*barLength)/_instance.nNegatives();
      if(negHeights[error] > maxNeg) {
        maxNeg = negHeights[error];
      }
    }
  }

  uint maxHeight = max(maxPos, maxNeg);

  /* display the bars */
  for(uint line = maxHeight+1; line >= 1; line --) {
    if(_type != Type::NEGATIVE) {
      for(uint error = 0; error < posHeights.size(); error ++) {
        if(posHeights[error] >= line) {
          str << "## ";
        } else {
          str << "   ";
        }
      }
    }
    if(_type == Type::POSITIVE_NEGATIVE) {
      str << "   |   ";
    }
    if(_type != Type::POSITIVE) {
      for(uint error = 0; error < negHeights.size(); error ++) {
        if(negHeights[error] >= line) {
          str << "## ";
        } else {
          str << "   ";
        }
      }
    }
    str << endl;
  }

  /* display the axis */
  if(_type != Type::NEGATIVE) {
    for(uint error = 0; error < posHeights.size(); error ++) {
      if(error < 10) {
        str << "0";
      }
      str << error << " ";
    }
  }

  if(_type == Type::POSITIVE_NEGATIVE) {
    str << "   |   ";
  }

  if(_type != Type::POSITIVE) {
    for(uint error = 0; error < negHeights.size(); error ++) {
      if(error < 10) {
        str << "0";
      }
      str << error << " ";
    }
    str << endl;
  }

  /* histogram names */
  if(_type != Type::NEGATIVE) {
    str << "POSITIVE HISTOGRAM";
    for(uint ind = 1; ind <= posHeights.size()*2 + posHeights.size()-1 - 17; ind ++) {
      str << " ";
    }
  }

  if(_type == Type::POSITIVE_NEGATIVE) {
    str << "   |   ";
  }

  if(_type != Type::POSITIVE) {
    str << "NEGATIVE HISTOGRAM";
  }


  return str;
}


/*----------------------------------------------------------------------------*/
ostream& operator<<(ostream& str, const Histogram& histogram) {
  return histogram.display2(str);
}
