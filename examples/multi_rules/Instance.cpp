#include "Instance.hpp"

/*-----------------------------------------------------------------------------*/
double LocalSearch::Instance::computeAtomScore(std::pair<int,int> atom) {

  int positiveError = 0;
  int negativeError = 0;

  for(auto positive : positives) {
    int dataVal = dataset.getData(positive, atom.first);
    if(dataVal != atom.second) {
      positiveError += 1;
    }
  }

  for(auto negative : negatives) {
    int dataVal = dataset.getData(negative, atom.first);
    if(dataVal != atom.second) {
      negativeError += 1;
    }
  }

  double score = static_cast<double>(negativeError) / static_cast<double>(negatives.size());
  score -= static_cast<double>(positiveError) / static_cast<double>(positives.size());

  return score;
}
