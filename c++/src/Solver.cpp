/*!
 * \file Solver.cpp
 * \author S. Buchet
 * \brief implementation of class Solver
 */

#include "Solver.hpp"

#include <cmath>
#include <algorithm>
#include <stack>

#include "Utils.hpp"

#include "Histogram.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
Solver::Solver(const Instance& instance):
_instance(instance)
{

}


/*----------------------------------------------------------------------------*/
pair<uint,uint> Solver::selectEquivalentScores(const vector<uint>& indexes, uint bodyLength) {
  pair<uint,uint> res(bodyLength-1, bodyLength-1);
  while(res.first > 0 && _instance.getAtomScore(indexes[res.first-1]) == _instance.getAtomScore(indexes[bodyLength-1])) {
    res.first -= 1;
  }
  while(res.second < indexes.size()-1 && _instance.getAtomScore(indexes[res.second+1]) == _instance.getAtomScore(indexes[bodyLength-1])) {
    res.second += 1;
  }
  return res;
}



/*----------------------------------------------------------------------------*/
pair<uint,uint> Solver::selectEquivalentScores(const vector<uint>& indexes, uint bodyLength, const vector<int>& score) {
  pair<uint,uint> res(bodyLength-1, bodyLength-1);
  while(res.first > 0 && score[indexes[res.first-1]] == score[indexes[bodyLength-1]]) {
    res.first -= 1;
  }
  while(res.second < indexes.size()-1 && score[indexes.at(res.second+1)] == score[indexes.at(bodyLength-1)]) {
    res.second += 1;
  }
  return res;
}



/*----------------------------------------------------------------------------*/
void Solver::createIndexList(vector<uint>& indexes, Objective obj) {
  indexes.resize(_instance.nFeatures());
  uint indFeature = 0;
  for(const auto& indVar : _instance.getFeatures()) {
    uint bestVal = 0;
    /* select the best atom from the current variable */
    for(uint value = 0; value < _instance.nValues(indVar); value ++) {
      if(value > 0) {
        auto score = _instance.getAtomScore(Atom(indVar,value));
        auto bestScore = _instance.getAtomScore(Atom(indVar,bestVal));
        if(obj == Objective::POSITIVE) {
          if(score.first < bestScore.first || ( (score.first == bestScore.first) && score.second > bestScore.second) ) {
            bestVal = value;
          }
        } else {
          if(score.second > bestScore.second || ( (score.second == bestScore.second) && score.first < bestScore.first) ) {
            bestVal = value;
          }
        }
      }
    }
    indexes[indFeature] = _instance.getAtomIndex(Atom(indVar, bestVal));
    indFeature ++;
  }
}



/*----------------------------------------------------------------------------*/
void Solver::createIndexList(vector<uint>& indexes, const vector<int>& score) {

  indexes.resize(_instance.nFeatures());
  uint indFeature = 0;
  for(uint indVar : _instance.getFeatures()) {
    uint bestVal = 0;
    /* select the best atom from the current variable */
    for(uint value = 0; value < _instance.nValues(indVar); value ++) {
      if(value > 0) {
        if(score[_instance.getAtomIndex(Atom(indVar, value))] < score[_instance.getAtomIndex(Atom(indVar, bestVal))]) {
          bestVal = value;
        }
      }
    }
    indexes[indFeature] = _instance.getAtomIndex(Atom(indVar, bestVal));
    indFeature ++;
  }
}


/*----------------------------------------------------------------------------*/
void Solver::createAssociationList(vector<uint>& indexes, pair<uint,uint> equiv, vector<vector<uint>>& result) {
   /* map between bi-objective score and index in the result vector */
  map<Score,uint> values;
  for(uint ind = equiv.first; ind <= equiv.second; ind ++) {
    auto score = _instance.getAtomScore(indexes[ind]);
    if(values.count(score) > 0) {
      result[values[score]].push_back(indexes[ind]);
    } else {
      values[score] = static_cast<int>(result.size());
      result.push_back(vector<uint>(1,indexes[ind]));
    }
  }
}


/*----------------------------------------------------------------------------*/
Score Solver::solveMono(uint bodyLength, Objective obj, Body* solution) {


  vector<uint> indexes; /* list of atom indexes */

  // indexes.resize(_instance.nAtoms());
  // for(uint ind = 0; ind  < _instance.nAtoms(); ind ++) {
  //   indexes[ind] = ind;
  // }

  /* create the list of index with no duplication from the same variable */
  createIndexList(indexes, obj);

  /* select the right objective to optimize */
  if(obj == Objective::POSITIVE) {
    auto comp = [this](uint ind1, uint ind2)
    { return _instance.getAtomPosScore(ind1) < _instance.getAtomPosScore(ind2) || ( (_instance.getAtomPosScore(ind1) == _instance.getAtomPosScore(ind2)) && _instance.getAtomNegScore(ind1) > _instance.getAtomNegScore(ind2)); };
    sort(indexes.begin(), indexes.end(), comp);


  } else { /* sort along the negative objective */
    auto comp = [this](uint ind1, uint ind2)
    { return _instance.getAtomNegScore(ind1) > _instance.getAtomNegScore(ind2) || ( (_instance.getAtomNegScore(ind1) == _instance.getAtomNegScore(ind2)) && _instance.getAtomPosScore(ind1) < _instance.getAtomPosScore(ind2)); };
    sort(indexes.begin(), indexes.end(), comp);
  }



  /* compute all possible subsets of elements having the same score */
  /*pair<uint,uint> equiv = selectEquivalentScores(indexes, bodyLength);
  vector<uint> elements(equiv.second-equiv.first);
  for(uint ind = 0; ind < elements.size(); ind ++) {
    elements[ind] = ind;
  }

  vector<vector<uint>> subsets;

  Combinatorics::choose_kn(elements, bodyLength-equiv.first, subsets);

  for(auto& subset : subsets) {
    for(uint ind = 0; ind < subset.size(); ind ++) {
      cout << subset[ind] << " ";
    }
    cout << endl;
  }*/

  /* compute one optimal solution */
  if(solution != nullptr) {
    solution->assign(indexes.begin(), indexes.begin()+bodyLength);
  }


  /* compute optimal value */
  Score score(0,0);
  for(uint ind = 0; ind < bodyLength; ind ++) {
    auto atomScore = _instance.getAtomScore(indexes[ind]);
    score.first += atomScore.first;
    score.second += atomScore.second;
  }

  // cout << "final score: (" << score.first << ", " << score.second << ")" << endl;

  return score;

}



/*----------------------------------------------------------------------------*/
void Solver::solveScalarization(uint bodyLength, Weight weight, vector<Score>* points, vector<Body>* solutions) {

  // cout << "solving " << weight.first << "*pos_score + " << weight.second << "*neg_score" << endl;

  /* compute the scalarization */
  vector<int> score(_instance.nAtoms());
  for(uint atomInd = 0; atomInd < _instance.nAtoms(); atomInd ++) {
    auto atomScore = _instance.getAtomScore(atomInd);
    score[atomInd] = weight.first * atomScore.first + weight.second * atomScore.second;
  }


  vector<uint> indexes;

  /* create a list of all atom indexes */
  /*indexes.resize(_instance.nAtoms());
  for(uint ind = 0; ind  < _instance.nAtoms(); ind ++) {
    indexes[ind] = ind;
  }*/

  createIndexList(indexes, score); /* create the list of index with no duplication from the same variable */

  /* minimization of the scalarization */
  sort(indexes.begin(), indexes.end(), [&score](uint ind1, uint ind2) { return score[ind1] < score[ind2]; }); /* compute atom ordering based on the score */

  /* compute tie elements */
  pair<uint,uint> equiv = selectEquivalentScores(indexes, bodyLength, score);

  // cout << "equivalent elements: " << equiv << endl;

  /* little hack to prevent the algorithm spending too much time on atoms that are not interesting at all */
  if(equiv.second <= bodyLength-1 || equiv.second-equiv.first >= 10) { /* all atoms selected have a different score */


    /* compute the unique point */
    Score point(0,0);
    for(uint ind = 0; ind < bodyLength; ind ++) {
      point.first += _instance.getAtomPosScore(indexes[ind]);
      point.second += _instance.getAtomNegScore(indexes[ind]);
    }
    points->push_back(point);

    /* compute the unique body of solutions is not null */
    if(solutions != nullptr) {
      solutions->push_back(Body(indexes.begin(), indexes.begin()+bodyLength));
    }


  } else { /* some atoms have the same scalarized score */

    /* assemble atoms with the same bi-objective score */
    vector<vector<uint>> associations;
    createAssociationList(indexes, equiv, associations);

    /* create a list of occurences to compute al possible bi-objective points */
    vector<uint> occurences(associations.size());
    for(uint ind = 0; ind < occurences.size(); ind ++) {
      occurences[ind] = static_cast<uint>(associations[ind].size());
    }

    /* compute all the different points that can be obtained */
    vector<vector<pair<uint,uint>>> subsets;
    // cout << "start kn occurencess" << endl;
    Combinatorics::choose_kn_occurences(occurences, bodyLength-equiv.first, subsets);
    // cout << "end kn occurencess" << endl;


    // cout << "n subsets: " << subsets.size() << endl;

    /* compute base score */
    Score baseScore(0,0);
    for(uint ind = 0; ind < equiv.first; ind ++) {
      baseScore.first += _instance.getAtomPosScore(indexes[ind]);
      baseScore.second += _instance.getAtomNegScore(indexes[ind]);
    }

    // cout << "n subsets: " << subsets.size() << endl;

    /* compute all sub scores */
    for(auto& subset : subsets) {

      /* compute the corresponding point */
      Score subScore(0,0);
      for(uint ind = 0; ind < subset.size(); ind ++) {
        auto atomScore = _instance.getAtomScore(associations[subset[ind].first].front());
        subScore.first += atomScore.first*static_cast<int>(subset[ind].second);
        subScore.second += atomScore.second*static_cast<int>(subset[ind].second);
      }
      points->push_back(Score(baseScore.first+subScore.first, baseScore.second+subScore.second));

      /* if required compute the first corresponding body */
      if(solutions != nullptr) {
        solutions->push_back(Body());
        solutions->back().insert(solutions->back().end(), indexes.begin(), indexes.begin()+equiv.first); /* insert tha base body */
        for(uint ind = 0; ind < subset.size(); ind ++) {
          solutions->back().insert(solutions->back().end(), associations[subset[ind].first].begin(), associations[subset[ind].first].begin()+subset[ind].second);
        }
      }
    }


  }

  // cout << "time sort: " << timeSort << endl;

}

/*----------------------------------------------------------------------------*/
void Solver::filterPoints(Score left, Score right, vector<Score>* points, vector<Body>* solutions) {

  /* create a list of indexes, excluding left and right if necessary */
  vector<uint> indexes;
  for(uint ind = 0; ind < points->size(); ind ++) {
    Score point = (*points)[ind];
    if(point.first > left.first && point.first < right.first && point.second > left.second && point.second < right.second) {
      indexes.push_back(ind);
    }
  }

  if(indexes.size() > 1) { /* no need for sorting on one point */
    /* points are sorted */
    sort(indexes.begin(), indexes.end(), [&points](uint ind1, uint ind2) { return (*points)[ind1].first < (*points)[ind2].first; });

    auto end = unique(indexes.begin(), indexes.end(), [points](uint ind1, uint ind2) {return (*points)[ind1].first == (*points)[ind2].first; });
    indexes.resize(std::distance(indexes.begin(), end));
  }

  /* rearrenge the final list of points */
  Combinatorics::rearrenge(*points, indexes);

  if(solutions != nullptr) { /* if necessary, initialize the final list of solutions */
    Combinatorics::rearrenge(*solutions, indexes);
  }

}


/*----------------------------------------------------------------------------*/
void Solver::computeSupported(uint bodyLength, vector<Score>* points, vector<Weight>* weights, vector<Body>* solutions) {

  bool computeSolutions = (solutions != nullptr);
  bool computeWeight = (weights != nullptr);

  stack<pair<Score,Score>> pending;

  /* solve the mono objective problems */
  Score pos, neg;


  if(computeSolutions) {

    solutions->push_back(Body());
    pos = solveMono(bodyLength, Objective::POSITIVE, &solutions->back());

    solutions->push_back(Body());
    neg = solveMono(bodyLength, Objective::NEGATIVE, &solutions->back());

    if(pos == neg) {
      solutions->pop_back();
    }
  } else {
    pos = solveMono(bodyLength, Objective::POSITIVE);
    neg = solveMono(bodyLength, Objective::NEGATIVE);
  }

  if(pos != neg) {
    points->resize(points->size()+2);
    (*points)[points->size()-2] = pos;
    (*points)[points->size()-1] = neg;
  } else {
    points->resize(points->size()+1);
    (*points)[points->size()-1] = pos;
  }

  if(computeWeight) {
    if(pos != neg) {
      weights->push_back(Weight(1,0));
      weights->push_back(Weight(0,-1));
    } else {
      weights->push_back(Weight(1,0));
    }
  }



  // cout << "pos: " << pos << endl;
  // cout << "neg: " << neg << endl;

  /* make sure the two points are different */

  // cout << "pos: " << pos << endl;
  // cout << "neg: " << neg << endl;
  // cout << endl;

  /* First points: left - right */
	if(pos != neg) {
		pending.push(pair<Score,Score>(pos,neg));
	}


  while(!pending.empty()) {

    /* get the top element and clear it */
    auto top = pending.top();
    pending.pop();

    // cout << "pending top: " << top.first << " ; " << top.second << endl;

    /* compute the weights */
    Weight weight(static_cast<int>(top.second.second)-static_cast<int>(top.first.second), static_cast<int>(top.first.first)-static_cast<int>(top.second.first));

    // cout << "weights: " << static_cast<int>(top.second.second)-static_cast<int>(top.first.second) << " ; " << static_cast<int>(top.first.first)-static_cast<int>(top.second.first) << endl;

    vector<Score> tempPoints;
    vector<Body> tempSolutions; /* optional */

    // cout << "begining scalarization" << endl;

    if(computeSolutions) {
      solveScalarization(bodyLength, weight, &tempPoints, &tempSolutions);
      /* remove points outside [first, second] and order them */
      filterPoints(top.first, top.second, &tempPoints, &tempSolutions);
    } else {
      solveScalarization(bodyLength, weight, &tempPoints, nullptr);
      filterPoints(top.first, top.second, &tempPoints, nullptr);
    }

    // cout << "end scalarization" << endl;


    /* insert the newly computed points */
    points->insert(points->end(), tempPoints.begin(), tempPoints.end());
    if(computeSolutions) {
      solutions->insert(solutions->end(), tempSolutions.begin(), tempSolutions.end());
    }
    if(computeWeight) { /* add one weigth per point (all points have the same weight) */
      weights->resize(weights->size() + tempPoints.size());
      for(uint ind = static_cast<uint>(weights->size()-tempPoints.size()); ind < weights->size(); ind ++) {
        (*weights)[ind] = weight;
      }
    }

    if(!tempPoints.empty()) { /* add the new scalarizations to the stack */
      pending.push(pair<Score,Score>(top.first, tempPoints.front()));
      for(int ind = 0; ind < static_cast<int>(tempPoints.size()-1); ind ++) {
        pending.push(pair<Score,Score>(tempPoints[ind], tempPoints[ind+1]));
      }
      pending.push(pair<Score,Score>(tempPoints.back(), top.second));
    }

  }

  /* sort on indexes and then modify all computed vectors */
  vector<uint> indexes;
  Combinatorics::generateRange(static_cast<uint>(points->size()), indexes);

  sort(indexes.begin(), indexes.end(), [&points] (uint ind1, uint ind2) { return (*points)[ind1].first < (*points)[ind2].first; } );

  Combinatorics::rearrenge(*points, indexes);
  if(computeWeight) {
    Combinatorics::rearrenge(*weights, indexes);
  }
  if(computeSolutions) {
    Combinatorics::rearrenge(*solutions, indexes);
  }

}

/*----------------------------------------------------------------------------*/
bool Solver::computeTargetSolution(uint bodyLength, uint targetCover, uint threshold, Score& targetPoint, Weight& targetWeight, Body& targetSolution) {

  bool res = true;

  Body solPos, solNeg;
  Score pos = solveMono(bodyLength, Objective::POSITIVE, &solPos);
  Score neg = solveMono(bodyLength, Objective::NEGATIVE, &solNeg);

  /* check if the target is feasible */
  Histogram histoPos(_instance, solPos, Histogram::Type::POSITIVE);
  if(histoPos.positiveCovered(threshold) < targetCover) {
    // cout << "goal is impossible" << endl;
    res = false;
  }

  /* first initialization: pos is the "worst" candidate point */
  targetPoint = pos;
  targetWeight = Weight(1,0);
  targetSolution = solPos;

  /* make sure there is not only one solution */
  if(res && pos.first != neg.first) {


    // cout << "pos: " << pos << endl;
    // cout << "neg: " << neg << endl;
    // cout << endl;

    /* First points: left - right */
    bool stop = false;

    /* two points define the scalarization to solve */
    pair<Score, Score> bounds(pos,neg);

    while(!stop) {

      /* compute the weights */
      Weight weight(bounds.second.second-bounds.first.second, bounds.first.first-bounds.second.first);

      vector<Score> tempPoints;
      vector<Body> tempSolutions;

      solveScalarization(bodyLength, weight, &tempPoints, &tempSolutions);

      /* remove points outside [first, second] and order them */
      filterPoints(bounds.first, bounds.second, &tempPoints, &tempSolutions);

      if(!tempPoints.empty()) {

        /* bound.first is guaranteed to verify the covering */
        /* look for the point with the larger index that also verify the target */

        /* look for the best point verifying the covering rate */
        int index = static_cast<int>(tempPoints.size()-1);
        bool found = false;

        while(!found && index >= 0) {
          const Body& solTemp = tempSolutions[index];
          Histogram histogram(_instance, solTemp, Histogram::Type::POSITIVE);
          uint cover = histogram.positiveCovered(threshold);

          // cout << "cover: " << cover << endl;

          if(cover >= targetCover) {
            found = true;
            // cout << "best cover: " << cover << endl;
          } else {
            index -= 1;
          }
        }

        if(found) { /* better candidate found */
          /* compute the new sub problem with a better candaidate */
          bounds = pair<Score, Score>(tempPoints[index], index+1 < static_cast<int>(tempPoints.size()) ? tempPoints[index+1] : bounds.second);
          targetPoint = tempPoints[index];
          targetWeight = weight;
          targetSolution = tempSolutions[index];
        } else { /* no better point has been found */
          bounds = pair<Score, Score>(bounds.first, tempPoints.front());
        }

      } else { /* the "almost optimal" supported is found */
        stop = true;
      }

    }

  } else {
    /* ill conditioned instance */
  }


  return res;
}


/*----------------------------------------------------------------------------*/
bool Solver::computeTargetSolutionSlow(uint bodyLength, uint targetCover, uint threshold, Score& targetPoint, Weight& targetWeight, Body& targetSolution) {

  bool res = true;

  Body solPos, solNeg;
  Score pos = solveMono(bodyLength, Objective::POSITIVE, &solPos);
  Score neg = solveMono(bodyLength, Objective::NEGATIVE, &solNeg);

  /* check if the target is feasible */
  Histogram histoPos(_instance, solPos, Histogram::Type::POSITIVE);
  if(histoPos.positiveCovered(threshold) < targetCover) {
    // cout << "goal is impossible" << endl;
    res = false;
  }

  /* first initialization: pos is the "worst" candidate point */
  targetPoint = pos;
  targetWeight = Weight(1,0);
  targetSolution = solPos;

  /* make sure there is not only one solution */
  if(res && pos.first != neg.first) {

    /* compute all the supported solutions */
    vector<Score> points;
    vector<Weight> weights;
    vector<Body> solutions;
    computeSupported(bodyLength, &points, &weights, &solutions);

    uint bestSol = 0;
    uint nNegBest = 0;

    /* iterate and compute the histograms -> select the best rule verifying the goal */
    for(uint indSol = 0; indSol < solutions.size(); indSol ++) {

        Histogram histo(_instance, solutions[indSol], Histogram::Type::POSITIVE_NEGATIVE);
        if(histo.positiveCovered(threshold) >= targetCover) {
          uint nNeg = histo.negativeCovered(threshold);
          if(indSol == 0 || nNeg < nNegBest) {
            bestSol = indSol;
            nNegBest = nNeg;
          }
        }

    }

    targetPoint = points[bestSol];
    targetWeight = weights[bestSol];
    targetSolution = solutions[bestSol];


  } else {
    /* ill conditioned instance */
  }


  return res;
}

/*----------------------------------------------------------------------------*/
void Solver::computeBestAtoms(vector<uint>& sortedAtoms, vector<double>& sortedScores) {

  Score leftAtomScore(0,0);
  Score rightAtomScore(_instance.nPositives(), _instance.nNegatives());

  /* build "diagonal" equation */
  double a = 1.;
  double b = (static_cast<double>(leftAtomScore.first)-rightAtomScore.first)/(static_cast<double>(rightAtomScore.second)-leftAtomScore.second);
  double c = -a*leftAtomScore.first - b*leftAtomScore.second;

  /* farthest point possible used to compute relative distances */
  double maxDist = b*static_cast<double>(_instance.nNegatives()) + c;

  /* compute all the relative distances */
  vector<double> relativeDistances(_instance.nAtoms(), 0.);
  for(uint indAtom = 0; indAtom < _instance.nAtoms(); indAtom ++) {
    Score score = _instance.getAtomScore(indAtom);
    relativeDistances[indAtom] = a*score.first + b*score.second + c;
    relativeDistances[indAtom] = relativeDistances[indAtom]/maxDist;
  }

  /* sort all the atoms based on the distances */
  Combinatorics::generateRange(_instance.nAtoms(), sortedAtoms);
  sort(sortedAtoms.begin(), sortedAtoms.end(), [&relativeDistances](uint ind1, uint ind2) { return relativeDistances[ind1] > relativeDistances[ind2]; });

  /* create a list of sorted scores (relative distances) */
  sortedScores.resize(sortedAtoms.size());
  for(uint ind = 0; ind < sortedScores.size(); ind ++) {
    sortedScores[ind] = relativeDistances[sortedAtoms[ind]];
  }

}

/*----------------------------------------------------------------------------*/
void Solver::computekBestAtoms(uint nAtoms, vector<uint>& selectedAtoms, vector<double>& atomScores) {

  vector<uint> sortedAtoms;
  vector<double> sortedScores;
  computeBestAtoms(sortedAtoms, sortedScores);

  selectedAtoms.assign(sortedAtoms.begin(), sortedAtoms.begin()+nAtoms);
  atomScores.assign(sortedScores.begin(), sortedScores.begin()+nAtoms);

}

/*----------------------------------------------------------------------------*/
void Solver::computeBestAtomsThreshold(double threshold, vector<uint>& selectedAtoms, vector<double>& atomScores) {

  vector<uint> sortedAtoms;
  vector<double> sortedScores;
  computeBestAtoms(sortedAtoms, sortedScores);

  /* compute the index */
  uint ind = 0;
  bool stop = false;
  while(!stop && ind < sortedAtoms.size()) {
    if(sortedScores[ind] >= threshold) {
      ind += 1;
    } else {
      stop = true;
    }
  }

  selectedAtoms.assign(sortedAtoms.begin(), sortedAtoms.begin()+ind);
  atomScores.assign(sortedScores.begin(), sortedScores.begin()+ind);

}


/*----------------------------------------------------------------------------*/
void Solver::computeNonDominatedAtoms(const vector<uint>& atoms, vector<uint>& nondominatedAtoms) {

  vector<uint> sortedAtoms = atoms;

  /* atoms are sorted: increasing in positive scores, decreasing in negative score */

  auto comp = [this](uint ind1, uint ind2)
  { return _instance.getAtomPosScore(ind1) < _instance.getAtomPosScore(ind2) || ( (_instance.getAtomPosScore(ind1) == _instance.getAtomPosScore(ind2)) && _instance.getAtomNegScore(ind1) > _instance.getAtomNegScore(ind2)); };

  sort(sortedAtoms.begin(), sortedAtoms.end(), comp);


  nondominatedAtoms.push_back(sortedAtoms.front());

  for(uint ind = 1; ind < sortedAtoms.size(); ind ++) {
    uint indAtomInsert = sortedAtoms[ind];

    uint indAtomLast = nondominatedAtoms.back();

    bool dominated = true;

    if(_instance.getAtomPosScore(indAtomLast) == _instance.getAtomPosScore(indAtomInsert) && _instance.getAtomNegScore(indAtomLast) == _instance.getAtomNegScore(indAtomInsert)) {
      dominated = false;
    }

    if(_instance.getAtomPosScore(indAtomLast) < _instance.getAtomPosScore(indAtomInsert) && _instance.getAtomNegScore(indAtomLast) < _instance.getAtomNegScore(indAtomInsert)) {
      dominated = false;
    }

    /* not dominated by its predecessor -> not dominated, can be inserted */
    if(!dominated) {
      nondominatedAtoms.push_back(indAtomInsert);
    }

  }


}

/*----------------------------------------------------------------------------*/
double Solver::computeRelativeAreaAtoms(const std::vector<uint>& atoms) {


  Score precScore = _instance.getAtomScore(atoms.front());

  double relativeArea = precScore.second*(_instance.nPositives()-precScore.first);

  for(uint ind = 1; ind < atoms.size(); ind ++) {
    Score score = _instance.getAtomScore(atoms[ind]);
    relativeArea += (score.second - precScore.second)*(_instance.nPositives()-score.first);
    precScore = score;
  }

  relativeArea /= static_cast<double>(_instance.nPositives()*_instance.nNegatives());

  return relativeArea;

}

/*----------------------------------------------------------------------------*/
double Solver::relativeArea(/*const Instance& instance, uint bodyLength, */const vector<Score>& supported) {

  // uint maxPos = instance.nPositives()*bodyLength;
  // uint maxNeg = instance.nNegatives()*bodyLength;
  //
  // unsigned long completeArea = maxPos*maxNeg;
  //
  // unsigned long dominatedArea = 0;
  //
  // /* first point */
  // dominatedArea += (maxPos - supported.front().first) * supported.front().second;
  //
  // for(uint ind = 1; ind < supported.size(); ind ++) {
  //   uint width = maxPos - supported[ind].first;
  //   uint height = supported[ind].second - supported[ind-1].second;
  //   dominatedArea += width*height;
  // }
  //
  // return static_cast<double>(dominatedArea) / static_cast<double>(completeArea);

  if(supported.size() <= 2) {
    return 1.;
  }

  unsigned long completeArea = supported.back().first - supported.front().first;
  completeArea *= (supported.back().second - supported.front().second);

  unsigned long dominatedArea = 0;

  for(uint ind = 1; ind < supported.size()-1; ind ++) {
    uint width = supported.back().first - supported[ind].first;
    uint height = supported[ind].second - supported[ind-1].second;
    dominatedArea += width*height;
  }

  return static_cast<double>(dominatedArea) / static_cast<double>(completeArea);

}


/*----------------------------------------------------------------------------*/
void Solver::maxMinDist(const vector<Score>& supported, double& maxDist, double& relativePos) {

  if(supported.size() <= 2) {
    maxDist = 1.0;
    relativePos = 0.5;
    return;
  }

  const Score& p1 = supported.front();
  const Score& p2 = supported.back();

  // cout << "left: " << p1.first << " ; " << p1.second << endl;
  // cout << "right: " << p2.first << " ; " << p2.second << endl;

  /* equation a*x + b*y + c = 0 */
  double a, b, c;
  if(p2.second-p1.second > p2.first-p1.first) {
    b = 1.;
    a = -static_cast<double>(p2.second-p1.second)/(p2.first-p1.first);
  } else {
    b = -static_cast<double>(p2.first-p1.first)/(p2.second-p1.second);
    a = 1.;
  }
  c = -a*p1.first - b*p1.second;

  double k = sqrt(a*a+b*b);

  maxDist = 0.;
  Score projectedPoint;

  for(uint ind = 1; ind < supported.size()-1; ind ++) {

    const Score& point = supported[ind];

    /* compute the projected point */
    Score p3;
    p3.first = static_cast<uint>((b*(b*point.first - a*point.second) -a*c)/(a*a+b*b));
    p3.second = static_cast<uint>((a*(-b*point.first + a*point.second) -b*c)/(a*a+b*b));

    double dist = abs(a*point.first + b*point.second+c)/k;

    if(ind == 1 || dist > maxDist) {
      maxDist = dist;
      projectedPoint = p3;
      // cout << "better point: " << point.first << " ; " << point.second << endl;
    }

  }

  // cout << "best maxDist: " << maxDist << endl;

  /* compute a relative distance by using the distance to utopia point */
  Score utopia(supported.front().first, supported.back().second);
  double distUtopia = abs(a*utopia.first + b*utopia.second+c)/k;
  maxDist /= distUtopia;

  relativePos = sqrt( (projectedPoint.first-p1.first)*(projectedPoint.first-p1.first) + (projectedPoint.second-p1.second)*(projectedPoint.second-p1.second) );
  relativePos = relativePos /  sqrt( (p2.first-p1.first)*(p2.first-p1.first) + (p2.second-p1.second)*(p2.second-p1.second) );

}

/*----------------------------------------------------------------------------*/
void Solver::computeDistances(const vector<Score>& supported, vector<double>& distances) {

  distances.resize(supported.size(), 0);

  if(supported.size() <= 2) {
    distances.push_back(0);
    return;
  }

  const Score& p1 = supported.front();
  const Score& p2 = supported.back();

  distances.front() = 0;
  distances.back() = 0;

  // cout << "left: " << p1.first << " ; " << p1.second << endl;
  // cout << "right: " << p2.first << " ; " << p2.second << endl;

  /* equation a*x + b*y + c = 0 */
  double a, b, c;
  if(p2.second-p1.second > p2.first-p1.first) {
    b = 1.;
    a = -static_cast<double>(p2.second-p1.second)/(p2.first-p1.first);
  } else {
    b = -static_cast<double>(p2.first-p1.first)/(p2.second-p1.second);
    a = 1.;
  }
  c = -a*p1.first - b*p1.second;

  double k = sqrt(a*a+b*b);

  for(uint ind = 1; ind < supported.size()-1; ind ++) {

    const Score& point = supported[ind];

    /* compute the projected point */
    Score p3;
    p3.first = static_cast<uint>((b*(b*point.first - a*point.second) -a*c)/(a*a+b*b));
    p3.second = static_cast<uint>((a*(-b*point.first + a*point.second) -b*c)/(a*a+b*b));

    distances[ind] = abs(a*point.first + b*point.second+c)/k;
  }

}

/*----------------------------------------------------------------------------*/
void Solver::medianMean(const vector<Score>& supported, const Instance& instance, double& medPosMean, double& medNegMean) {

  uint ind = static_cast<uint>(floor(static_cast<double>(supported.size())/2.));

  if(supported.size() % 2 == 0) {

    medPosMean = static_cast<double>(supported[ind-1].first);
    medPosMean += static_cast<double>(supported[ind].first);
    medPosMean /= (2.*instance.nPositives());

    medNegMean = static_cast<double>(supported[ind-1].second);
    medNegMean += static_cast<double>(supported[ind].second);
    medNegMean /= 2.*instance.nNegatives();


  } else {
    medPosMean = static_cast<double>(supported[ind].first)/instance.nPositives();
    medNegMean = static_cast<double>(supported[ind].second)/instance.nNegatives();
  }

}

/*----------------------------------------------------------------------------*/
void Solver::meanVar(const vector<Score>& supported, const vector<Body>& bodies, const Instance& instance, pair<double, double>& posVarMean, pair<double, double>& negVarMean) {

  posVarMean = pair<double,double>(0.,0.);
  negVarMean = pair<double,double>(0.,0.);

  for(uint ind = 0; ind < bodies.size(); ind ++) {

    Histogram histo(instance, bodies[ind], Histogram::Type::POSITIVE_NEGATIVE);

    double posMean = static_cast<double>(supported[ind].first)/instance.nPositives();
    double negMean = static_cast<double>(supported[ind].second)/instance.nNegatives();
    double posVar = 0, negVar = 0;
    for(uint score = 0; score <= bodies.front().size(); score ++) {
      posVar += (score-posMean)*(score-posMean)*static_cast<double>(histo.getPosSampleScore(score).size());
      negVar += (score-negMean)*(score-negMean)*static_cast<double>(histo.getNegSampleScore(score).size());
    }
    posVar /= instance.nPositives();
    negVar /= instance.nNegatives();

    posVarMean.first += posMean;
    posVarMean.second += posVar;

    negVarMean.first += negMean;
    negVarMean.second += negVar;

  }

  posVarMean.first /= static_cast<double>(bodies.size());
  posVarMean.second /= static_cast<double>(bodies.size());

  negVarMean.first /= static_cast<double>(bodies.size());
  negVarMean.second /= static_cast<double>(bodies.size());
}

/*----------------------------------------------------------------------------*/
void Solver::meanVarMed(const vector<Score>& supported, const vector<Body>& bodies, const Instance& instance, pair<double, double>& posMed, pair<double, double>& negMed) {

  vector<double> posMeans(supported.size()), posVars(supported.size());
  vector<double> negMeans(supported.size()), negVars(supported.size());

  for(uint ind = 0; ind < bodies.size(); ind ++) {

    Histogram histo(instance, bodies[ind], Histogram::Type::POSITIVE_NEGATIVE);

    double posMean = static_cast<double>(supported[ind].first)/instance.nPositives();
    double negMean = static_cast<double>(supported[ind].second)/instance.nNegatives();
    double posVar = 0, negVar = 0;
    for(uint score = 0; score <= bodies.front().size(); score ++) {
      posVar += (score-posMean)*(score-posMean)*static_cast<double>(histo.getPosSampleScore(score).size());
      negVar += (score-negMean)*(score-negMean)*static_cast<double>(histo.getNegSampleScore(score).size());
    }
    posVar /= instance.nPositives();
    negVar /= instance.nNegatives();

    posMeans[ind] = posMean;
    posVars[ind] = posVar;

    negMeans[ind] = negMean;
    negVars[ind] = negVar;

  }

  sort(posMeans.begin(), posMeans.end());
  sort(posVars.begin(), posVars.end());
  sort(negMeans.begin(), negMeans.end());
  sort(negVars.begin(), negVars.end());

  uint ind = static_cast<int>(supported.size()/2);

  if(supported.size()%2 == 0) {
    posMed.first = (posMeans[ind-1]+posMeans[ind])/2.;
    posMed.second = (posVars[ind-1]+posVars[ind])/2.;
    negMed.first = (negMeans[ind-1]+negMeans[ind])/2.;
    negMed.second = (negVars[ind-1]+negVars[ind])/2.;

  } else {
    posMed.first = posMeans[ind];
    posMed.second = posVars[ind];
    negMed.first = negMeans[ind];
    negMed.second = negVars[ind];
  }

}


/*----------------------------------------------------------------------------*/
void Solver::tuneInstance(const Instance& instance, uint bodyLength, double ratio, vector<bool>& classVector,  vector<bool>& classVectorModified) {

  vector<Body> solutions;
  vector<Score> points;

  Solver solver(instance);
  solver.computeSupported(bodyLength, &points, nullptr, &solutions);

  // cout << "Relative area (tuning) : " << Solver::relativeArea(points) << endl;

  vector<uint> posScores(instance.nPositives(), 0), negScores(instance.nNegatives(), 0);

  /* compute the score for each sample: sum of score of each rule */
  for(uint indSol = 0; indSol < solutions.size(); indSol ++) {
    Histogram histogram(instance, solutions[indSol], Histogram::Type::POSITIVE_NEGATIVE);
    for(uint score = 0; score <= bodyLength; score ++) {
      auto& posSamplesScore = histogram.getPosSampleScore(score);
      for(auto& ind : posSamplesScore) {
        posScores[ind] += score;
      }
      auto& negSamplesScore = histogram.getNegSampleScore(score);
      for(auto& ind : negSamplesScore) {
        negScores[ind] += score;
      }
    }
  }

  vector<uint> posSamples, negSamples;
  Combinatorics::generateRange(instance.nPositives(), posSamples);
  Combinatorics::generateRange(instance.nNegatives(), negSamples);

  sort(posSamples.begin(), posSamples.end(), [&posScores](const uint& left, const uint& right) { return posScores[left] > posScores[right]; } );
  sort(negSamples.begin(), negSamples.end(), [&negScores](const uint& left, const uint& right) { return negScores[left] < negScores[right]; } );

  /*cout << "pos score: " << endl;
  for(uint ind = 0; ind < 5; ind ++) {
    cout << posScores[posSamples[ind]] << " ";
  }
  cout << endl;

  cout << "neg score: " << endl;
  for(uint ind = 0; ind < 5; ind ++) {
    cout << negScores[negSamples[ind]] << " ";
  }
  cout << endl;*/

  classVectorModified = classVector;
  uint posLimit = static_cast<uint>(static_cast<double>(posSamples.size())*ratio);
  if(posLimit < instance.nPositives()-1) {
    for(uint ind = 0; ind < posLimit; ind ++) {
      classVectorModified[instance.getPosClassIndex(posSamples[ind])] = false;
    }
  }

  uint negLimit = static_cast<uint>(static_cast<double>(negSamples.size())*ratio);
  if(negLimit < instance.nNegatives()-1) {
    for(uint ind = 0; ind < negLimit; ind ++) {
      classVectorModified[instance.getNegClassIndex(negSamples[ind])] = true;
    }
  }
}



/*----------------------------------------------------------------------------*/
ostream& operator<<(ostream& str, const pair<int,int>& tuple) {
  str << "(" << tuple.first << "," << tuple.second << ")";
  return str;
}
