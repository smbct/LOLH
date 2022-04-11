/*!
 * \file Induction.cpp
 * \author S. Buchet
 * \brief implementation of class Induction
 */

#include "Induction.hpp"
#include "Solver.hpp"
#include "Histogram.hpp"
#include "LogicProgram.hpp"
#include "Constants.hpp"

#include <algorithm>
#include <set>

using namespace std;

/*----------------------------------------------------------------------------*/
Induction::Induction() {

}


/*----------------------------------------------------------------------------*/
void Induction::ruleInduction(const Instance& instance, Parameters& param, vector<Body>& bodies, bool clusterInstance) {

  bool stop = false;

  /* the first sub-instance is the original one */
  const Instance* subInstance = &instance;

  uint targetCover = instance.nPositives()/param.nRules;
  if(targetCover < 1) {
    targetCover ++;
  }

  if(targetCover < 50) {
    targetCover = std::min(static_cast<uint>(50), subInstance->nPositives());
  }

  #if DEBUG_LOG == 1
  string debugStr;
  debugStr += "Cover target: " + to_string(targetCover) + "\n";
  #endif

  vector<uint> covered, coveredTemp;

  while(!stop) {

    // cout << covered.size() << " out of " << instance.nPositives() << endl;

    if(covered.empty()) {
      subInstance = &instance;
    } else {
      subInstance = new Instance(instance, covered, 2);
    }

    // cout << "sub instance n pos: " << subInstance->nPositives() << endl;


    Solver solver(*subInstance);

    Score targetPoint; /* target solution in the objective space */
    Weight targetWeight; /* weights used for the scalarized objective */
    Body targetSolution; /* solutions computed */

    bool targetFound;
    if(clusterInstance) {
      targetFound = solver.computeTargetSolutionSlow(param.clusterBodyLength, targetCover, param.matchingThreshold, targetPoint, targetWeight, targetSolution);
    } else {
      targetFound  = solver.computeTargetSolutionSlow(param.bodyLength, targetCover, param.matchingThreshold, targetPoint, targetWeight, targetSolution);
    }

    /* delete the sub instance */
    if(subInstance != &instance) {
      delete subInstance;
      subInstance = nullptr;
    }

    if(targetFound) { /* check if the target was possible */

      /* compute covered samples */
      Histogram histogram(instance, targetSolution, Histogram::Type::POSITIVE_NEGATIVE);
      cout << histogram << endl;

      // cout << histogram << endl;

      uint nFalsePositive = histogram.negativeCovered(param.matchingThreshold);
      double falsePositiveRate = static_cast<double>(nFalsePositive)/instance.nNegatives();

      if(falsePositiveRate <= param.falsePositiveRate/*nFalsePositive == 0*/) {
        bodies.push_back(targetSolution);
      } else {
        #if DEBUG_LOG == 1
          debugStr += "rule is rejected\n";
        #endif
      }

      #if DEBUG_LOG == 1
      debugStr += "n false positives: " + to_string(nFalsePositive);
      debugStr += " -> " + to_string(100.*falsePositiveRate) + "%\n";
      #endif


      histogram.getPositiveCovered(param.matchingThreshold, coveredTemp);

      // cout << histogram << endl;

      int cptCov = 0;
      for(auto& elt : coveredTemp) {
        auto it = find(covered.begin(), covered.end(), elt);
        if(it == covered.end()) {
          covered.push_back(elt);
          cptCov ++;
        }
      }

      // cout << "n covered: " << covered.size() << endl;

      // #if DEBUG_LOG == 1
      // cout << "n covered new: " << cptCov << " (-> " << coveredTemp.size() << endl;
      // cout << "target goal: " << targetCover << endl << endl << endl;
      // #endif

      if(covered.size() >= instance.nPositives() - 1) {
        stop = true;
      } else if(covered.size() + targetCover > instance.nPositives()) {
        targetCover = static_cast<uint>(instance.nPositives()-covered.size());
      }

    } else {
      stop = true;

      #if DEBUG_LOG == 1
      debugStr += "Goal cannot be reached anymore";
      #endif

    }

  }



  #if DEBUG_LOG == 1
  debugStr += "Solutions: \n";
  for(auto& body: bodies) {
    cout << body << endl;
    debugStr += bodyToString(body) + "\n";
    // Histogram histogram(instance, body, Histogram::Type::POSITIVE_NEGATIVE);
    // cout << histogram << endl;
  }
  cout << debugStr << endl;
  #endif

}

/*----------------------------------------------------------------------------*/
void Induction::createLFITProgram(string datasetFileName, string transitionsFileName, string cellLabelFileName, Parameters& param) {

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  DataFrame<string> labels;
  labels.loadFromCSV(cellLabelFileName);
  /* create the list of labels */
  set<string> labels_set;
  for(uint ind = 0; ind < labels.nRows(); ind ++) {
    labels_set.insert(labels.getData(ind, 0));
  }

  DataFrame<string>* transitionsPtr = nullptr;

  if(!param.coExpression) {
    transitionsPtr = new DataFrame<string>();
    transitionsPtr->loadFromCSV(transitionsFileName);
  }

  /* create the logic program */
  LogicProgram program;
  vector<string> col;
  dataset.getColumnLabels(col);
  program.create(col, 2);

  /* list of atoms to test */
  vector<Atom> target;
  for(uint ind = 0; ind < dataset.nColumns()/*10*/; ind ++) {
    target.push_back(pair<uint,uint>(ind, 0));
    target.push_back(pair<uint,uint>(ind, 1));
  }

  /* compute a list of successors (indexes) from hte transitions */
  NGraph* successorsPtr = nullptr;
  if(!param.coExpression) {
    successorsPtr = new NGraph();
    DataFrame<uint>::createNeighbourhoodGraph(dataset, *transitionsPtr, param.graphDelay, *successorsPtr);
  }

  /* list of features if co-expression learning */
  vector<uint>* featuresPtr = nullptr;
  if(param.coExpression) {
    featuresPtr = new vector<uint>();
  }

  for(auto& learnedAtom: target) {

    #if DEBUG_LOG == 1
    cout << "Learning " << dataset.getColLabel(learnedAtom.first) << "_" << learnedAtom.second << endl;
    cout << "var index: " << learnedAtom.first << endl;
    #endif

    /*------------------------------------------------------------------------*/
    /* create the classification instance based on the transition list */

    /* create the class vector */
    vector<bool> classVector(dataset.nRows(), false);
    uint nPos;
    if(!param.coExpression) {
      nPos = Instance::initRegulationClass(dataset, *successorsPtr, learnedAtom.first, learnedAtom.second, param.posRate, false, classVector);
    } else {
      nPos = Instance::initCoexprClass(dataset, learnedAtom.first, learnedAtom.second, *featuresPtr, classVector);
    }

    /* make sure the instance is not ill conditioned */
    if(nPos > 0 && nPos < classVector.size()) {

      Instance *instancePtr = nullptr;
      if(!param.coExpression) {
        instancePtr = new Instance(dataset, classVector, 2);
      } else {
        instancePtr = new Instance(dataset, *featuresPtr, classVector, 2);
      }

      #if DEBUG_LOG == 1
      cout << "instance nPos: " << instancePtr->nPositives() << endl;
      cout << "instance nNeg: " << instancePtr->nNegatives() << endl;
      #endif

      vector<Body> solutions;

      ruleInduction(*instancePtr, param, solutions, false);

      /* add the body to the logic program */
      for(auto& body: solutions) {
        vector<Atom> lpBody;
        for(auto& atomInd: body) {
          lpBody.push_back(instancePtr->getAtom(atomInd));
        }
        program.addRule(learnedAtom, lpBody);
      }

      delete instancePtr;
      instancePtr = nullptr;

    }

  }

  /* create the rules for the clusters */
  for(auto   label: labels_set) {

    if(label == "cluster 7") {
      program.addVariable("cluster_7", 2);
    } else {
      program.addVariable(label, 2);
    }

    cout << "Learning label: " << label << endl;

    /* create the class vector */
    vector<bool> classVector(dataset.nRows(), false);

    uint nPos = Instance::initClusterClass(dataset, labels, label, classVector);

    /* make sure the instance is not ill conditioned */
    if(nPos > 0 && nPos < classVector.size()) {

      Instance *instancePtr = new Instance(dataset, classVector, 2);

      #if DEBUG_LOG == 1
      cout << "instance nPos: " << instancePtr->nPositives() << endl;
      cout << "instance nNeg: " << instancePtr->nNegatives() << endl;
      #endif

      vector<Body> solutions;

      ruleInduction(*instancePtr, param, solutions, true);

      /* add the body to the logic program */
      Atom head(program.nVariables()-1, 1);
      for(auto& body: solutions) {
        vector<Atom> lpBody;

        for(auto& atomInd: body) {
          lpBody.push_back(instancePtr->getAtom(atomInd));
        }
        program.addRule(head, lpBody);

      }

      delete instancePtr;
      instancePtr = nullptr;

    }



  }

  /* export the logic program */
  program.exportToFile("program_MAIT_reduced_clusters_b10_.lp");


  if(param.coExpression) {
    delete featuresPtr;
    featuresPtr = nullptr;
  } else {
    delete successorsPtr;
    successorsPtr = nullptr;
    delete transitionsPtr;
    transitionsPtr = nullptr;
  }

}
