/*!
 * \file NetworkInduction.cpp
 * \author S. Buchet
 * \brief implementation of class NetworkInduction
 */

#include "NetworkInduction.hpp"

#include "Constants.hpp"

#include <chrono>
#include <fstream>
#include <algorithm>

#include "Solver.hpp"


using namespace std;


/*----------------------------------------------------------------------------*/
void NetworkInduction::computeNetwork(DataFrame<uint>& dataset, double selectionThreshold, string output_filename, vector<bool>* negativeCells) {

  // cout << dataset.toString() << endl;

  uint nVar = dataset.nColumns();

  #if DEBUG_REDUCTION == 1
  nVar = NVAR_DEBUG;
  #endif

  cout << "hello network :)" << endl;
  cout << "nvar to do: " << nVar << endl;

  /* list of atoms to test */
  vector<pair<uint,uint>> target;

  for(uint ind = 0; ind < nVar; ind ++) {
    /* adaptive discretization */
    uint nVal = dataset.nUnique(ind);
    for(uint val = 0; val < nVal; val ++) {
      target.push_back(pair<uint,uint>(ind, val));
    }
  }

  /* output for the file */
  vector<string> output(target.size());

  /* computation time */
  // #if DEBUG_LOG == 1
  auto startTime = chrono::steady_clock::now();
  // #endif

  /* nb atoms done */
  #if DEBUG_LOG == 1
  uint nAtomsDone = 0;
  #endif

  int cpt = 0;

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif



  for(uint targetIndex = 0; targetIndex < target.size(); targetIndex ++) {

    auto& elt = target[targetIndex];

    /* index of the feature to exclude */
    #if DEBUG_LOG == 1
    string debugStr;
    debugStr += "target feature: " + dataset.getColLabel(elt.first) + "\n";
    debugStr += "target feature index: " + to_string(elt.first) + "\n";
    debugStr += "target value: " + to_string(elt.second) + "\n";
    #endif


    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : debut" << endl;

    /* create the class vector for the classification instance */
    vector<bool> classVector(dataset.nRows(), false);
    vector<uint> features;

    uint nPos =  Instance::initCoexprClass(dataset, elt.first, elt.second, features, classVector);

    /* additional step: if negative cells are given, modify the instance to go away from a subset of the cells */
    if(negativeCells != nullptr) {
      for(uint indCell = 0; indCell < negativeCells->size(); indCell ++) {
        if((*negativeCells)[indCell] && classVector[indCell]) {
          classVector[indCell] = false;
        }
      }
    }

    // debugStr += "nPos: " + to_string(nPos) + "\n";

    uint nNeg = classVector.size()-nPos;

    #if DEBUG_LOG == 1
      debugStr += "nPos: " + to_string(nPos) + "\n";
      debugStr += "nNeg: " + to_string(nNeg) + "\n";
      debugStr += "percentage: " + to_string(min(nPos, nNeg)/static_cast<double>(nPos+nNeg)) + "\n";
    #endif

    double pct = min(nPos, nNeg)/static_cast<double>(nPos+nNeg);

    /* make sure the instance is well conditioned -> non constant feature */
    if(nPos > 0 && nPos < classVector.size() && pct >= 0.03) {

      cpt += 1;

      #if DEBUG_LOG == 1
      debugStr += "* creation of the instance\n";
      auto start = chrono::steady_clock::now();
      #endif

      Instance instance(dataset, features, classVector);

      /*------------------------------------------------------------------------*/
      /* instance tuning */
      // vector<bool> classVectorModified;
      // Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);
      // Instance instanceModified(dataset, features, classVectorModified, 2);


      #if DEBUG_LOG == 1
      auto stop = chrono::steady_clock::now();
      debugStr += "instance creation execution time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      debugStr += instance.toString();
      debugStr += "* solving the instance\n";
      #endif


      // Solver solver(instanceModified);
      Solver solver(instance);

      /*--------------------------------------------------------------------------*/
      /* computation of the best atoms for thats instance */

      vector<uint> atomsNetwork;
      vector<double> atomsNetworkScores;

      #if DEBUG_LOG == 1
      debugStr += "Computation of the best atoms:\n";
      start = chrono::steady_clock::now();
      #endif

      // solver.computekBestAtoms(nGenes, atomsNetwork, atomsNetworkScores);
      // selectionThreshold = 0.3;
      solver.computeBestAtomsThreshold(selectionThreshold, atomsNetwork, atomsNetworkScores);

      #if DEBUG_LOG == 1
      debugStr += "n atom selected: " + to_string(atomsNetwork.size()) + "\n";
      #endif

      #if DEBUG_LOG == 1
      stop = chrono::steady_clock::now();
      debugStr += "best atoms computation time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      #endif

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
      for(uint ind = 0; ind < atomsNetwork.size(); ind ++) {
        Atom atom = instance.getAtom(atomsNetwork[ind]);
        // cout << dataset.getColLabel(atom.first) << "_" << atom.second << ": " << atomsNetworkScores[ind] << endl;
        output[targetIndex] += " " + dataset.getColLabel(atom.first) + " " + to_string(atom.second);
        output[targetIndex] += " " + to_string(atomsNetworkScores[ind]) + " " + to_string(instance.getAtomPosScore(atomsNetwork[ind])) + " " + to_string(instance.getAtomNegScore(atomsNetwork[ind]));
      }

    } else { /* otherwise register a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(nPos) + " " + to_string(classVector.size()-nPos);


    }

    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
    #endif

    #if DEBUG_LOG == 1
    nAtomsDone ++;
    cout << "n atoms done yet: " << nAtomsDone << endl;
    #endif

    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : fin" << endl;

  }

  #if USE_OPENMP == 1
  }
  #endif

  /* output file containing the area */
  ofstream file(output_filename);
  for(auto& line: output) {
    file << line << endl;
  }
  file.close();

  // #if DEBUG_LOG == 1
  auto stopTime = chrono::steady_clock::now();
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(stopTime-startTime).count() << endl;
  // #endif

  cout << "cpt= " << cpt << endl;

}


/*----------------------------------------------------------------------------*/
void NetworkInduction::computeRegulationNetwork(DataFrame<uint>& dataset, NGraph successors, double trRate, uint predNeq, double selectionThreshold, string output_filename) {

  uint nVar = dataset.nColumns();

  #if DEBUG_REDUCTION == 1
  nVar = NVAR_DEBUG;
  #endif

  /* list of atoms to test */
  vector<pair<uint,uint>> target;

  for(uint ind = 0; ind < nVar; ind ++) {
    /* adaptive discretization */
    uint nVal = dataset.nUnique(ind);
    for(uint val = 0; val < nVal; val ++) {
      target.push_back(pair<uint,uint>(ind, val));
    }
  }


  /* output for the file */
  vector<string> output(target.size());

  /* computation time */
  // #if DEBUG_LOG == 1
  auto startTime = chrono::steady_clock::now();
  // #endif

  /* nb atoms done */
  #if DEBUG_LOG == 1
  uint nAtomsDone = 0;
  #endif

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif

  for(uint targetIndex = 0; targetIndex < target.size(); targetIndex ++) {

    auto& elt = target[targetIndex];

    /* index of the feature to exclude */
    #if DEBUG_LOG == 1
    string debugStr;
    debugStr += "target feature: " + dataset.getColLabel(elt.first) + "\n";
    debugStr += "target feature index: " + to_string(elt.first) + "\n";
    debugStr += "target value: " + to_string(elt.second) + "\n";
    #endif


    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : debut" << endl;

    /* create the class vector for the classification instance */
    vector<bool> classVector(dataset.nRows(), false);

    uint nPos = Instance::initRegulationClass(dataset, successors, elt.first, elt.second, trRate, predNeq, classVector);

    // debugStr += "nPos: " + to_string(nPos) + "\n";

    /* make sure the instance is well conditioned -> non constant feature */
    if(nPos > 0 && nPos < classVector.size()) {

      #if DEBUG_LOG == 1
      debugStr += "* creation of the instance\n";
      auto start = chrono::steady_clock::now();
      #endif

      Instance instance(dataset, classVector);

      /*------------------------------------------------------------------------*/
      /* instance tuning */
      // vector<bool> classVectorModified;
      // Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);
      // Instance instanceModified(dataset, features, classVectorModified, 2);


      #if DEBUG_LOG == 1
      auto stop = chrono::steady_clock::now();
      debugStr += "instance creation execution time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      debugStr += instance.toString();
      debugStr += "* solving the instance\n";
      #endif


      // Solver solver(instanceModified);
      Solver solver(instance);

      /*--------------------------------------------------------------------------*/
      /* computation of the best atoms for thats instance */

      vector<uint> atomsNetwork;
      vector<double> atomsNetworkScores;

      #if DEBUG_LOG == 1
      debugStr += "Computation of the best atoms:\n";
      start = chrono::steady_clock::now();
      #endif

      solver.computeBestAtomsThreshold(selectionThreshold, atomsNetwork, atomsNetworkScores);

      #if DEBUG_LOG == 1
      debugStr += "n atoms selected: " + to_string(atomsNetwork.size()) + "\n";
      #endif

      #if DEBUG_LOG == 1
      stop = chrono::steady_clock::now();
      debugStr += "best atoms computation time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      #endif

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
      for(uint ind = 0; ind < atomsNetwork.size(); ind ++) {
        Atom atom = instance.getAtom(atomsNetwork[ind]);
        // cout << dataset.getColLabel(atom.first) << "_" << atom.second << ": " << atomsNetworkScores[ind] << endl;
        output[targetIndex] += " " + dataset.getColLabel(atom.first) + " " + to_string(atom.second);
        output[targetIndex] += " " + to_string(atomsNetworkScores[ind]) + " " + to_string(instance.getAtomPosScore(atomsNetwork[ind])) + " " + to_string(instance.getAtomNegScore(atomsNetwork[ind]));
      }

    } else { /* otherwise record a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
    }


    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
    #endif

    #if DEBUG_LOG == 1
    nAtomsDone ++;
    cout << "n atoms done yet: " << nAtomsDone << endl;
    #endif

    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : fin" << endl;

  }

  #if USE_OPENMP == 1
  }
  #endif

  /* output file containing the area */
  ofstream file(output_filename);
  for(auto& line: output) {
    file << line << endl;
  }
  file.close();

  // #if DEBUG_LOG == 1
  auto stopTime = chrono::steady_clock::now();
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(stopTime-startTime).count() << endl;
  // #endif

}


/*----------------------------------------------------------------------------*/
void NetworkInduction::computeNetworkOld(DataFrame<uint>& dataset, uint bodyLength, uint nGenes, double rate, string output_filename) {

  // cout << dataset.toString() << endl;

  /* list of atoms to test */
  vector<pair<uint,uint>> target;

  uint nVar = dataset.nColumns();

  #if DEBUG_REDUCTION == 1
  nVar = NVAR_DEBUG;
  #endif

  for(uint ind = 0; ind < nVar; ind ++) {

    /* adaptive discretization */
    uint nVal = dataset.nUnique(ind);

    for(uint val = 0; val < nVal; val ++) {
      target.push_back(pair<uint,uint>(ind, val));
    }

    // target.push_back(pair<uint,uint>(ind, 0));
    // target.push_back(pair<uint,uint>(ind, 1));
    // target.push_back(pair<uint,uint>(ind, 2));
  }

  /* output for the file */
  vector<string> output(target.size());

  /* computation time */
  #if DEBUG_LOG == 1
  auto startTime = chrono::steady_clock::now();
  #endif

  /* nb atoms done */
  #if DEBUG_LOG == 1
  uint nAtomsDone = 0;
  #endif

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif

  for(uint targetIndex = 0; targetIndex < target.size(); targetIndex ++) {


    auto& elt = target[targetIndex];

    /* index of the feature to exclude */
    #if DEBUG_LOG == 1
    string debugStr;
    debugStr += "target feature: " + dataset.getColLabel(elt.first) + "\n";
    debugStr += "target feature index: " + to_string(elt.first) + "\n";
    debugStr += "target value: " + to_string(elt.second) + "\n";
    #endif


    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : debut" << endl;

    /* create the class vector for the classification instance */
    vector<bool> classVector(dataset.nRows(), false);
    vector<uint> features;

    uint nPos =  Instance::initCoexprClass(dataset, elt.first, elt.second, features, classVector);

    // debugStr += "nPos: " + to_string(nPos) + "\n";

    /* make sure the instance is well conditioned -> non constant feature */
    if(nPos > 0 && nPos < classVector.size()) {

      #if DEBUG_LOG == 1
      debugStr += "* creation of the instance\n";
      auto start = chrono::steady_clock::now();
      #endif

      Instance instance(dataset, features, classVector);


      /*------------------------------------------------------------------------*/
      /* instance tuning */
      // vector<bool> classVectorModified;
      // Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);
      // Instance instanceModified(dataset, features, classVectorModified, 2);


      #if DEBUG_LOG == 1
      auto stop = chrono::steady_clock::now();
      debugStr += "instance creation execution time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      debugStr += instance.toString();
      debugStr += "* solving the instance\n";
      #endif


      // Solver solver(instanceModified);
      Solver solver(instance);

      /*--------------------------------------------------------------------------*/
      /* computation of all supported points */

      vector<Body> solutions;
      vector<Score> points;
      vector<Weight> weights;

      #if DEBUG_LOG == 1
      debugStr += "Computation of the supported solutions:\n";
      start = chrono::steady_clock::now();
      #endif

      solver.computeSupported(bodyLength, &points, &weights, &solutions);

      #if DEBUG_LOG == 1
      debugStr += "N supported points: " + to_string(solutions.size()) + "\n";
      debugStr += "first point : " + to_string(points.front().first) + ", " + to_string(points.front().second) + "\n";
      debugStr += "last point : " + to_string(points.back().first) + ", " + to_string(points.back().second) + "\n";
      #endif

      #if DEBUG_LOG == 1
      stop = chrono::steady_clock::now();
      debugStr += "computation supported time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      #endif

      double area = Solver::relativeArea(points);

      /* compute the distances of all points to the diagonal */
      vector<double> distances;
      Solver::computeDistances(points, distances);
      /* sort: farthest first */
      vector<uint> selectedSol;
      Combinatorics::generateRange(static_cast<uint>(points.size()), selectedSol);
      sort(selectedSol.begin(), selectedSol.end(), [&distances](const uint& left, const uint& right) { return distances[left] > distances[right]; });
      uint limitIndex = static_cast<int>(rate*static_cast<double>(selectedSol.size()));
      selectedSol.resize(limitIndex);

      /* create a list of the atoms considered, scored by their solution ranked */
      map<uint, uint> selectedAtoms;
      for(auto& indSol: selectedSol) {
        for(uint atomRank = 0; atomRank < solutions[indSol].size(); atomRank ++) {
          uint atomInd = solutions[indSol][atomRank];
          if(selectedAtoms.count(atomInd) == 0) {
            selectedAtoms.insert(pair<uint,uint>(atomInd, bodyLength-atomRank));
          } else {
            selectedAtoms[atomInd] += bodyLength - atomRank;
          }
        }
      }

      /* create the final list of genes */
      vector<uint> atomsNetwork(selectedAtoms.size());
      {
        uint ind = 0;
        for(auto it = selectedAtoms.begin(); it != selectedAtoms.end(); it ++) {
          atomsNetwork[ind] = it->first;
          ind += 1;
        }
      }

      sort(atomsNetwork.begin(), atomsNetwork.end(), [&selectedAtoms](const uint& left, const uint& right) {return selectedAtoms[left] > selectedAtoms[right]; });
      atomsNetwork.resize(nGenes);

      #if DEBUG_LOG == 1
      debugStr += "Relative area: " + to_string(area) + "\n";
      #endif

      // file << dataset.getColLabel(elt.first) << " " << elt.second << " " << area << " " << instance.nPositives()  << " " << instance.nNegatives() << endl;


      /* update the ouptut */
      // output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(maxDist) + " " + to_string(relativeDist) + " ";
      // output[targetIndex] += to_string(medPosMean) + " " + to_string(medNegMean) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
      for(uint ind = 0; ind < atomsNetwork.size(); ind ++) {
        Atom atom = instance.getAtom(atomsNetwork[ind]);
        output[targetIndex] += " " + dataset.getColLabel(atom.first) + " " + to_string(atom.second);
        output[targetIndex] += " " + to_string(selectedAtoms[atomsNetwork[ind]]) + " " + to_string(instance.getAtomPosScore(atomsNetwork[ind])) + " " + to_string(instance.getAtomNegScore(atomsNetwork[ind]));
      }

    } else { /* otherwise register a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);

    }

    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
    #endif

    #if DEBUG_LOG == 1
    nAtomsDone ++;
    cout << "n atoms done yet: " << nAtomsDone << endl;
    #endif

    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : fin" << endl;

  }

  #if USE_OPENMP == 1
  }
  #endif

  /* output file containing the area */
  ofstream file(output_filename);
  for(auto& line: output) {
    file << line << endl;
  }
  file.close();

  #if DEBUG_LOG == 1
  auto stopTime = chrono::steady_clock::now();
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(stopTime-startTime).count() << endl;
  #endif

}
