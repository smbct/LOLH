/*!
 * \file ClassificationQuality.cpp
 * \author S. Buchet
 * \brief implementation of class ClassificationQuality
 */

#include "ClassificationQuality.hpp"

#include <fstream>
#include <iostream>
#include <cstdio>

/* time */
#include <chrono>
#include <unistd.h>

#include <set>


#include "Constants.hpp"
#include "Solver.hpp"
#include "Utils.hpp"

using namespace std;



/*----------------------------------------------------------------------------*/
void ClassificationQuality::computeCoexprQuality(DataFrame<uint>& dataset, uint bodyLength, string output_filename) {

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

  #if DEBUG_LOG == 1
  uint nAtomsDone = 0;
  #endif

  /* output for the file */
  vector<string> output(target.size());

  /* computation time */
  // #if DEBUG_LOG == 1
  auto startTime = chrono::steady_clock::now();
  // #endif

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

      // Instance instance(dataset, features, classVector, 2);
      Instance instance(dataset, features, classVector); // adaptive discretization


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
      stop = chrono::steady_clock::now();
      debugStr += "computation supported time: " + to_string(chrono::duration_cast<chrono::milliseconds>(stop-start).count()) + "\n";
      #endif

      double area = Solver::relativeArea(points);

      /* print the points */
      // for(uint ind = 0; ind < solutions.size(); ind ++) {
      //   cout << points[ind] << endl;
      // }


      /* extra infos: max min dist, median of pos and neg means*/
      double maxDist, relativeDist;
      double medPosMean, medNegMean;
      Solver::maxMinDist(points, maxDist, relativeDist);

      Solver::medianMean(points, instance, medPosMean, medNegMean);


      #if DEBUG_LOG == 1
      debugStr += "Relative area: " + to_string(area) + "\n";
      debugStr += "Max min dist: " + to_string(maxDist) + "\n";
      debugStr += "Relative pos max dist: " + to_string(relativeDist) + "\n";
      debugStr += "Median pos mean: " + to_string(medPosMean) + "\n";
      debugStr += "Median neg mean: " + to_string(medNegMean) + "\n";
      #endif

      #if DEBUG_LOG == 1
      nAtomsDone ++;
      cout << nAtomsDone << " atoms done over " << target.size() << endl;
      #endif

      // file << dataset.getColLabel(elt.first) << " " << elt.second << " " << area << " " << instance.nPositives()  << " " << instance.nNegatives() << endl;


      /* update the ouptut */
      // output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());
      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(maxDist) + " " + to_string(relativeDist) + " ";
      output[targetIndex] += to_string(medPosMean) + " " + to_string(medNegMean) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());

    } else { /* otherwise register a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 -1 -1 -1 -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
    }

    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
    #endif

    // cout << dataset.getColLabel(elt.first) << " - " << elt.first << " - " << elt.second << " : fin" << endl;

    // #if DEBUG_LOG == 1
    // nAtomsDone ++;
    // if(nAtomsDone % 500 == 0) {
    //   cout << nAtomsDone << " atoms done over " << target.size() << endl;
    // }
    // #endif


  }

  #if USE_OPENMP == 1
  }
  #endif

  auto stopTime = chrono::steady_clock::now();

  /* output file containing the area */
  ofstream file(output_filename);
  for(auto& line: output) {
    file << line << endl;
  }
  // file << "time: " << to_string(chrono::duration_cast<chrono::milliseconds>(stopTime-startTime).count()) << endl;
  file.close();

  #if DEBUG_LOG == 1
  cout << "computation time: " << chrono::duration_cast<chrono::milliseconds>(stopTime-startTime).count() << endl;
  #endif
}

/*----------------------------------------------------------------------------*/
void ClassificationQuality::computeFastCoexprQuality(DataFrame<uint>& dataset, uint bodyLength, string output_filename) {

  // cout << dataset.toString() << endl;

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

    #if DEBUG_LOG == 1
    debugStr += "nPos: " + to_string(nPos) + "\n";
    #endif

    /* make sure the instance is well conditions -> non constant feature */
    if(nPos > 0 && nPos < classVector.size()) {

      #if DEBUG_LOG == 1
      debugStr += "* creation of the instance\n";
      auto start = chrono::steady_clock::now();
      #endif

      Instance instance(dataset, features, classVector, 2);

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
      /* computation the extreme points and the farthest */
      Score rPlusScore;
      Body rPlusBody(bodyLength);

      rPlusScore = solver.solveMono(bodyLength, Objective::POSITIVE, &rPlusBody);

      Score rMinusScore;
      Body rMinusBody;

      rMinusScore = solver.solveMono(bodyLength, Objective::NEGATIVE, &rMinusBody);

      // debugStr += "Left point: " + to_string(rPlusScore.first) + " ; " + to_string(rPlusScore.second) + "\n";
      // debugStr += "Right point: " + to_string(rMinusScore.first) + " ; " + to_string(rMinusScore.second) + "\n";


      std::vector<Score> points;
      std::vector<Body> solutions;

      // cout << "weights: " << static_cast<int>(rMinusScore.second)-static_cast<int>(rPlusScore.second) << " ; " << static_cast<int>(rPlusScore.first)-static_cast<int>(rMinusScore.first) << endl;

      solver.solveScalarization(bodyLength, Weight(static_cast<int>(rMinusScore.second)-static_cast<int>(rPlusScore.second), static_cast<int>(rPlusScore.first)-static_cast<int>(rMinusScore.first)), &points, &solutions);

      // debugStr += "N points first scalarization: " + to_string(points.size()) + "\n";
      // debugStr += "Farthest suported point: " + to_string(points.front().first) + " ; " + to_string(points.front().second) + "\n";

      /*debugStr += to_string(points.size()) + "\n";
      for(auto& elt: points) {
        debugStr += to_string(elt.first) + " " + to_string(elt.second) + "\n";
      }*/

      #if DEBUG_LOG == 1
      debugStr += "Computation of the quality:\n";
      start = chrono::steady_clock::now();
      #endif

      /* extra infos: max min dist, median of pos and neg means*/
      double maxDist, relativeDist;

      /* equation a*x + b*y + c = 0 */
      double a, b, c;
      if(rMinusScore.second-rPlusScore.second > rMinusScore.first-rPlusScore.first) {
        b = 1.;
        a = -static_cast<double>(rMinusScore.second-rPlusScore.second)/(rMinusScore.first-rPlusScore.first);
      } else {
        b = -static_cast<double>(rMinusScore.first-rPlusScore.first)/(rMinusScore.second-rPlusScore.second);
        a = 1.;
      }
      c = -a*rPlusScore.first - b*rPlusScore.second;

      double k = sqrt(a*a+b*b);

      maxDist = abs(a*points.front().first + b*points.front().second+c)/k;

      /* distance from utopia */
      Score utopian(rPlusScore.first, rMinusScore.second);
      double distUtopian = abs(a*utopian.first + b*utopian.second+c)/k;

      relativeDist = maxDist/distUtopian;

      /* point projected on the diagonal */
      Score projected;
      projected.first = static_cast<uint>((b*(b*points.front().first - a*points.front().second) -a*c)/(a*a+b*b));
      projected.second = static_cast<uint>((a*(-b*points.front().first + a*points.front().second) -b*c)/(a*a+b*b));

      double relativePos = sqrt( (projected.first-rPlusScore.first)*(projected.first-rPlusScore.first) + (projected.second-rPlusScore.second)*(projected.second-rPlusScore.second) );
      relativePos = relativePos /  sqrt( (rMinusScore.first-rPlusScore.first)*(rMinusScore.first-rPlusScore.first) + (rMinusScore.second-rPlusScore.second)*(rMinusScore.second-rPlusScore.second) );


      #if DEBUG_LOG == 1
      debugStr += "Max min dist: " + to_string(maxDist) + "\n";
      debugStr += "Relative dist: " + to_string(relativeDist) + "\n";
      debugStr += "Relative pos max dist: " + to_string(relativePos) + "\n";
      #endif

      // file << dataset.getColLabel(elt.first) << " " << elt.second << " " << area << " " << instance.nPositives()  << " " << instance.nNegatives() << endl;

      /* update the ouptut */
      // output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());
      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(relativeDist) + " " + to_string(relativePos) + " ";
      output[targetIndex] += to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());

    } else { /* otherwise register a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
    }

    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
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

/*----------------------------------------------------------------------------*/
void ClassificationQuality::computeAtomCoexprQuality(DataFrame<uint>& dataset, string output_filename) {

  // cout << dataset.toString() << endl;

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
  #if DEBUG_LOG == 1
  auto startTime = chrono::steady_clock::now();
  #endif

  uint nvar = 0;

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

    #if DEBUG_LOG == 1
    debugStr += "nPos: " + to_string(nPos) + "\n";
    #endif

    /* make sure the instance is well conditions -> non constant feature */
    if(nPos > 0 && nPos < classVector.size()) {

      #if DEBUG_LOG == 1
      debugStr += "* creation of the instance\n";
      auto start = chrono::steady_clock::now();
      #endif

      Instance instance(dataset, features, classVector, 2);


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

      #if DEBUG_LOG == 1
      debugStr += "Computation of the quality:\n";
      start = chrono::steady_clock::now();
      #endif

      /*--------------------------------------------------------------------------*/
      /* computation the extreme points and the farthest */
      vector<uint> atoms;
      vector<uint> nonDominatedAtoms;
      Combinatorics::generateRange(instance.nAtoms(), atoms);
      solver.computeNonDominatedAtoms(atoms, nonDominatedAtoms);
      double relativeArea = solver.computeRelativeAreaAtoms(nonDominatedAtoms);


      #if DEBUG_LOG == 1
      debugStr += "Relative atom area: " + to_string(relativeArea)  + "\n";
      #endif

      // file << dataset.getColLabel(elt.first) << " " << elt.second << " " << area << " " << instance.nPositives()  << " " << instance.nNegatives() << endl;

      /* update the ouptut */
      // output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());
      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(relativeArea) + " ";
      output[targetIndex] += to_string(instance.nPositives()) + " " + to_string(instance.nNegatives());

    } else { /* otherwise register a special value */

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
    }

    #if DEBUG_LOG == 1
    cout << debugStr << endl << endl;
    #endif

    nvar ++;
    if(nvar % 500 == 0) {
      cout << nvar << " atoms done over " << target.size() << endl;
    }

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


/*----------------------------------------------------------------------------*/
void ClassificationQuality::computeRegulQuality(DataFrame<uint>& dataset, NGraph successors, uint bodyLength, double trRate, uint predNeq, string output_filename) {

  #if DEBUG_LOG == 1
  cout << "n genes: " << dataset.nColumns() << endl;
  #endif

  uint nVar = dataset.nColumns();

  #if DEBUG_REDUCTION == 1
  nVar = NVAR_DEBUG;
  #endif

  /* list of atoms to compute */
  vector<pair<uint,uint>> target;

  // for(uint ind = 0; ind <= dataset.nColumns()/*10*/; ind ++) {
  //   target.push_back(pair<uint,uint>(ind, 0));
  //   target.push_back(pair<uint,uint>(ind, 1));
  // }

  for(uint ind = 0; ind < nVar; ind ++) {
    /* adaptive discretization */
    uint nVal = dataset.nUnique(ind);
    for(uint val = 0; val < nVal; val ++) {
      target.push_back(pair<uint,uint>(ind, val));
    }
  }


  /* output */
  vector<string> output(target.size());


  #if DEBUG_LOG == 1
  uint nN = 0;
  for(auto& elt: successors) {
    nN += static_cast<uint>(elt.size());
  }
  cout << "mean n neighbours: " << static_cast<double>(nN)/static_cast<double>(successors.size()) << endl;
  cout << "n transitions: " << nN << endl;
  #endif

  /* computation time */
  auto start = chrono::steady_clock::now();

  #if DEBUG_LOG == 1
  unsigned int k = 0;
  #endif

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif

  for(unsigned int targetIndex = 0; targetIndex < target.size(); targetIndex ++) {

    auto& elt = target[targetIndex];

    // auto start = chrono::steady_clock::now();

    /* index of the feature to exclude */
    #if DEBUG_LOG == 1
    string debug_str;
    debug_str += "target feature: " + dataset.getColLabel(elt.first) + "\n";
    debug_str += "target feature index: " + to_string(elt.first) + "\n";
    debug_str += "target value: " + to_string(elt.second) + "\n";
    #endif


    /*------------------------------------------------------------------------*/
    /* create the classification instance based on the transition list */

    /* create the class vector */
    vector<bool> classVector(dataset.nRows(), false);
    // auto start2 = chrono::steady_clock::now();

    uint nPos = Instance::initRegulationClass(dataset, successors, elt.first, elt.second, trRate, predNeq, classVector);
    // auto stop2 = chrono::steady_clock::now();
    // cout << "init instance time: " << chrono::duration_cast<chrono::milliseconds>(stop2-start2).count() << endl;

    #if DEBUG_LOG == 1
    debug_str += "Instance positives/negatives: " + to_string(nPos) + " " + to_string(classVector.size()-nPos) + "\n";
    #endif

    /* output class vector */
    /*if(elt.second > 0) {
      ofstream file("classes_" + dataset.getColLabel(elt.first) + "_" + to_string(elt.second) + ".txt");
      file << nPos << endl;
      for(uint ind = 0; ind < classVector.size(); ind ++) {
        if(classVector[ind]) {
          file << dataset.getRowLabel(ind) << endl;
        }
      }
      for(uint ind = 0; ind < classVector.size(); ind ++) {
        if(!classVector[ind]) {
          file << dataset.getRowLabel(ind) << endl;
        }
      }
      file.close();
    }*/


    /*------------------------------------------------------------------------*/
    /* compute cluster frequencies */
    // map<string, int> labelFreq;
    // for(auto& label : labels) {
    //   labelFreq.insert(pair<string, int>(label, 0));
    // }
    // for(unsigned int ind = 0; ind < classVector.size(); ind ++) {
    //   if(classVector[ind]) {
    //     labelFreq[ cellLabel[dataset.getRowLabel(ind)] ] += 1;
    //   }
    // }

    /* make sure the instance is not ill conditioned */
    if(nPos > 0 && nPos < classVector.size()) {

      // Instance instance(dataset, classVector, 2);
      Instance instance(dataset, classVector); /* previously mistake here: not pre-defined n values */

      /*------------------------------------------------------------------------*/
      /* instance tuning */
      // vector<bool> classVectorModified;
      // Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);
      // Instance instanceModified(dataset, classVectorModified, 2);
      /*------------------------------------------------------------------------*/


      #if DEBUG_LOG == 1
      debug_str +=  "instance nPos: " + to_string(instance.nPositives()) + "\n";
      debug_str +=  "instance nNeg: " + to_string(instance.nNegatives()) + "\n";
      #endif

      Solver solver(instance);

      vector<Body> solutions;
      vector<Score> points;
      vector<Weight> weights;

      solver.computeSupported(bodyLength, &points, &weights, &solutions);

      double area = Solver::relativeArea(points);

      #if DEBUG_LOG == 1
      debug_str += "Relative area: " + to_string(area) + "\n";
      #endif

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(area) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives()) /*+ " #"*/;
      // for(auto& label: labels) {
      //   output[targetIndex] += label + ":" + to_string(labelFreq[label]) + ";";
      // }


    } else {

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
      // for(auto& label: labels) {
      //   output[targetIndex] += ":" + to_string(labelFreq[label]) + ";";
      // }

    }


    // auto stop = chrono::steady_clock::now();
    // cout << "gene computation time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

    #if DEBUG_LOG == 1
    debug_str += "\n";
    cout << debug_str << endl;
    #endif

    #if DEBUG_LOG == 1
    k += 1;
    cout << endl << endl << k << " atoms over " << target.size() << endl << endl << endl;
    #endif

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


  auto stop = chrono::steady_clock::now();
  cout << "total computation time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

}

/*----------------------------------------------------------------------------*/
void ClassificationQuality::computeAtomRegulQuality(DataFrame<uint>& dataset, NGraph successors, double trRate, uint predNeq, string output_filename) {

  #if DEBUG_LOG == 1
  cout << "n genes: " << dataset.nColumns() << endl;
  #endif

  uint nVar = dataset.nColumns();

  #if DEBUG_REDUCTION == 1
  nVar = NVAR_DEBUG;
  #endif

  /* list of atoms to compute */
  vector<pair<uint,uint>> target;

  // for(uint ind = 0; ind <= dataset.nColumns()/*10*/; ind ++) {
  //   target.push_back(pair<uint,uint>(ind, 0));
  //   target.push_back(pair<uint,uint>(ind, 1));
  // }

  for(uint ind = 0; ind < nVar; ind ++) {
    /* adaptive discretization */
    uint nVal = dataset.nUnique(ind);
    for(uint val = 0; val < nVal; val ++) {
      target.push_back(pair<uint,uint>(ind, val));
    }
  }


  /* output */
  vector<string> output(target.size());


  #if DEBUG_LOG == 1
  uint nN = 0;
  for(auto& elt: successors) {
    nN += static_cast<uint>(elt.size());
  }
  cout << "mean n neighbours: " << static_cast<double>(nN)/static_cast<double>(successors.size()) << endl;
  cout << "n transitions: " << nN << endl;
  #endif

  /* computation time */
  auto start = chrono::steady_clock::now();

  // #if DEBUG_LOG == 1
  unsigned int k = 0;
  // #endif

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif

  for(unsigned int targetIndex = 0; targetIndex < target.size(); targetIndex ++) {

    auto& elt = target[targetIndex];

    // auto start = chrono::steady_clock::now();

    /* index of the feature to exclude */
    #if DEBUG_LOG == 1
    string debug_str;
    debug_str += "target feature: " + dataset.getColLabel(elt.first) + "\n";
    debug_str += "target feature index: " + to_string(elt.first) + "\n";
    debug_str += "target value: " + to_string(elt.second) + "\n";
    #endif


    /*------------------------------------------------------------------------*/
    /* create the classification instance based on the transition list */

    /* create the class vector */
    vector<bool> classVector(dataset.nRows(), false);
    // auto start2 = chrono::steady_clock::now();

    uint nPos = Instance::initRegulationClass(dataset, successors, elt.first, elt.second, trRate, predNeq, classVector);
    // auto stop2 = chrono::steady_clock::now();
    // cout << "init instance time: " << chrono::duration_cast<chrono::milliseconds>(stop2-start2).count() << endl;

    #if DEBUG_LOG == 1
    debug_str += "Instance positives/negatives: " + to_string(nPos) + " " + to_string(classVector.size()-nPos) + "\n";
    #endif

    /* output class vector */
    /*if(elt.second > 0) {
      ofstream file("classes_" + dataset.getColLabel(elt.first) + "_" + to_string(elt.second) + ".txt");
      file << nPos << endl;
      for(uint ind = 0; ind < classVector.size(); ind ++) {
        if(classVector[ind]) {
          file << dataset.getRowLabel(ind) << endl;
        }
      }
      for(uint ind = 0; ind < classVector.size(); ind ++) {
        if(!classVector[ind]) {
          file << dataset.getRowLabel(ind) << endl;
        }
      }
      file.close();
    }*/


    /*------------------------------------------------------------------------*/
    /* compute cluster frequencies */
    // map<string, int> labelFreq;
    // for(auto& label : labels) {
    //   labelFreq.insert(pair<string, int>(label, 0));
    // }
    // for(unsigned int ind = 0; ind < classVector.size(); ind ++) {
    //   if(classVector[ind]) {
    //     labelFreq[ cellLabel[dataset.getRowLabel(ind)] ] += 1;
    //   }
    // }

    /* make sure the instance is not ill conditioned */
    if(nPos > 0 && nPos < classVector.size()) {

      // Instance instance(dataset, classVector, 2);
      Instance instance(dataset, classVector); /* previously mistake here: not pre-defined n values */

      /*------------------------------------------------------------------------*/
      /* instance tuning */
      // vector<bool> classVectorModified;
      // Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);
      // Instance instanceModified(dataset, classVectorModified, 2);
      /*------------------------------------------------------------------------*/


      #if DEBUG_LOG == 1
      debug_str +=  "instance nPos: " + to_string(instance.nPositives()) + "\n";
      debug_str +=  "instance nNeg: " + to_string(instance.nNegatives()) + "\n";
      #endif

      Solver solver(instance);

      /*--------------------------------------------------------------------------*/
      /* computation the extreme points and the farthest */
      vector<uint> atoms;
      vector<uint> nonDominatedAtoms;
      Combinatorics::generateRange(instance.nAtoms(), atoms);
      solver.computeNonDominatedAtoms(atoms, nonDominatedAtoms);
      double relativeArea = solver.computeRelativeAreaAtoms(nonDominatedAtoms);

      #if DEBUG_LOG == 1
      debug_str += "Relative area: " + to_string(relativeArea) + "\n";
      #endif

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " " + to_string(relativeArea) + " " + to_string(instance.nPositives()) + " " + to_string(instance.nNegatives()) /*+ " #"*/;
      // for(auto& label: labels) {
      //   output[targetIndex] += label + ":" + to_string(labelFreq[label]) + ";";
      // }


    } else {

      output[targetIndex] = dataset.getColLabel(elt.first) + " " + to_string(elt.second) + " -1 " + to_string(nPos) + " " + to_string(classVector.size()-nPos);
      // for(auto& label: labels) {
      //   output[targetIndex] += ":" + to_string(labelFreq[label]) + ";";
      // }

    }


    // auto stop = chrono::steady_clock::now();
    // cout << "gene computation time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

    #if DEBUG_LOG == 1
    debug_str += "\n";
    cout << debug_str << endl;
    #endif

    // #if DEBUG_LOG == 1
    // k += 1;
    // cout << endl << endl << k << " atoms over " << target.size() << endl << endl << endl;
    // #endif

    k += 1;
    if(k % 500 == 0) {
      cout << "n atoms done yet: " << k << " over " << target.size() << endl;
    }

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


  auto stop = chrono::steady_clock::now();
  cout << "total computation time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

}
