/*!
 * \file LinearGraph.cpp
 * \author S. Buchet
 * \brief implementation of class LinearGraph
 */

#include "LinearGraph.hpp"

#include <cmath>

#include <chrono>

#include <iostream>

#include <fstream>

#include "Constants.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
LinearGraph::LinearGraph(DataFrame<double>& dataset) : _dataset(dataset)
{

}

/*----------------------------------------------------------------------------*/
void LinearGraph::computeGraph(double threshold, string fileName) {

  vector<pair<uint,uint>> edges;
  for(uint indVar1 = 0; indVar1 < _dataset.nColumns(); indVar1 ++) {
    for(uint indVar2 = indVar1+1; indVar2 < _dataset.nColumns(); indVar2 ++) {
      edges.push_back(pair<uint,uint>(indVar1, indVar2));
    }
  }

  vector<Edge> output(edges.size());


  // cout << "n edges: " << edges.size() << endl;

  uint nEdges = 0;

  auto start = chrono::steady_clock::now();

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif

  for(uint ind = 0; ind < edges.size()/*1000000*/; ind ++) {

    string debugStr;

    auto& edge = edges[ind];
    double cor = computeCor(edge.first, edge.second);

    output[ind] = Edge(edge.first, edge.second, cor);

    // debugStr += to_string(edge.first) + " <-> " + to_string(edge.second) + ": ";
    // debugStr += to_string(cor) + "\n";

    #if DEBUG_LOG == 1
      if(nEdges % 500000 == 0) {
        debugStr += to_string(nEdges) + " over " + to_string(edges.size()) + "\n";
        cout << debugStr << endl;
      }
    #endif


    nEdges += 1;
  }

  #if USE_OPENMP == 1
  }
  #endif

  // #if DEBUG_LOG == 1
  auto stop = chrono::steady_clock::now();
  // #endif

  cout << "Computational time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;


  /* output of the computation -> store into a file */
  ofstream file(fileName);
  for(uint ind = 0; ind < output.size(); ind ++) {
    auto edge = output[ind];
    if(abs(output[ind].weight) >= threshold) {
      file << _dataset.getColLabel(edge.indVar1) << " " << _dataset.getColLabel(edge.indVar2) << " " << edge.weight << endl;
    }
  }
  file.close();


}


/*----------------------------------------------------------------------------*/
double LinearGraph::computeCor(uint indVar1, uint indVar2) {

  /* http://www.biostat.ulg.ac.be/pages/Site_r/corr_pearson.html */


  double mean1 = 0.;
  for(uint indVal = 0; indVal < _dataset.nRows(); indVal ++) {
    mean1 += _dataset.getData(indVal, indVar1);
  }
  mean1 /= static_cast<double>(_dataset.nRows());

  double mean2 = 0.;
  for(uint indVal = 0; indVal < _dataset.nRows(); indVal ++) {
    mean2 += _dataset.getData(indVal, indVar2);
  }
  mean2 /= static_cast<double>(_dataset.nRows());

  double num = 0;
  double den1 = 0., den2 = 0;

  for(uint indVal = 0; indVal < _dataset.nRows(); indVal ++) {
    num += (_dataset.getData(indVal, indVar1) - mean1)*(_dataset.getData(indVal, indVar2) - mean2);
    den1 += (_dataset.getData(indVal, indVar1) - mean1)*(_dataset.getData(indVal, indVar1) - mean1);
    den2 += (_dataset.getData(indVal, indVar2) - mean2)*(_dataset.getData(indVal, indVar2) - mean2);
  }

  return num/(sqrt(den1)*sqrt(den2));

}
