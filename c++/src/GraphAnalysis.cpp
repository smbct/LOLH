/*!
 * \file GraphAnalysis.cpp
 * \author S. Buchet
 * \brief implementation of class GraphAnalysis
 */

#include "GraphAnalysis.hpp"

#include <iostream>
#include <algorithm>
#include <stack>

using namespace std;

/*----------------------------------------------------------------------------*/
void GraphAnalysis::test() {

  cout << "test graph analysis" << endl;

  string datasetFileName = "dataset/MAIT/binary_logcpm_reduced.csv";

  string transitionsFileName = "dataset/MAIT/transitions_n2_bidir.csv"; /* 5 neighbours */
  // string transitionsFileName = "dataset/MAIT/transitions_bidir_pt.csv"; /* 15 neighbours */

  cout << "graph file: " << transitionsFileName << endl;

  /* load the dataset */
  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  DataFrame<string> transitions;
  transitions.loadFromCSV(transitionsFileName);

  NGraph graph; /* list of successors */

  DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, 1, graph);

  cout << "n vertices: " << graph.size() << endl;
  cout << "n transitions: " << transitions.nRows() << endl;

  /* display the graph */
  // for(uint ind = 0; ind < graph.size(); ind ++) {
  //   cout << ind << ": ";
  //   for(auto& suc: graph[ind]) {
  //     cout << " " << suc;
  //   }
  //   cout << endl;
  // }


  /* compute shortest distance from vert to any other vertex */
  vector<vector<int>> distances(graph.size(), vector<int>(graph.size(), -1));

  for(uint ind = 0; ind < graph.size(); ind ++) {
    GraphAnalysis::dijkstra(ind, graph, distances[ind]);
  }

  int max = 1;
  for(uint ind = 0; ind < graph.size(); ind ++) {
    for(uint ind2 = ind+1; ind2 < graph.size(); ind2 ++) {
      if(distances[ind][ind2] != -1 && distances[ind][ind2] > max) {
        max = distances[ind][ind2];
      }
    }
  }

  cout << "max dist: " << max << endl;

  vector<int> occ(max+1, 0);
  for(uint ind = 0; ind < graph.size(); ind ++) {
    for(uint ind2 = ind+1; ind2 < graph.size(); ind2 ++) {
      if(distances[ind][ind2] >= 0) {
        occ[distances[ind][ind2]] += 1;
      }
    }
  }

  /* display distacnes occurances (histogram) */
  for(uint ind = 0; ind < occ.size(); ind ++) {
    cout << ind << ": " << occ[ind] << endl;
  }

  cout << "unreachable pairs: " << endl;
  uint count = 0;
  for(uint ind = 0; ind < graph.size(); ind ++) {
    for(uint ind2 = ind+1; ind2 < graph.size(); ind2 ++) {
      if(distances[ind][ind2] < 0) {
        // cout << ind << "," << ind2 << " ";
        count += 1;
      }
    }
  }
  // cout << endl;
  cout << count << " pairs" << endl;

  // for(uint ind = 0; ind < graph.size(); ind ++) {
  //   cout << ind << ": " << distances[0][ind] << endl;
  // }

  vector<vector<int>> conComponents;
  GraphAnalysis::connectedComponents(graph, conComponents);

  cout << "n connected components: " << conComponents.size() << endl;

}


/*----------------------------------------------------------------------------*/
void GraphAnalysis::dijkstra(int origin, NGraph& graph,  vector<int>& distances) {

  /* compute shortest distance from vert to any other vertex */
  distances[origin] = 0;

  vector<int> elements(graph.size());
  for(int ind = 0; ind < static_cast<int>(graph.size()); ind ++) {
    elements[ind] = ind;
  }

  const auto& comp = [&distances](const int& left, const int& right) { return distances[left] >= 0 && (distances[right] < 0 || distances[left] < distances[right]); };

  while(!elements.empty()) {

    /* find and remove the smallest element */
    const auto& it = min_element(elements.begin(), elements.end(), comp);
    int elt = *it;

    elements.erase(it);

    for(auto& suc: graph[elt]) {
      int dist2 = 1+distances[elt]; /* distance of the alternative path going though elt */
      if(distances[suc] < 0 || dist2 < distances[suc]) {
        distances[suc] = dist2; /* shorter path found from vert to suc */
      }
    }
  }


}

/*----------------------------------------------------------------------------*/
void GraphAnalysis::connectedComponents(NGraph& graph, vector<vector<int>>& comp) {

  vector<bool> visited(graph.size(), false);
  int n_visited = 0;

  while(n_visited < static_cast<int>(graph.size())) {

    comp.push_back(vector<int>());
    stack<int> pending;
    int ind = 0;
    while(pending.empty() && ind < static_cast<int>(graph.size())) { /* find the first unvisited vertex */
      if(!visited[ind]) {
        pending.push(ind);
      } else {
        ind += 1;
      }
    }

    while(!pending.empty()) { /* compute all its connected components */
      auto top = pending.top();
      pending.pop();
      if(!visited[top]) {
        visited[top] = true;
        n_visited ++;
        comp.back().push_back(top);
        for(auto& suc: graph[top]) {
          if(!visited[suc]) {
            pending.push(suc);
          }
        }
      }
    }

  }



}
