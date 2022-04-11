/*!
 * \file TransitionGraph.cpp
 * \author S. Buchet
 * \brief implementation of class TransitionGraph
 */

#include "TransitionGraph.hpp"

#include <iostream>
#include <fstream>

using namespace std;


/*----------------------------------------------------------------------------*/
TransitionGraph::TransitionGraph(const LogicProgram& program):
_program(program)
{

}

/*----------------------------------------------------------------------------*/
void TransitionGraph::exportCSV(string fileName) {

  /* create the vertices file */
  ofstream vertFile(fileName + "_vert.csv");

  if(vertFile) {

    for(uint varInd = 0; varInd < _program.nVariables(); varInd ++) {
      vertFile << "," << _program.varName(varInd);
    }
    vertFile << endl;

    for(uint stateInd = 0; stateInd < _states.size(); stateInd ++) {
      vertFile << "s_" + to_string(stateInd);
      for(uint varInd = 0; varInd < _program.nVariables(); varInd ++) {
        vertFile << "," << _states[stateInd][varInd];
      }
      vertFile << endl;
    }


    vertFile.close();
  }

  /* create the edges file */
  ofstream vertFile2(fileName + "_edge.csv");

  if(vertFile2) {

    vertFile2 << ",T-1,T" << endl;

    uint ind_tr = 0;
    for(uint ind = 0; ind < _transitions.size(); ind ++) {
      for(uint ind2 = 0; ind2 < _transitions[ind].size(); ind2 ++) {
        vertFile2 << "t_" + to_string(ind_tr) << "," << ind << "," << _transitions[ind][ind2] << endl;
        ind_tr ++;
      }
    }

    vertFile2.close();
  }

}

/*----------------------------------------------------------------------------*/
void TransitionGraph::generateCompleteGraph(Semantics semantics) {



  /* generate all the states */
  Simulator simulator(_program);
  simulator.generateAllStates(_states);

  /* prepare list of transitions */
  _transitions.resize(_states.size());

  for(uint ind = 0; ind < _states.size(); ind ++) {
    _stateIndexes[_states[ind]] = ind;
  }

  // for(auto state : _states) {
  //   for(uint ind = 0; ind < static_cast<uint>(state.size()); ind ++) {
  //     cout << state[ind];
  //   }
  //   cout << endl;
  // }

  /* generate all the transitions */
  std::vector<State> suc;
  for(const State& state : _states) {

    uint ind1 = _stateIndexes[state];

    suc.clear();

    simulator.generateAllSuccessors(state, semantics, suc);

    /* add the edges to the graph */
    for(auto& state_suc : suc) {
      uint ind2 = _stateIndexes[state_suc];
      _transitions[ind1].push_back(ind2);
      // cout << "(" << ind1 << "," << ind2 << ")" << endl;
    }

  }

  // for(uint ind1 = 0; ind1 < _states.size(); ind1 ++) {
  //   cout << "nb " << ind1 << " : " << _transitions[ind1].size() << endl;
  // }

}
