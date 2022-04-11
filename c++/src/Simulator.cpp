/*!
 * \file Simulator.cpp
 * \author S. Buchet
 * \brief implementation of class Simulator
 */

#include "Simulator.hpp"

#include <iostream>

using namespace std;

/*----------------------------------------------------------------------------*/
Simulator::Simulator(const LogicProgram& program):
_program(program)
{

}

/*----------------------------------------------------------------------------*/
void Simulator::generateAllStates(vector<State>& states) {

  State pending(_program.nVariables(), 0);

  bool stop = false;

  while(!stop) {

    states.push_back(pending);

    /* find the first 0 from the right */
    uint ind = static_cast<uint>(pending.size())-1;
    bool found = false;

    while(!stop && !found) {

      pending[ind] = 1-pending[ind];

      if(pending[ind] == 1) {
        found = true;
      } else if(ind > 0) {
        ind -= 1;
      } else {
        stop = true;
      }

    }

    if(!found) {
      stop = true;
    }

  }

}

/*----------------------------------------------------------------------------*/
void Simulator::generateAllSuccessors(const State& state, Semantics semantics, vector<State>& successors) {

  if(semantics == Semantics::SYNCHRONOUS) {
    generateSynchronousSuc(state, successors);
  } else if(semantics == Semantics::ASYNCHRONOUS) {
    generateAsynchronousSuc(state, successors);
  }


}

/*----------------------------------------------------------------------------*/
void Simulator::generateSynchronousSuc(const State& state, vector<State>& successors) {

  vector<vector<uint>> atoms;

  _program.getConclusions(state, atoms);

  bool empty = false;
  for(uint ind = 0; ind < atoms.size(); ind ++) {
    if(atoms[ind].size() == 0) {
      empty = true;
    }
  }

  if(!empty) {

    vector<uint> pending(atoms.size(), 0);

    bool stop = false;

    uint ind;

    while(!stop) {

      /* add the last successor to the list */
      successors.push_back(vector<uint>(atoms.size()));
      for(uint pind = 0; pind < successors.back().size(); pind ++) {
        successors.back()[pind] = atoms[pind][pending[pind]];
      }

      // cout << "last: ";
      // for(uint ind = 0; ind < successors.back().size(); ind ++) {
      //   cout << successors.back()[ind];
      // }
      // cout << endl;

      ind = static_cast<uint>(pending.size())-1;
      bool found = false;
      while(!stop && !found) {

        if(atoms[ind].size() > 1) {
          if(pending[ind] > 0 && pending[ind] == atoms[ind].size()-1) {
            pending[ind] = 0;
            found = true;
            if(ind == 0) {
              stop = true;
            }
          } else {
            pending[ind] += 1;
          }

          if(ind > 0) {
            ind -= 1;
          } else {
            stop = true;
          }
        } else if(ind > 0) {
          ind -= 1;
        } else {
          stop = true;
        }

      }
      // cout << endl;

    }

    // for(uint ind = 0; ind < atoms.size(); ind ++) {
    //   cout << ind << " : ";
    //   for(auto val : atoms[ind]) {
    //     cout << val << " ";
    //   }
    //   cout << endl;
    // }
    // cout << endl << endl;
  }
  // cout << "n suc: " << successors.size() << endl;

}

/*----------------------------------------------------------------------------*/
void Simulator::generateAsynchronousSuc(const State& state, vector<State>& successors) {

  vector<vector<uint>> atoms;

  _program.getConclusions(state, atoms);

  for(uint ind = 0; ind < atoms.size(); ind ++) {

    /* create one successor per conclusion */
    for(auto value : atoms[ind]) {
      successors.push_back(state);
      successors.back()[ind] = value;
    }

  }

}
