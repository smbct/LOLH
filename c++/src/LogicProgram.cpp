/*!
 * \file LogicProgram.cpp
 * \author S. Buchet
 * \brief implementation of class LogicProgram
 */

#include "LogicProgram.hpp"
#include "Utils.hpp"

#include <fstream>
#include <iostream>

using namespace std;

/*----------------------------------------------------------------------------*/
LogicProgram::LogicProgram() {

}

/*----------------------------------------------------------------------------*/
uint LogicProgram::nVariables() const {
  return static_cast<uint>(_variables.size());
}

/*----------------------------------------------------------------------------*/
string LogicProgram::varName(uint index) const {
  return _variables[index];
}

/*----------------------------------------------------------------------------*/
uint LogicProgram::varIndex(string varName) const {
  return _varIndex.at(varName);
}

/*----------------------------------------------------------------------------*/
uint LogicProgram::nValues(uint varIndex) const {
  return _nValues[varIndex];
}

/*----------------------------------------------------------------------------*/
void LogicProgram::create(vector<string>& variables, uint nValues) {
  _nValues.resize(variables.size(), nValues);
  _variables = variables;
  _ruleIndexes.resize(variables.size());
  /* variables indexes */
  for(uint ind = 0; ind < _variables.size(); ind ++) {
    _varIndex.insert(pair<string,uint>(_variables[ind], ind));
    _ruleIndexes[ind].resize(nValues);
  }

}

/*----------------------------------------------------------------------------*/
void LogicProgram::addRule(Atom& head, vector<Atom>& body) {
  _ruleIndexes[head.first][head.second].push_back(static_cast<uint>(_rules.size()));
  _rules.push_back(Rule(head, body));
}

/*----------------------------------------------------------------------------*/
void LogicProgram::addVariable(string varName, uint nValues) {
  _variables.push_back(varName);
  _nValues.push_back(nValues);
  _varIndex.insert(pair<string, uint>(varName, _variables.size()-1));
  _ruleIndexes.push_back(vector<vector<uint>>(nValues));
}

/*----------------------------------------------------------------------------*/
void LogicProgram::loadFromFile(string filename) {

  cout << "loading logic program " << filename << endl;

  ifstream file(filename);

  if(file) {

      string line;

      vector<string> tokens;

      while(std::getline(file, line)) {

        tokens.clear();

        if(line.find("VAR") != string::npos) {
          // cout << line << endl;
          Utils::getInstance().split(line, " ",tokens);

          /* variable name */
          _variables.push_back(tokens[1]);
          _nValues.push_back(static_cast<uint>(tokens.size())-2);
          _varIndex.insert(pair<string, uint>(tokens[1], _variables.size()-1));

          /* prepare rule indexes vector */
          _ruleIndexes.push_back(vector<vector<uint>>(_nValues.back()));

          /* variable values */
          // for(uint ind = 2; ind < tokens.size(); ind ++) {
          // }

        } else if(line.find(":-") != string::npos) {

          Utils::getInstance().split(line, " ",tokens);

          Atom head = extractAtom(tokens[0]);

          vector<Atom> body;

          if(tokens[2].size() > 1) {
            for(uint ind = 2; ind < tokens.size(); ind ++) {
              body.push_back(extractAtom(tokens[ind]));
            }
          }

          _ruleIndexes[head.first][head.second].push_back(static_cast<uint>(_rules.size()));
          _rules.push_back(Rule(head, body));

        }


      }

      file.close();
  }

}

/*----------------------------------------------------------------------------*/
void LogicProgram::exportToFile(string fileName) {

  ofstream file(fileName);

  if(file) {

    /* write variables names and values */
    for(uint indVar = 0; indVar < _variables.size(); indVar ++) {
      file << "VAR " << _variables[indVar];
      for(uint val = 0; val < _nValues[indVar]; val ++) {
        file << " " << to_string(val);
      }
      file << endl;
    }

    file << endl;

    /* write rules */
    for(auto& rule: _rules) {

      /* rule head */
      file << _variables[rule.getHead().first] << "(" << rule.getHead().second << ",T)" << " :- ";

      /* rule body */
      for(uint indBody = 0; indBody < rule.bodyLength(); indBody ++) {
        const auto& atom = rule.getBody(indBody);
        file << _variables[atom.first] << "(" << to_string(atom.second) << ",T-1)";
        if(indBody < rule.bodyLength()-1) {
          file << ", ";
        }
      }
      file << "." << endl;

    }

    file.close();

    cout << "writing done!" << endl;

  } else {
    cout << "error" << endl;
  }

}


/*----------------------------------------------------------------------------*/
void LogicProgram::createSubProgram(const LogicProgram& program, uint nVar) {

  vector<bool> selected(program.nVariables(), false);
  for(uint ind = 0; ind < nVar; ind ++) {
    selected[ind] = true;
  }

  _variables.resize(nVar);
  _nValues.resize(nVar);
  for(uint ind = 0; ind < nVar; ind ++) {
    _variables[ind] = program._variables[ind];
    _nValues[ind] = program._nValues[ind];
    _varIndex[_variables[ind]] = ind;
  }

  /* only keep rules for these variables, and make sur there is no condition from another variable */
  for(auto rule : program._rules) {

    if(selected[rule.getHead().first]) { /* make sure this is a selected variable */
      auto body = vector<Atom>();

      /* create a new body containing only selected variables */
      for(auto atom : rule.getBody()) {
        if(selected[atom.first]) {
          body.push_back(atom);
        }
      }

      /* add the new rule */
      _rules.push_back(Rule(rule.getHead(), body));
    }

  }

  /* simplify the new rules */
  simplify();

}



/*----------------------------------------------------------------------------*/
void LogicProgram::getConclusions(const State& state, vector<vector<uint>>& atoms) const {
  atoms.resize(nVariables());
  for(const Rule& rule : _rules) {
    if(rule.match(state)) {
      auto head = rule.getHead();
      vector<uint>& list = atoms.at(head.first);
      if(find(list.begin(), list.end(), head.second) == list.end()) {
        list.push_back(head.second);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
void LogicProgram::simplify() {

  vector<bool> to_keep(_rules.size(), true);

  for(uint indRule = 0; indRule < _rules.size(); indRule ++) {
    uint indRule2 = 0;
    while(indRule2 < _rules.size() && to_keep[indRule]) {
      if(indRule != indRule2 && _rules[indRule2].subsumes(_rules[indRule])) {
        /* if the rules are equal, decide on indexes which one to keep */
        if(_rules[indRule2].bodyLength() != _rules[indRule].bodyLength() || indRule2 < indRule) {
          to_keep[indRule] = false;
        }
      }
      indRule2 += 1;
    }
  }

  /* keep only the selected rules */
  vector<Rule> newRules;
  for(uint ind = 0; ind < _rules.size(); ind ++) {
    if(to_keep[ind]) {
      newRules.push_back(_rules[ind]);
    }
  }

  _rules.clear();
  _rules.insert(_rules.end(), newRules.begin(), newRules.end());


}

/*----------------------------------------------------------------------------*/
string LogicProgram::toString() const {
  string res;

  for(uint indRule = 0; indRule < _rules.size(); indRule ++) {

    const Rule& rule = _rules[indRule];

    res += _variables[rule.getHead().first] + "_" + to_string(rule.getHead().second);
    res += " :-";

    for(uint ind = 0; ind < rule.bodyLength(); ind ++) {
      res += " " + _variables[rule.getBody(ind).first] + "_" + to_string(rule.getBody(ind).second);
    }

    res += "\n";

  }

  return res;
}


/*----------------------------------------------------------------------------*/
pair<uint,uint> LogicProgram::extractAtom(string str) {
  pair<uint,uint> res;
  uint ind1 = static_cast<uint>(str.find("("));
  res.first = _varIndex[str.substr(0,ind1)];
  uint ind2 = static_cast<uint>(str.find(","));
  res.second = std::atoi(str.substr(ind1+1, ind2-ind1).c_str());
  return res;
}
