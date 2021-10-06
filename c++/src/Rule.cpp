/*!
 * \file Rule.cpp
 * \author S. Buchet
 * \brief implementation of class Rule
 */

#include "Rule.hpp"

#include <algorithm>

using namespace std;

/*----------------------------------------------------------------------------*/
Rule::Rule(Atom head, std::vector<Atom>& body):
_head(head), _body(body)
{

}


/*----------------------------------------------------------------------------*/
bool Rule::subsumes(const Rule& rule) {

  bool result = true;

  if(_body.size() > rule._body.size() || _head != rule._head) {
    result = false;
  } else {
    /* make sure that each atom in the body is also in the other body */
    for(auto atom : _body) {
      if(find(rule._body.begin(), rule._body.end(), atom) == rule._body.end()) {
        result = false;
      }
    }
  }

  return result;

}

/*----------------------------------------------------------------------------*/
bool Rule::match(const State& state) const {

  bool res = true;

  uint index = 0;
  while(res && index < bodyLength()) {
    if(_body[index].second != state[_body[index].first]) {
      res = false;
    } else {
      index += 1;
    }
  }

  return res;
}


/*----------------------------------------------------------------------------*/
uint Rule::matchingScore(const State& state) const {
  uint score = 0;
  for(uint index = 0; index < bodyLength(); index ++) {
    if(_body[index].second != state[_body[index].first]) {
      score += 1;
    }
    index += 1;
  }
  return score;
}


/*----------------------------------------------------------------------------*/
uint Rule::bodyLength() const {
  return static_cast<uint>(_body.size());
}

/*----------------------------------------------------------------------------*/
Atom Rule::getHead() const {
  return _head;
}


/*----------------------------------------------------------------------------*/
Atom Rule::getBody(uint index) const {
  return _body.at(index);
}

/*----------------------------------------------------------------------------*/
const vector<Atom>& Rule::getBody() const {
  return _body;
}
