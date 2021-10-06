/*!
 * \file Combinatorics.cpp
 * \author S. Buchet
 * \brief implementation of Combinatorics
 */


#include "Combinatorics.hpp"

#include <list>
#include <algorithm>

#include <iostream>

using namespace std;

/*----------------------------------------------------------------------------*/
void Combinatorics::choose_kn(const vector<uint>& elements, uint k, vector<std::vector<uint>>& subsets) {

  uint n = static_cast<uint>(elements.size());
  list<uint> pending;
  pending.push_back(0);
  bool backtrack = false;

  while(pending.size() > 0) {

    /* backtracking: find the last elements that can still be inserted */
    while(backtrack) {

      if(pending.empty()) {
        backtrack = false;
      } else {
        if(pending.back() < n-1 && n-pending.back()-1 > k-pending.size()) {
          backtrack = false;
          pending.back() += 1;
        } else {
          pending.pop_back();
        }
      }

    }

    if(pending.size() == k) { /* new subset found */

      /* add the new subset to the list */
      subsets.push_back(vector<uint>(k));
      int ind = 0;
      for(auto it = pending.begin(); it != pending.end(); it ++) {
        subsets.back()[ind] = *it;
        ind ++;
      }
      backtrack = true;

    } else if(!pending.empty()) {

      /* add the next element in the list */
      uint val = pending.back();
      if(n-val >= k-pending.size()) {
        pending.push_back(val+1);
      }

    } else {
      backtrack = true;
    }

  }

}


/*----------------------------------------------------------------------------*/
void Combinatorics::choose_kn_occurences(const vector<uint>& occurences, uint k, vector<vector<pair<uint,uint>>>& subsets) {

  bool backtrack = false;

  /* number of items still available */
  vector<uint> remaining(occurences.size(), 0);
  for(int ind = static_cast<int>(remaining.size())-2; ind >= 0; ind --) {
    remaining[ind] = remaining[ind+1]+occurences[ind+1];
  }

  int new_count = max((int)1, static_cast<int>(k-remaining[0]));

  list<pair<uint,uint>> pending; /* stack used for oredered exploration */
  pending.push_back(pair<uint,uint>(0,new_count));
  uint n_selected = new_count;

  while(!pending.empty()) {

    /*cout << "current stack: ";
    for(auto it = pending.begin(); it != pending.end(); it ++) {
      cout << "(" << it->first << "," << it->second << ")   ";
    }
    cout << endl;
    cout << "n selected: " << n_selected << endl;
    cout << "backtrack ? " << backtrack << endl;*/

    if(backtrack) {

      /* backtracking phase */
      while(backtrack) {
        if(!pending.empty()) {
          pair<uint,uint> top = pending.back();

          if(n_selected < k && top.second < occurences[top.first]) { /* found one element still available */
            pending.back().second ++;
            n_selected ++;
            backtrack = false;
          } else { /* get back to the previous available element */
            pending.pop_back();
            n_selected -= top.second;

            if(top.first < occurences.size()-1 && remaining[top.first] >= k-n_selected) {
              new_count = max(1, static_cast<int>(k-n_selected-remaining[top.first+1]));
              n_selected += new_count;
              pending.push_back(pair<uint, uint>(top.first+1, new_count));
              backtrack = false;
            }
          }

        } else {
          backtrack = false;
        }

      }

    } else if(n_selected == k) {
      subsets.push_back(vector<pair<uint,uint>>(pending.begin(), pending.end()));
      backtrack = true;
    } else { /* exploration phase */

      pair<uint,uint> top = pending.back();

      if(top.first < occurences.size()-1 && remaining[top.first] >= k-n_selected) {
        new_count = max(1, static_cast<int>(k-n_selected-remaining[top.first+1]));
        pending.push_back(pair<uint,uint>(top.first+1, new_count));
        n_selected += new_count;
      } else {
        backtrack = true;
      }

    }

  }

}

/*----------------------------------------------------------------------------*/
void Combinatorics::generatePermutations(unsigned int size, vector<vector<bool>>& permutations) {
  vector<bool> perm(size, false);
  genereatePermutationsRec(0, perm, permutations);
}

/*----------------------------------------------------------------------------*/
void Combinatorics::genereatePermutationsRec(unsigned int index, vector<bool>& perm, vector<vector<bool>>& permutations) {
  if(index >= perm.size()) { /* base case */
    permutations.push_back(perm);
  } else { /* recursive case */
    perm.at(index) = false;
    genereatePermutationsRec(index+1, perm, permutations);
    perm.at(index) = true;
    genereatePermutationsRec(index+1, perm, permutations);
  }
}

/*----------------------------------------------------------------------------*/
void Combinatorics::generateRange(uint length, vector<uint>& index) {
  index.resize(length);
  for(uint ind = 0; ind < length; ind ++) {
    index[ind] = ind;
  }
}
