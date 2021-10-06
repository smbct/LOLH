/*!
 * \file DataFrame.cpp
 * \author S. Buchet
 * \brief implementation of class DataFrame
 */

#include <set>
#include <stack>

#include <iostream>
#include <fstream>
#include <sstream>

#include "Utils.hpp"
#include "DataFrame.hpp"

using namespace std;

template class DataFrame<uint>;
template class DataFrame<string>;
template class DataFrame<double>;

/*----------------------------------------------------------------------------*/
template <class T>
DataFrame<T>::DataFrame() {

}

/*----------------------------------------------------------------------------*/
template <class T>
string DataFrame<T>::getColLabel(uint index) const {
  return _columns[index];
}

/*----------------------------------------------------------------------------*/
template <class T>
void DataFrame<T>::getColumnLabels(vector<string>& labels) {
  labels = _columns;
}

/*----------------------------------------------------------------------------*/
template <class T>
bool DataFrame<T>::containsRowIndex(string label) {
  return _rowIndexes.count(label) > 0;
}


/*----------------------------------------------------------------------------*/
template <class T>
uint DataFrame<T>::getColumnIndex(string label) {
  return _colIndexes[label];
}

/*----------------------------------------------------------------------------*/
template <class T>
std::string DataFrame<T>::getRowLabel(uint index) {
  return _rows[index];
}

/*----------------------------------------------------------------------------*/
template <class T>
uint DataFrame<T>::getRowIndex(string label) {
  return _rowIndexes[label];
}

/*----------------------------------------------------------------------------*/
template <class T>
uint DataFrame<T>::nRows() const {
  return static_cast<uint>(_rows.size());
}

/*----------------------------------------------------------------------------*/
template <class T>
uint DataFrame<T>::nColumns() const {
  return static_cast<uint>(_columns.size());
}

/*----------------------------------------------------------------------------*/
template <class T>
T& DataFrame<T>::getData(uint rowIndex, uint colIndex) {
  return _data[rowIndex][colIndex];
}

/*----------------------------------------------------------------------------*/
template <class T>
const vector<T>& DataFrame<T>::getRow(uint rowIndex) const {
  return _data[rowIndex];
}

/*----------------------------------------------------------------------------*/
template <class T>
const std::vector<T>& DataFrame<T>::getRow(std::string rowLabel) const {
  uint index = _rowIndexes.at(rowLabel);
  return getRow(index);
}

/*----------------------------------------------------------------------------*/
template <>
uint DataFrame<uint>::nUnique(uint colInd) const {

  if(!_nVal.empty()) {
    return _nVal[colInd];
  } else {
    uint nUnique = 0;
    for(uint rowIndex = 0; rowIndex < nRows(); rowIndex ++) {
      uint val = _data[rowIndex][colInd]+1;
      if(val > nUnique) {
        nUnique = val;
      }
    }
    return nUnique;
  }

}

/*----------------------------------------------------------------------------*/
template <>
uint DataFrame<string>::nUnique(uint colInd) const {
  set<string> unique;
  for(uint rowIndex = 0; rowIndex < nRows(); rowIndex ++) {
    const string& val = _data[rowIndex][colInd];
    unique.insert(val);
  }
  return static_cast<uint>(unique.size());
}

/*----------------------------------------------------------------------------*/
template <typename T>
map<T, uint> DataFrame<T>::uniqueCount(uint colInd) const {
  map<T, uint> counts;
  for(uint rowIndex = 0; rowIndex < nRows(); rowIndex ++) {
    const T& val = _data[rowIndex][colInd];
    if(counts.find(val) == counts.end()) {
      counts.insert(pair<T,uint>(val,1));
    } else {
      counts[val] += 1;
    }
  }
  return counts;
}

/*----------------------------------------------------------------------------*/
template<typename T>
map<T, uint> DataFrame<T>::uniqueCount(uint colInd, const std::vector<uint>& rows) const {
  map<T, uint> counts;
  for(uint rowIndex : rows) {
    const T& val = _data[rowIndex][colInd];
    if(counts.find(val) == counts.end()) {
      counts.insert(pair<T,uint>(val,1));
    } else {
      counts[val] += 1;
    }
  }
  return counts;
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::uniqueCount(uint colInd, const vector<uint>& rows, vector<uint>& counts) const {
  for(uint ind = 0; ind < rows.size(); ind ++) {
    counts[_data[rows[ind]][colInd]] ++;
  }
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::uniqueCount(const vector<bool>& classVector, vector<vector<uint>>& posCounts, vector<vector<uint>>& negCounts) const {
  // for(uint indRow = 0; indRow < nRows(); indRow ++) {
  //   for(uint indCol = 0; indCol < nColumns(); indCol ++) {
  //     if(classVector[indRow]) {
  //       posCounts[indCol][_data[indRow][indCol]] ++;
  //     } else {
  //       negCounts[indCol][_data[indRow][indCol]] ++;
  //     }
  //   }
  // }

  // for(uint indCol = 0; indCol < nColumns(); indCol ++) {
  //   for(uint indRow = 0; indRow < nRows(); indRow ++) {
  //     if(classVector[indRow]) {
  //       posCounts[indCol][_data[indRow][indCol]] ++;
  //     } else {
  //       negCounts[indCol][_data[indRow][indCol]] ++;
  //     }
  //   }
  // }

  if(posCounts.empty()) {
    posCounts.resize(_columns.size(), vector<uint>(_maxVal, 0));
    negCounts.resize(_columns.size(), vector<uint>(_maxVal, 0));
  }

  uint indRow = 0;
  for(auto& elt : _data) {
    for(uint indCol = 0; indCol < nColumns(); indCol ++) {
      if(classVector[indRow]) {
        posCounts[indCol][elt[indCol]] ++;
      } else {
        negCounts[indCol][elt[indCol]] ++;
      }
    }
    indRow += 1;
  }

}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::uniqueCount(const vector<uint>& colIndexes, const vector<bool>& classVector, vector<vector<uint>>& posCounts, vector<vector<uint>>& negCounts) const {

  uint indRow = 0;
  for(auto& elt : _data) {
    for(const uint& indCol: colIndexes) {
      if(classVector[indRow]) {
        posCounts[indCol][elt[indCol]] ++;
      } else {
        negCounts[indCol][elt[indCol]] ++;
      }
    }
    indRow += 1;
  }

}


/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::uniqueCount(const vector<uint>& rows, vector<std::vector<uint>>& counts) const {

  if(counts.empty()) {
    counts.resize(_columns.size(), vector<uint>(_maxVal, 0));
  }

  for(const uint& row : rows) {
    for(uint indCol = 0; indCol < nColumns(); indCol ++) {
      counts[indCol][_data[row][indCol]] ++;
    }
  }
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::uniqueCount(vector<std::vector<uint>>& counts) const {

  if(counts.empty()) {
    counts.resize(_columns.size(), vector<uint>(_maxVal, 0));
  }

  for(uint rowIndex = 0; rowIndex < nRows(); rowIndex ++) {
    for(uint indCol = 0; indCol < nColumns(); indCol ++) {
      counts[indCol][_data[rowIndex][indCol]] ++;
    }
  }
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::computeUniqueVal() {
  _nVal.resize(_columns.size());
  _maxVal = 0;
  for(uint indCol = 0; indCol < _nVal.size(); indCol ++) {
    _nVal[indCol] = 0;
    for(uint rowIndex = 0; rowIndex < nRows(); rowIndex ++) {
      uint val = _data[rowIndex][indCol]+1;
      if(val > _nVal[indCol]) {
        _nVal[indCol] = val;
      }
      if(val > _maxVal) {
        _maxVal = val;
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
template <class T>
void DataFrame<T>::createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<string>& transitions, NGraph& graph) {
  /* compute a list of successors (indexes) from hte transitions */
  graph.resize(dataframe.nRows());
  /* create successor lists */
  for(uint trIndex = 0; trIndex < transitions.nRows(); trIndex ++) {
    uint indLeft = dataframe.getRowIndex(transitions.getRow(trIndex)[0]);
    uint indRight = dataframe.getRowIndex(transitions.getRow(trIndex)[1]);
    graph[indLeft].push_back(indRight);
  }
}

/*----------------------------------------------------------------------------*/
template <class T>
void DataFrame<T>::createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<string>& transitions, uint delay, NGraph& graph) {

  /* initialization of the transition graph structure */
  graph.resize(dataframe.nRows());

  /* compute a list of successors (indexes) from the transitions */
  vector<vector<uint>> successors(dataframe.nRows());
  for(uint trIndex = 0; trIndex < transitions.nRows(); trIndex ++) {
    uint indLeft = dataframe.getRowIndex(transitions.getRow(trIndex)[0]);
    uint indRight = dataframe.getRowIndex(transitions.getRow(trIndex)[1]);
    successors[indLeft].push_back(indRight);
  }

  stack<pair<uint,uint>> pending; /* data structure to explore the neighborhood */

  /* compute the transitions structure with some delay */
  for(uint indRow = 0; indRow < dataframe.nRows(); indRow ++) {

    /* compute all successors at a given distance (delay) */
    pending.push(pair<uint,uint>(indRow,0));

    while(!pending.empty()) {
      auto top = pending.top();
      pending.pop();
      if(top.second == delay) {
        graph[indRow].push_back(top.first); /* new successor found */
      } else if(top.second < delay) {
        for(auto& indSuc: successors[top.first]) { /* add all successors to the stack */
          pending.push(pair<uint,uint>(indSuc, top.second+1)); /* distance += 1 */
        }
      }
    }

    /* remove duplicates */

    /* keep only unique successors */
    sort(graph[indRow].begin(), graph[indRow].end());
    auto it = unique(graph[indRow].begin(), graph[indRow].end());
    graph[indRow].resize(distance(graph[indRow].begin(), it));

  }


}

/*----------------------------------------------------------------------------*/
template <class T>
void DataFrame<T>::createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<string>& transitions, DataFrame<string>& ordering, uint delay, NGraph& graph) {

  cout << "hello ordering :)" << endl;

  cout << "n theoretical transitions: " << dataframe.nRows()*(dataframe.nRows()-1) << endl;

  /* create a dictionary between original integer index of the cell and pseudtime ordering index */
  vector<uint> pt_index(dataframe.nRows());
  for(uint ind = 0; ind < ordering.nRows(); ind ++) {
    uint df_ind = dataframe.getRowIndex(ordering.getData(ind, 0));
    uint ordering_index = std::atoi(ordering.getData(ind, 1).c_str());
    pt_index[df_ind] = ordering_index;
  }

  // for(uint ind = 0; ind < 10; ind ++) {
  //   cout << ind << " vs " << pt_index[ind] << endl;
  // }

  /* initialization of the transition graph structure */
  graph.resize(dataframe.nRows());

  /* compute a list of successors (indexes) from the transitions */
  vector<vector<uint>> successors(dataframe.nRows());
  for(uint trIndex = 0; trIndex < transitions.nRows(); trIndex ++) {
    uint indLeft = dataframe.getRowIndex(transitions.getRow(trIndex)[0]);
    uint indRight = dataframe.getRowIndex(transitions.getRow(trIndex)[1]);
    successors[indLeft].push_back(indRight);
  }

  stack<pair<uint,uint>> pending; /* data structure to explore the neighborhood */

  /* compute the transitions structure with some delay */
  for(uint indRow = 0; indRow < dataframe.nRows(); indRow ++) {

    /* compute all successors at a given distance (delay) */
    pending.push(pair<uint,uint>(indRow,0));
    while(!pending.empty()) {
      auto top = pending.top();
      pending.pop();

      if(top.second <= delay) {

        /* here add a constraint on the pseudotime */
        int ordering_dist =  pt_index[top.first] - pt_index[indRow];
        if((ordering_dist >= 1000 && ordering_dist <= 1300)/* || (ordering_dist >= -200 && ordering_dist <= -50)*/) {
          graph[indRow].push_back(top.first); /* new successor found */
        }


      }

      if(top.second < delay) {
        for(auto& indSuc: successors[top.first]) { /* add all successors to the stack */
          pending.push(pair<uint,uint>(indSuc, top.second+1)); /* distance += 1 */
        }
      }
    }

    /* remove duplicates */
    sort(graph[indRow].begin(), graph[indRow].end());
    auto it = unique(graph[indRow].begin(), graph[indRow].end());
    graph[indRow].resize(distance(graph[indRow].begin(), it));


  }

  uint nTransitions = 0;
  for(auto& elt : graph) {
    nTransitions += static_cast<int>(elt.size());
  }
  cout << "n transitions: " << nTransitions << endl;

  /******************************************************/
  /* compute pseudotime difference over all transitions */
  ofstream pt_out("pt_out.txt");
  if(pt_out) {
    for(uint ind = 0; ind < graph.size(); ind ++) {
      for(auto& suc: graph[ind]) {
        pt_out << static_cast<int>(pt_index[suc])-static_cast<int>(pt_index[ind]) << endl;
      }
    }
    pt_out.close();
  }

}


/*----------------------------------------------------------------------------*/
template <class T>
void DataFrame<T>::loadFromCSV(string filename) {

  cout << "loading " << filename << "..." << endl;

  ifstream file(filename);

  if(file) {

    string line;
    string token;

    /* read columns names */
    getline(file, line);
    istringstream stream(line);
    std::getline(stream, token, ',');
    while(std::getline(stream, token, ',')) {
      _columns.push_back(token);
      _colIndexes.insert(pair<string, uint>(token, _columns.size()-1));
    }

    /* read each sample of the dataset */
    vector<T> sample(_columns.size());

    while(!file.eof()) {

      getline(file, line);

      // cout << "*" << line << "*" << endl;

      if(line.size() > 0) {

        istringstream stream2(line);
        std::getline(stream2, token, ',');

        /* read the row index label */
        _rows.push_back(token);
        _rowIndexes.insert(pair<string, uint>(_rows.back(), _rows.size()-1));

        readSample(stream2, sample);

        /* add a new line inot the matrix */
        _data.push_back(sample);

      }

    }

  } else {
    cout << "error opening the csv file" << endl;
  }


}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<uint>::readSample(istringstream& stream, vector<uint>& sample) {
  string token;
  uint colInd = 0;
  while(std::getline(stream, token, ',')) {
    sample[colInd] = std::atoi(token.c_str());
    colInd += 1;
  }
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<string>::readSample(istringstream& stream, vector<string>& sample) {
  string token;
  uint colInd = 0;
  while(std::getline(stream, token, ',')) {
    sample[colInd] = token;
    colInd += 1;
  }
}

/*----------------------------------------------------------------------------*/
template<>
void DataFrame<double>::readSample(istringstream& stream, vector<double>& sample) {
  string token;
  uint colInd = 0;
  while(std::getline(stream, token, ',')) {
    sample[colInd] = std::stod(token.c_str()) ;
    colInd += 1;
  }
}

/*----------------------------------------------------------------------------*/
template <class T>
string DataFrame<T>::toString() {

  string res = "";

  cout << "columns: ";

  for(unsigned int indCol = 0; indCol < _columns.size(); indCol ++) {
    cout << _columns.at(indCol) << " ; ";
  }

  cout << endl << endl << "rows: ";

  for(uint indRow = 0; indRow < _rows.size(); indRow ++) {
    cout << _rows.at(indRow) << " ; ";
  }

  cout << endl << endl;

  cout << "data: " << endl;

  for(uint indRow = 0; indRow < _rows.size(); indRow ++) {
    for(uint indCol = 0; indCol < _columns.size(); indCol ++) {
      cout << _data.at(indRow).at(indCol) << " ";
    }
    cout << endl;
  }

  return res;
}
