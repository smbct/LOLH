/*!
 * \file TransitionEmbeddingc.cpp
 * \author S. Buchet
 * \brief implementation of a class TransitionEmbedding
 */

#include "TransitionEmbedding.hpp"

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>

#include "Utils.hpp"
#include "Combinatorics.hpp"

#include "Constants.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
RegulatoryGraph::RegulatoryGraph() {

}

/*----------------------------------------------------------------------------*/
uint RegulatoryGraph::getAtomIndex(Atom atom) {
  uint index = 0;
  auto it = _atomIndexes.find(atom);
  if(it == _atomIndexes.end()) {
    _atoms.push_back(atom);
    index = static_cast<uint>(_atoms.size()-1);
    _atomIndexes.insert(pair<Atom, uint>(atom, index));

    _proportions.push_back(Score());
    _predecessors.push_back(vector<pair<uint,uint>>()); /* for each atom: a list of predecessors: atom index plus index of the edge for this atom */

    /* add the place for edges and reversed edges */
    _successors.push_back(vector<uint>());
    _discreteScores.push_back(vector<Score>());
    _relativeScores.push_back(vector<double>());

  } else {
    index = it->second;
  }
  return index;
}

/*----------------------------------------------------------------------------*/
void RegulatoryGraph::loadFromFile(string fileName) {

  ifstream file(fileName);

  if(file) {

    string line;
    while(getline(file, line)) {

      vector<string> tokens;
      Utils::getInstance().split(line, " ", tokens);

      if(tokens.size() >= 4) {

        /* read the atom */
        Atom atom(tokens[0], atoi(tokens[1].c_str()));
        uint atomIndex = getAtomIndex(atom);

        /* read the proportions */
        Score proportion(atoi(tokens[2].c_str()), atoi(tokens[3].c_str()));
        _proportions.at(atomIndex) = proportion;

        uint token_ind = 4;
        while(token_ind < tokens.size()) {

          /* extract data */
          Atom connectedAtom(tokens[token_ind], atoi(tokens[token_ind+1].c_str()));
          uint connectedAtomIndex = getAtomIndex(connectedAtom);

          /* read the relative edge score */
          double relativeScore = stod(tokens[token_ind+2]);

          /* read the disctete score */
          pair<uint,uint> discreteScore(atoi(tokens[token_ind+3].c_str()), atoi(tokens[token_ind+4].c_str()));

          /* add the edge: reversed direction */
          _successors[connectedAtomIndex].push_back(atomIndex);
          _discreteScores[connectedAtomIndex].push_back(discreteScore);
          _relativeScores[connectedAtomIndex].push_back(relativeScore);

          /* add the reversed edge */
          uint indEdge = static_cast<uint>(_successors[connectedAtomIndex].size())-1;
          _predecessors[atomIndex].push_back(pair<uint,uint>(connectedAtomIndex, indEdge));

          /********************************************************************/
          /* experimentation: add the normal edge as well */
          /* add the edge: reversed direction */
          // _successors[atomIndex].push_back(connectedAtomIndex);
          // _discreteScores[atomIndex].push_back(discreteScore);
          // _relativeScores[atomIndex].push_back(relativeScore);
          // indEdge = static_cast<uint>(_successors[atomIndex].size())-1;
          // _predecessors[connectedAtomIndex].push_back(pair<uint,uint>(atomIndex, indEdge));
          /********************************************************************/

          /* go to the next edge */
          token_ind += 5;

        }

      }

    }

    file.close();


    /* check the content of the graph */
    // Atom atom("CMIP", 0);
    // uint indAtom = _atomIndexes[atom];
    // cout << "atom: " << atom.first << "_" << atom.second << endl;
    // cout << "edges: " << endl;
    // for(uint indEdge = 0; indEdge < _successors[indAtom].size(); indEdge ++) {
    //   Atom suc = _atoms[_successors[indAtom][indEdge]];
    //   cout << suc.first << "_" << suc.second << " ";
    //   cout << _relativeScores[indAtom][indEdge] << " ";
    //   cout << endl;
    // }
    // cout << "reversed edges: " << endl;
    // for(uint redgeInd = 0; redgeInd < _predecessors[indAtom].size(); redgeInd ++) {
    //   uint indAtomPred = _predecessors[indAtom][redgeInd].first;
    //   uint indEdge = _predecessors[indAtom][redgeInd].second;
    //   Atom prec = _atoms[indAtomPred];
    //   cout << prec.first << "_" << prec.second << " ";
    //   cout << _discreteScores[indAtomPred][indEdge].first << "," << _discreteScores[indAtomPred][indEdge].second << " ";
    //   cout << _relativeScores[indAtomPred][indEdge] << endl;
    // }

  } else {
    cout << "unable to open " << fileName << endl;
  }

}

/*----------------------------------------------------------------------------*/
void RegulatoryGraph::computePartialSucessor(DataFrame<uint>& matrix, string& cell, double threshold, vector<Atom>& successors) {

    for(uint indAtom = 0; indAtom < _atoms.size(); indAtom ++) {

      /* make sure there are enough predecessors in the graph for an accurate prediction */

      if(_predecessors[indAtom].size() >= 2) {

        double score = 0;

        /* iterate over the predecessors */
        for(uint indPred = 0; indPred < _predecessors[indAtom].size(); indPred ++) {

          uint indAtomPred = _predecessors[indAtom][indPred].first;
          Atom atomPred = _atoms[indAtomPred];
          uint indEdgePred = _predecessors[indAtom][indPred].second;

          /* check if this predecessor matches with the content of the cell */
          uint cellValue = matrix.getData(matrix.getRowIndex(cell), matrix.getColumnIndex(atomPred.first));
          if(cellValue == atomPred.second) {
            score += _relativeScores[indAtomPred][indEdgePred];
          } else {
            score -= _relativeScores[indAtomPred][indEdgePred];
          }

        }


        score /= static_cast<double>(_predecessors[indAtom].size());

        if(score >= threshold) { /* add the atom as a partial successor */
          successors.push_back(_atoms[indAtom]);
        }

      }

    }

    // cout << "successors: ";
    // for(auto& elt: successors) {
    //   cout << elt.first << "_" << elt.second << ", ";
    // }
    // cout << endl;

}

/*----------------------------------------------------------------------------*/
TransitionEmbedding::TransitionEmbedding(string matrixFileName, string embeddingFileName, string graphFileName, string neighborhoodGraphFileName) {

  cout << "hello transition embedding :)" << endl;

  /* load the discrete expression matrix */
  _matrix.loadFromCSV(matrixFileName);

  loadEmbeddingCoord(embeddingFileName);

  loadNeighbours(neighborhoodGraphFileName);

  _regulatoryGraph.loadFromFile(graphFileName);

  computeTransitionProbabilities();

  computeTransitionEmbedding();


  saveTransitionEmbedding("../../Learning/IMAGINE_dataset/CMSB/Tcells/transitionEmbedding_dist4_pn0_-1.txt");



}

/*----------------------------------------------------------------------------*/
void TransitionEmbedding::computeTransitionProbabilities() {

  #if DEBUG_LOG == 1
    cout << "computation of the transition probabilities" << endl;
  #endif

  double sigma = 2.;

  uint nCells = 0;

  uint nCellsMax = _matrix.nRows();
  #if DEBUG_REDUCTION == 1
    nCellsMax = NVAR_DEBUG;
  #endif

  /* compute on a subset of the cells */
  vector<uint> cellIndexes;
  Combinatorics::generateRange(_matrix.nRows(), cellIndexes);
  #if SAMPLING == true
    nCellsMax = static_cast<uint>(floor(static_cast<double>(_matrix.nRows())*SAMPLE_SIZE));
    Utils::getInstance().shuffle(cellIndexes);
    cellIndexes.resize(nCellsMax);
  #else
    // nCellsMax = _matrix.nRows();
  #endif

  /* initialize the transition probabilities vector for all the cells */
  for(uint indCellIndexes = 0; indCellIndexes < nCellsMax; indCellIndexes ++) {
    /* get the cell label */
    string cell = _matrix.getRowLabel(cellIndexes[indCellIndexes]);
    /* insert the transition probability list associated */
     _transitionProbabilites.insert(pair<string, vector<double>>(cell, vector<double>()));
  }

  /* DEBUG */
  /* initialisation for parallel computation -> take care of empty strings */
  vector<string> file_content(nCellsMax);
  /* DEBUG */

  // geneIndexes.resize()

  #if USE_OPENMP == 1
  #pragma omp parallel num_threads(N_THREADS)
  {
  #pragma omp for
  #endif
  for(uint indCellIndexes = 0; indCellIndexes < nCellsMax; indCellIndexes ++) {

    nCells ++;
    if(nCells%100 == 0) {
      string line = "n cells done yet: " + to_string(nCells) + " over " + to_string(nCellsMax) +"\n";
      cout <<  line;
    }

    uint indRow = cellIndexes[indCellIndexes];

    vector<RegulatoryGraph::Atom> successor;
    string cell = _matrix.getRowLabel(indRow);

    // cout << "cell: " << cell << endl;

    /* create a discrete sorted vector */
    double threshold = 0.1; /* 0.3 */
    _regulatoryGraph.computePartialSucessor(_matrix, cell, threshold, successor); // -> most of the time spent here

    // cout << "n successors: " + to_string(successor.size()) << endl;

    /* TODO: check partial vector size */
    if(successor.size() >= 10) { /* enough genes have been predicted */

      /* successor is sorted according the gene indexes in the single cell matrix */
      auto comp = [this](RegulatoryGraph::Atom left, RegulatoryGraph::Atom right) { return _matrix.getColumnIndex(left.first) < _matrix.getColumnIndex(right.first); };
      sort(successor.begin(), successor.end(), comp);

      /* compute the partial vector of successor discrete values */
      vector<uint> successorVal(successor.size());
      for(uint indSuc = 0; indSuc < successor.size(); indSuc ++) {
        successorVal[indSuc] = successor[indSuc].second;
      }


      /* DEBUG */
      string file_line = cell + " " + to_string(successor.size());
      /* compute similarities with the same cell */
      {
        vector<uint> partialCell(successor.size());
        for(uint indSuc = 0; indSuc < successor.size(); indSuc ++) {
          partialCell[indSuc] = _matrix.getData(_matrix.getRowIndex(cell), _matrix.getColumnIndex(successor[indSuc].first));
        }
        uint similarity = successorVal.size() - Utils::getInstance().manDistance(successorVal, partialCell);
        file_line += " " + to_string(similarity);
      }
      /* DEBUG */


      /* iterate over all the neighbours */
      auto& neighbours = _neighbours[cell];
      // auto& probaIt = _transitionProbabilites.insert(pair<string, vector<double>>(cell, vector<double>(neighbours.size()))).first->second; /* warining: transitions are not inserted if the partial successor does not exist */
      /* TODO warining potentially dangerous on multi core exacution */

      auto transitionMapIterator = _transitionProbabilites.find(cell);
      auto& probaIt = transitionMapIterator->second;
      probaIt.resize(neighbours.size()); /* oneprobability per neighbour */

      double sumScore = 0;

      for(uint indNeighbour = 0; indNeighbour < neighbours.size(); indNeighbour ++) {

        string neighbourCell = neighbours[indNeighbour];
        uint indRowNeighbour = _matrix.getRowIndex(neighbourCell);

        /* make a partial neighbour vector */
        vector<uint> partialNeighbour(successor.size());
        for(uint indSuc = 0; indSuc < successor.size(); indSuc ++) {
          partialNeighbour[indSuc] = _matrix.getData(indRowNeighbour, _matrix.getColumnIndex(successor[indSuc].first));
        }

        /* compute similarities between the partial prediction and the neighbours */

        // double similarity = Utils::getInstance().dotProduct(successorVal, partialNeighbour)/(Utils::getInstance().euclideanNorm(successorVal)*Utils::getInstance().euclideanNorm(partialNeighbour));

        /* Manhattan distance */
        double similarity = static_cast<double>(successorVal.size()) - static_cast<double>(Utils::getInstance().manDistance(successorVal, partialNeighbour)); /* high probability when the number of difference is 0 */

        /* Uniform value */
        // double similarity = 1.;

        /* Random value */
        // double similarity = Utils::getInstance().randInt(0,1);

        /* DEBUG */
        file_line += " " + to_string(static_cast<uint>(similarity));
        /* DEBUG */

        similarity = similarity/static_cast<double>(successorVal.size());

        /* application of a kernel function */
        // similarity = exp( similarity/ (sigma*sigma));

        probaIt[indNeighbour] = similarity;

        /* compute the sum to normalize (<-> probability distribution) */
        sumScore += similarity;

      }

      // cout << "probabilities: ";

      /* normalize the probabilities */
      for(uint indNeighbour = 0; indNeighbour < neighbours.size(); indNeighbour ++) {
        probaIt[indNeighbour] /= sumScore;
        // probaIt[indNeighbour] = 1./neighbours.size();

      }

      /* DEBUG */
      file_content[indCellIndexes] = file_line;
      file_line = " ";
      /* DEBUG */

    }

  }

  #if USE_OPENMP == 1
  }
  #endif

  /* DEBUG */
  ofstream file("../../Learning/IMAGINE_dataset/CMSB/Tcells/transition_probabilities_dist4_pn0_-1_Tcells.txt");
  if(file) {
    for(auto& line: file_content) {
      if(!line.empty()) {
        file << line << endl;
      }
    }
    file.close();
  }
  /* DEBUG */

  // #if DEBUG_LOG == 1
  //   cout << "end transition probabilities computation" << endl;
  // #endif

}

/*----------------------------------------------------------------------------*/
void TransitionEmbedding::computeTransitionEmbedding() {

  #if DEBUG_LOG == 1
    cout << "computation of the transition embedding" << endl;
  #endif

  uint nCell = 0;

  /* replace by iterations over the transition probabilities */
  for(auto it = _transitionProbabilites.begin(); it != _transitionProbabilites.end(); it ++) {

    string cell = it->first;

    /* vector of transition probabilities for this cell */
    vector<double>& probabilities = it->second;

    /* make sure the probabilites are not empty, otherwise, no prediction could be performed for that cell */
    if(!probabilities.empty()) {

      // if(nCell%100 == 0) {
      //   cout << "n cells done yet: " << nCell << endl;
      // }

      // cout << cell << endl << endl;

      /* neighbours and transition probabilities iterators for this cell */
      map<string, vector<string>>::iterator neighboursIterator = _neighbours.find(cell);

      /* make sure the probability distribution over the neighbours has been computed */
      if(neighboursIterator != _neighbours.end()) {

        /* add the new transition vector */
        pair<double, double>& transitionVector = _transitionEmbedding[cell];
        transitionVector.first = 0.;
        transitionVector.second = 0.;

        /* get the embedding coordinates of the cell */
        pair<double, double> cellCoord = _embedding[cell];

        vector<string>& neighbours = neighboursIterator->second;

        double correction = 1./static_cast<double>(neighbours.size());

        /* iterate over all neighbours to compute the expected direction */
        for(uint indNeighbour = 0; indNeighbour < neighbours.size(); indNeighbour ++) {

          string neighbour = neighbours[indNeighbour];
          pair<double, double> neighbourCoord = _embedding[neighbour];

          /* compute and normalize the direction vector */
          pair<double, double> direction = pair<double, double>(neighbourCoord.first-cellCoord.first, neighbourCoord.second-cellCoord.second);
          Utils::getInstance().normalize(direction); /* the vector is normalized */

          transitionVector.first += direction.first*(probabilities[indNeighbour]-correction);
          transitionVector.second += direction.second*(probabilities[indNeighbour]-correction);

        }


        // cout << endl << endl << endl;

        /* test: select a random direction */
        // if(neighbours.size() > 0) {
        //   string neighbour = neighbours[Utils::getInstance().randUInt(0, neighbours.size()-1)];
        //   pair<double, double> neighbourCoord = _embedding[neighbour];
        //   pair<double, double> direction = pair<double, double>(neighbourCoord.first-cellCoord.first, neighbourCoord.second-cellCoord.second);
        //   transitionVector.first = direction.first;
        //   transitionVector.second = direction.second;
        //   // cout << atan2(direction.second, direction.first) << ", ";
        //   // cout << endl << endl;
        // } else {
        //     transitionVector.first = 0;
        //     transitionVector.second = 0;
        // }
        //
        // transitionVector.first = 1;
        // transitionVector.second = 1;

        // cout << transitionVector.first << " ; " << transitionVector.second << endl;

      }


    }


  }

}

/*----------------------------------------------------------------------------*/
void TransitionEmbedding::saveTransitionEmbedding(string fileName) {

  #if DEBUG_LOG == 1
    cout << "saving the transition embedding file..." << endl;
  #endif

  ofstream file(fileName);

  if(file) {
    for(auto it = _transitionEmbedding.begin(); it != _transitionEmbedding.end(); it ++) {
      file << it->first << "," << it->second.first << "," << it->second.second << endl;
    }
    file.close();
  }

}

/*----------------------------------------------------------------------------*/
void TransitionEmbedding::loadNeighbours(std::string fileName) {

  cout << "Loading the neighbours..." << endl;

  map<string,list<string>> vertices;

  DataFrame<string> df;
  df.loadFromCSV(fileName);

  // cout << df.toString() << endl;

  for(uint indRow = 0; indRow < df.nRows(); indRow ++) {

    string barcode1 = df.getData(indRow, 0);
    string barcode2 = df.getData(indRow, 1);


    /* insert the edge */
    map<string, list<string>>::iterator mit1;
    mit1 = vertices.find(barcode1);

    if(mit1 != vertices.end()) {
      mit1->second.push_back(barcode2);
    } else {
      auto mit1bis = vertices.insert(pair<string, list<string>>(barcode1, list<string>()));
      mit1bis.first->second.push_back(barcode2);
    }

    /* insert the reversed edge */
    // map<string, list<string>>::iterator mit2;
    // mit2 = vertices.find(barcode2);
    //
    // if(mit2 != vertices.end()) {
    //   mit2->second.push_back(barcode1);
    // } else {
    //   auto mit2bis = vertices.insert(pair<string, list<string>>(barcode2, list<string>()));
    //   mit2bis.first->second.push_back(barcode1);
    // }

  }

  /* final list of neighbours: add direct and un-direct neighbours */
  // map<string, vector<string>> neighbours;
  // for(auto& elt: vertices) {
  //
  //   pair<map<string, vector<string>>::iterator,bool> new_elt = neighbours.insert(pair<string, vector<string>>(elt.first, vector<string>()));
  //
  //   vector<string> neighbours_list = new_elt.first->second;
  //
  //   /* add the direct successors */
  //   neighbours_list.assign(elt.second.begin(), elt.second.end());
  //
  //   /* look for additional elements */
  //   for(auto& suc: elt.second) {
  //     auto& suc_suc_list = vertices.find(suc)->second;
  //     cout << suc_suc_list.size() << " ";
  //     neighbours_list.insert(neighbours_list.end(), suc_suc_list.begin(), suc_suc_list.end());
  //   }
  //
  //   sort(neighbours_list.begin(), neighbours_list.end());
  //   auto it = unique(neighbours_list.begin(), neighbours_list.end());
  //   neighbours_list.resize( std::distance(neighbours_list.begin(),it) );
  // }


  for(auto& elt: vertices) {

    pair<map<string, vector<string>>::iterator,bool> new_elt = _neighbours.insert(pair<string, vector<string>>(elt.first, vector<string>()));

    vector<string>& neighbours_list = new_elt.first->second;

    /* add the direct successors */
    neighbours_list.assign(elt.second.begin(), elt.second.end());

    /* make sure the key is not in the successors */
    auto it = remove(neighbours_list.begin(), neighbours_list.end(), elt.first);
    neighbours_list.resize( std::distance(neighbours_list.begin(),it) );

    sort(neighbours_list.begin(), neighbours_list.end());
    it = unique(neighbours_list.begin(), neighbours_list.end());
    neighbours_list.resize( std::distance(neighbours_list.begin(),it) );

  }

  // for(uint indRow = 0; indRow < /*df.nRows()*/100; indRow ++) {
  //   string barcode = _matrix.getRowLabel(indRow);
  //   auto& neighbours = _neighbours[barcode];
  //   cout << barcode << ":";
  //   for(auto& elt: neighbours) {
  //     cout << elt << "(" << _embedding[elt].first << "," << _embedding[elt].second << ") ";
  //   }
  //   cout << neighbours.size() << " " << _embedding[barcode].first << " " << _embedding[barcode].second << ", ";
  //   cout << endl;
  // }

}

/*----------------------------------------------------------------------------*/
void TransitionEmbedding::loadEmbeddingCoord(std::string fileName) {

  cout << "Loading the embedding..." << endl;

  DataFrame<double> embeddingFile;
  embeddingFile.loadFromCSV(fileName);

  for(uint indRow = 0; indRow < embeddingFile.nRows(); indRow ++) {

    string barcode = embeddingFile.getRowLabel(indRow);

    pair<double, double> coord;
    coord.first = embeddingFile.getData(indRow, 0);
    coord.second = embeddingFile.getData(indRow, 1);

    _embedding.insert(std::pair<string, std::pair<double,double>>(barcode, coord));

  }

  /* just check the data
  for(auto& elt: _embedding) {
    cout << elt.first << ": ";
    auto coord = elt.second;
    cout << coord.first << " " << coord.second << endl;
  }
  */

}
