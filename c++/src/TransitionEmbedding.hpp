/*!
 * \file TransitionEmbedding.hpp
 * \author S. Buchet
 * \brief definition of a class TransitionEmbedding
 */

#pragma once

#ifndef TRANSITION_EMBEDDING_HPP
#define TRANSITION_EMBEDDING_HPP

#include <string>
#include <map>

#include "DataFrame.hpp"

#include "Instance.hpp" /* typedef Score */


/*!
 * \class RegulatoryGraph
 * \brief data structure to load and store a regulatory graph infered from the transitions
 */
class RegulatoryGraph {

  public: /* public methods */

    typedef std::pair<std::string, uint> Atom; /* alias for atoms: warning, not the same definition in Instance */

    /*
     * \brief default constructor
     */
    RegulatoryGraph();

    /*
     * \brief load the graph from a file
     * \param fileNmae the name of the file
     */
    void loadFromFile(std::string fileName);

    /*
     * \brief compute the partial successor of a cell regarding the regulatory graph
     * \param matrix the discrete gene expression matrix
     * \param cell the cell barcode
     * \param threshold a score threshold to accept the atom as predecessor
     * \param successors the list of (partial) successors computed
     */
    void computePartialSucessor(DataFrame<uint>& matrix, std::string& cell, double threshold, std::vector<Atom>& successors);

  private: /* private methods */

    /*
     * \brief get the index of an atom: the atom is created if necessary
     * \return the index
     */
    uint getAtomIndex(Atom atom);

  private: /* private attributes */

    std::vector<Atom> _atoms;
    std::map<Atom, uint> _atomIndexes;

    std::vector<std::vector<uint>> _successors; /* atom indexes of the successors */
    std::vector<std::vector<Score>> _discreteScores; /* for each edge, positive and negative score of the edge */
    std::vector<std::vector<double>> _relativeScores;

    std::vector<std::vector<std::pair<uint,uint>>> _predecessors; /* for each atom, access the edge of the predecessor atom */

    std::vector<Score> _proportions; /* proportions of positive/negative examples for each atom */

};

/*!
 * \class TransitionEmbedding
 * \brief representation of the transitions in a 2d embedding, using the discrete matrix and a neghborhood graph
 */
class TransitionEmbedding {

  public: /* public methods */

    /*!
     * \brief constructor
     * \param matrixFileName name of the file containing the discrete expression
     * \param embeddingFileName name of the file containing the coordinates of the cells in the embedding
     * \param graphFileName name of the file containing the regulatory graph
     * \param neighborhoodGraphFileame name of the file containing the neighborhood graph
     */
    TransitionEmbedding(std::string matrixFileName, std::string embeddingFileName, std::string graphFileName, std::string neighborhoodGraphFileName);

    /*!
     * \brief compute a probability distribution of transition from one cell to another one
     */
    void computeTransitionProbabilities();

    /*!
     * \brief compute the coordinates of the transition vector in the embedding from the cordinates and the transition probabilities
     */
    void computeTransitionEmbedding();

    /*!
     * \brief save the transition coordinates of the embedding into a file
     * \param fileName name of the file
     */
    void saveTransitionEmbedding(std::string fileName);

  private: /* private methods */

    /*!
     * \brief load the embedding coordinates of the cells
     * \param fileName file containing the coordinates
     */
    void loadEmbeddingCoord(std::string fileName);

    /*!
     * \brief load the file containing the neighborhood graph
     * \param fileName name of the file containing the neighbours
     */
    void loadNeighbours(std::string fileName);


  private: /* private attributes */

    std::map<std::string, std::pair<double,double>> _embedding; /* 2d coordinates of the cells */

    DataFrame<uint> _matrix; /* the discrete single cell expression matrix */

    std::map<std::string, std::vector<std::string>> _neighbours; /* list of neighbours for each cell barcode */

    RegulatoryGraph _regulatoryGraph; /* the graph infered from the transitions, representing dynamical relation between the genes */

    std::map<std::string, std::vector<double>> _transitionProbabilites; /* a probability distribution of the transitions from one cell to its neighbours */

    std::map<std::string, std::pair<double, double>> _transitionEmbedding; /* transition vector in the embedding space for each cell */
};

#endif
