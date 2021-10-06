/*!
 * \file LinearGraph.hpp
 * \author S. Buchet
 * \brief definition of LinearGraph class: linear correlation between the genes
 */

#ifndef LINEAR_GRAPH_HPP
#define LINEAR_GRAPH_HPP

#include "DataFrame.hpp"

class LinearGraph {

  public:

    /*!
     * \class Edge
     * \brief data structure to store the wieghted edges
     */
    class Edge {

      public:

        /*!
         * \brief constructor
         */
        Edge(uint pIndVar1, uint pIndVar2, double pweight) : indVar1(pIndVar1), indVar2(pIndVar2), weight(pweight) {

        }

        /*!
         * \brief default constructor
         */
        Edge(): indVar1(0), indVar2(0), weight(0.) {

        }

        uint indVar1;
        uint indVar2;
        double weight;

    };

  public:

    /*!
     * \brief constructor
     */
    LinearGraph(DataFrame<double>& dataset);

    /*!
     * \brief compute the graph
     * \param threshold correlation thresold to restrict the edges
     * \param fileName name of the output file
     */
    void computeGraph(double threshold, std::string fileName);

  private:

    /*!
     * \brief compute the pearson correlation coefficient between two variables
     * \return the correlation coefficient
     */
    double computeCor(uint indVar1, uint indVar2);

  private:

    DataFrame<double>& _dataset;

};

#endif
