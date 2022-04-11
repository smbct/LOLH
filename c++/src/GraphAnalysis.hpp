/*!
 * \file GraphAnalysis.hpp
 * \author S. Buchet
 * \brief definition of class GraphAnalysis
 */

#pragma once

#ifndef GRAPH_ANALYSIS_HPP
#define GRAPH_ANALYSIS_HPP

#include "DataFrame.hpp"

/*!
 * \class GraphAnalysis
 * \brief functions to compute indicators on the neiborhood graph
 */
class GraphAnalysis {

  public: /* public static functions */

    static void test();

    static void dijkstra(int origin, NGraph& graph,  std::vector<int>& distances);

    static void connectedComponents(NGraph& graph, std::vector<std::vector<int>>& comp);

  private: /* private static functions */

  private: /* private constructor (static class) */
    GraphAnalysis();

};

#endif
