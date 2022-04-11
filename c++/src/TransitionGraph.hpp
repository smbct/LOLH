/*!
 * \file TransitionGraph.hpp
 * \author S. Buchet
 * \brief definition of class TransitionGraph
 */

#pragma once

#ifndef TRANSITION_GRAPH_HPP
#define TRANSITION_GRAPH_HPP

#include "LogicProgram.hpp"
#include "Simulator.hpp"

/*!
 * \class TransitionGraph
 * \brief graph containing the states and transitions (edges) generated from a model
 */
class TransitionGraph {

  public:

    /*!
     * \brief constructor
     * \param program program represented by the graoh
     */
    TransitionGraph(const LogicProgram& program);

    /*!
     * \brief generate the complete graph
     * \param semantics the semantics used
     */
    void generateCompleteGraph(Semantics semantics);

    /*!
     * \brief export the graph in two csv files, a vertex file an edge file
     * \param fileName name of the base file
     */
    void exportCSV(std::string fileName);

  private:

    const LogicProgram& _program;
    std::vector<State> _states;
    std::map<State, uint> _stateIndexes;
    std::vector<std::vector<uint>> _transitions;

};



#endif
