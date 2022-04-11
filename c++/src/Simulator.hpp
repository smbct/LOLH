/*!
 * \file Simulator.hpp
 * \author S. Buchet
 * \brief definition of class Simulator
 */

#pragma once

#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "LogicProgram.hpp"

enum Semantics { ASYNCHRONOUS, SYNCHRONOUS, GENERALIZED};

/*!
 * \class Simulator
 * \brief Manage the execution of a logic program
 */
class Simulator {

  public:

    /*!
     * \brief constructor
     * \param program logic program to execute
     */
    Simulator(const LogicProgram& program);

    /*!
     * \brief generate all states of the program
     * \param states list of states generated
     */
    void generateAllStates(std::vector<State>& states);

    /*!
     * \brief generate all successors of a state given a semantics
     * \param state the origin state
     * \param semantics the semantics simulated
     * \param successors the successors of the state
     */
    void generateAllSuccessors(const State& state, Semantics semantics, std::vector<State>& successors);

  private: /* private methods */

    /*!
     * \brief generate all successors according to the synchronous semantics
     * \param state the origin state
     * \param successors the list of successors of the state
     */
    void generateSynchronousSuc(const State& state, std::vector<State>& successors);

    /*!
     * \brief generate all successors according to the asynchronous semantics
     * \param state the origin state
     * \param successors the list of successors of the state
     */
    void generateAsynchronousSuc(const State& state, std::vector<State>& successors);

  private:

    const LogicProgram& _program;

};

#endif
