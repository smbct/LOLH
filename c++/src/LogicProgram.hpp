/*!
 * \file LogicProgram.hpp
 * \author S. Buchet
 * \brief definition of class LogicProgram
 */

#pragma once

#ifndef LOGIC_PROGRAM_HPP
#define LOGIC_PROGRAM_HPP

#include <string>

#include "Rule.hpp"

/*!
 * \class LogicProgram
 * \brief data structure representing a logic program
 */
class LogicProgram {

  public:

    /*!
     * \brief constructor
     */
    LogicProgram();

    /*!
     * \brief load a logic program from a file
     * \param filename name of the file
     */
    void loadFromFile(std::string filename);

    /*!
     * \brief export the program to a file
     * \param fileName name of the file
     */
    void exportToFile(std::string fileName);

    /*!
     * \brief create an empty logic program
     * \param variables variable labels of the program
     * \param nValues number of values for each variable
     */
    void create(std::vector<std::string>& variables, uint nValues);

    /*!
     * \brief add a new rule to the body
     * \param head head of the rule
     * \param body body of the rule
     */
    void addRule(Atom& head, std::vector<Atom>& body);

    /*!
     * \brief add a new variable to a logic program
     * \param varName name of the variable
     * \param nValues number of logic values of the program
     */
    void addVariable(std::string varName, uint nValues);

    /*!
     * \brief initialize the program as a sub program of another program
     * \param program the original program
     * \param nVar number of variable to select
     */
    void createSubProgram(const LogicProgram& program, uint nVar);

    /*!
     * \brief number of variables in the program
     * \return the number of variable
     */
    uint nVariables() const;

    /*!
     * \brief return the name of a variable
     * \param index index of the variable
     */
    std::string varName(uint index) const;

    /*!
     * \brief return the index of a variable
     * |param varName name of the variable
     * \return index of the variable in the program
     */
    uint varIndex(std::string varName) const;

    /*!
     * \brief number of discrete values of a given variable
     * \param varIndex index of the variable
     * \return the number of values of this variable
     */
    uint nValues(uint varIndex) const;

    /*!
     * \brief get atoms that are conclusion of matching rule on a state
     * \param state a state of the program
     * \param list of values per variable, possible conclusions on the state
     */
    void getConclusions(const State& state, std::vector<std::vector<uint>>& atoms) const;

    /*!
     * \brief remove rules that are subsumed by other rules
     */
    void simplify();

    /*!
     * \brief display the program into a string
     * \return the string
     */
    std::string toString() const;

  private: /* private static methods */

    /*!
     * \brief extract an atom from a string, with form "label(value,T[-1])"
     * \param str string containing the atom
     * \return pair (index,value)
     */
    std::pair<uint,uint> extractAtom(std::string str);


  private:

    std::vector<std::string> _variables; /* variable labels */

    std::vector<uint> _nValues; /* number of values per variable */

    std::map<std::string, uint> _varIndex; /* dictionary between variable labels and indexes */

    std::vector<Rule> _rules; /* rules of the program */

    std::vector<std::vector<std::vector<uint>>> _ruleIndexes; /* index of rules classifies per variable index and value index */


};

#endif
