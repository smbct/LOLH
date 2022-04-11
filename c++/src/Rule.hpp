/*!
 * \file Rule.hpp
 * \author S. Buchet
 * \brief definition of class Rule
 */

#pragma once

#ifndef RULE_HPP
#define RULE_HPP

#include <utility> /* std::pair */

#include <cstdlib> /* uint */
#include <vector>
#include <map>

typedef std::pair<uint, uint> Atom;

typedef std::vector<uint> State; /* alias for state of the system */

class Rule {

  public:

    /*!
     * \brief constructor
     * \param head head of the rule
     * \param body body of the rule
     */
    Rule(Atom head, std::vector<Atom>& body);

    /*!
     * \brief predicate true if the body matches a state
     * \param state the state to match
     * \return true iff there is a match
     */
    bool match(const State& state) const;

    /*!
     * \brief test if the rule subsumes another rule: i.e. the rule macthes more generally
     * \param rule the rule tested
     * \return true iff the rule subsumes rule
     */
    bool subsumes(const Rule& rule);

    /*!
     * \brief compute matching score (nb of error) of the body on a state
     * \param state the state tested
     * \return the matching score: number of atoms missmatching
     */
    uint matchingScore(const State& state) const;

    /*!
     * \brief return the length of the body (number of atoms)
     * \return length of the body
     */
    uint bodyLength() const;

    /*!
     * \brief return the head of the rule
     * \return Atom: head of the rule
     */
    Atom getHead() const;

    /*!
     * \brief return an atom of the body
     * \param index position in the body
     * \return Atom: atom in the body, at specified position
     */
    Atom getBody(uint index) const;

    /*!
     * \brief get a constant reference to the body
     * \return the body
     */
    const std::vector<Atom>& getBody() const;

  private:

    Atom _head; /* rule head */
    std::vector<Atom> _body; /* rule body */

};

#endif /* RULE_HPP */
