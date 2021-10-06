/*!
 * \file Solver.hpp
 * \author S. Buchet
 * \brief definition of class Solver
 */

#pragma once

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <vector>

#include "Instance.hpp"

#include "Combinatorics.hpp"

enum Objective { POSITIVE, NEGATIVE };

typedef std::pair<int,int> Weight; /* alias for a pair of weight: scalarization of the bi-objective function */

/*!
 * \class Solver
 * \brief solver algorithm that computes supported solutions of un-constrained bi-objective linear program
 */
class Solver {

  public:

    /*!
     * \brief constructor
     */
    Solver(const Instance& instance);

    /*!
     * \brief solve the instance considering one objective
     * \param bodyLength length of the corresponding body
     * \param obj objective solved
     * \param solution one solution associated to the objective value, if not null
     * \return the optimal score obtained
     */
    Score solveMono(uint bodyLength, Objective obj, Body* solution = nullptr);

    /*!
     * \brief solve a scalarization of the two objectives
     * \param bodyLength length of the body
     * \param weight a pair of weight: obj = pos*score * weight.firsts + neg_score * weight.second
     * \param points list of points computed
     * \param solutions list of bodies computes
     */
    void solveScalarization(uint bodyLength, Weight weight, std::vector<Score>* points, std::vector<Body>* solutions = nullptr);

    /*!
     * \brief compute the list of supported solutions of the instance
     * \param bodyLength length of the body
     * \param points list of points computed
     * \param solutions a list of solutions associated to the points
     */
    void computeSupported(uint bodyLength, std::vector<Score>* points, std::vector<Weight>* weights = nullptr, std::vector<Body>* solutions = nullptr);

    /*!
     * \brief compute a solution covering at least a given number of positive samples
     * \param bodyLength length of the body solution
     * \param targetCover number of positive samples covered as target
     * \param targetPoint the best point found
     * \param targetWeigth corresponding weight
     * \param targetSolution corresponding solution
     * \return true iff the target has been achieved (if it is not possible, the function returns false)
     */
    bool computeTargetSolution(uint bodyLength, uint targetCover, uint threshold, Score& targetPoint, Weight& targetWeight, Body& targetSolution);

    /*!
     * \brief compute a solution covering at least a given number of positive samples: slow version which compute all the supported solutions
     * \param bodyLength length of the body solution
     * \param targetCover number of positive samples covered as target
     * \param targetPoint the best point found
     * \param targetWeigth corresponding weight
     * \param targetSolution corresponding solution
     * \return true iff the target has been achieved (if it is not possible, the function returns false)
     */
    bool computeTargetSolutionSlow(uint bodyLength, uint targetCover, uint threshold, Score& targetPoint, Weight& targetWeight, Body& targetSolution);

    /*!
     * \brief computation of the k best atoms for this instance
     * \param nAtom number of atoms to compute
     * \param selectedAtoms list of atoms returned
     * \param atomScores scores of the selected atoms
     */
    void computekBestAtoms(uint nAtoms, std::vector<uint>& selectedAtoms, std::vector<double>& atomScores);

    /*!
     * \brief computation of the best atoms for this instance, based on a threshold
     * \param threshold the threshold, atoms are selected if their score is >= threshold
     * \param selectedAtoms list of atoms returned
     * \param atomScores scores of the selected atoms
     */
    void computeBestAtomsThreshold(double threshold, std::vector<uint>& selectedAtoms, std::vector<double>& atomScores);


    /*!
     * \brief compute a list of non dominated atoms from a list of give atoms
     * \param atoms atoms that are given
     * \param nondominatedAtoms atoms from the initial list that are not dominated
     */
    void computeNonDominatedAtoms(const std::vector<uint>& atoms, std::vector<uint>& nondominatedAtoms);

    /*!
     * \brief compute the relative area formed by a list of ordered atoms (potentially not dominated)
     * \param atoms list of atom indexes to compute the area from
     * \return the relative area
     */
    double computeRelativeAreaAtoms(const std::vector<uint>& atoms);

  public: /* public static methods */

    /*!
     * \brief compute the relative domination area of supported solutions
     * \param supported set of ordered supported points: must contain optimum on each individual objective
     * \return the relative area: % percentage of dominated space
     */
    static double relativeArea(/*const Instance& instance, uint bodyLength, */const std::vector<Score>& supported);

    /*!
     * \brief compute the max distance of a point to the diagonal of the pareto front
     * \param supported set of ordered supported points: must contain optimum on each individual objective
     * \param maxDist the maximum distance found
     * \param relativePos the relative position of the projection on the diagonal of the farthest point, in [0,1]
     */
    static void maxMinDist(const std::vector<Score>& supported, double& maxDist, double& relativePos);

    /*!
     * \brief compute the distances of all points from the diagonal
     * \param the supported points
     * \param distances the distances that are computes
     */
    static void computeDistances(const std::vector<Score>& supported, std::vector<double>& distances);

    /*!
     * \brief compute the mean of the mean histogram error and the mean histogram variance on pos and neg samples
     * \param supported supported solutions in the objective space
     * \param bodies bodies of the supported solutions
     * \param instance instance of the problem
     * \param posVarMean mean and variance on positive samples
     * \param negVarMean mean and variance on negative samples
     */
    static void meanVar(const std::vector<Score>& supported, const std::vector<Body>& bodies, const Instance& instance, std::pair<double, double>& posVarMean, std::pair<double, double>& negVarMean);

    /*!
     * \brief compute the median of the mean histogram error and the mean histogram variance on pos and neg samples
     * \param supported supported solutions in the objective space
     * \param bodies bodies of the supported solutions
     * \param instance instance of the problem
     * \param posVarMean mean and variance on positive samples
     * \param negVarMean mean and variance on negative samples
     */
    static void meanVarMed(const std::vector<Score>& supported, const std::vector<Body>& bodies, const Instance& instance, std::pair<double, double>& posMed, std::pair<double, double>& negMed);

    /*!
     * \brief compute the median of the positive and the negative mean of the points
     * \param supported the list of the supported points
     * \param instance the instance of the problem
     * \param medPosMean the median of the positive mean [output]
     * \param medNegMean the median of the negative mea [output]
     */
    static void medianMean(const std::vector<Score>& supported, const Instance& instance, double& medPosMean, double& medNegMean);

    /*!
     * \brief modify an instance based on the pareto front to improve its correlation score
     * \param instance the instance to modify
     * \param bodyLnegth length of the bodues (solutions)
     * \param ratio of samples to transfer in the instance (poitive <-> negative samples)
     * \param classVector the class vector of the original instance
     * \param classVetcorModified the class vector of the new instance
     */
    static void tuneInstance(const Instance& instance, uint bodyLength, double ratio, std::vector<bool>& classVector, std::vector<bool>& classVectorModified);

  private: /* private methods */

    /*!
     * \brief select atoms that have the same multi-objective score
     * \param indexes list of ordered indexes
     * \param bodyLength length of the body
     * \return range of indexes containing bodyLength-1 such that all indexes in its range have the same score
     */
    std::pair<uint,uint> selectEquivalentScores(const std::vector<uint>& indexes, uint bodyLength);

    /*!
     * \brief select atoms that have the same scalarized score
     * \param indexes list of ordered indexes
     * \param bodyLength length of the body
     * \param score the score of each atom in indexes
     * \return range of indexes containing bodyLength-1 such that all indexes in its range have the same score
     */
    std::pair<uint,uint> selectEquivalentScores(const std::vector<uint>& indexes, uint bodyLength, const std::vector<int>& score);


    /*!
     * \brief create a list of atom indexes by selecting the best atom for each variable
     * \param indexes list of indexes
     * \param obj the objective considered
     */
    void createIndexList(std::vector<uint>& indexes, Objective obj);

    /*!
     * \brief create a list of atom indexes by selecting the best atom for each variable
     * \param indexes list of indexes
     * \param score the score used to select the best index
     */
    void createIndexList(std::vector<uint>& indexes, const std::vector<int>& score);

    /*!
     * \brief create subsets of elements with the same bi-objective score
     * \param indexes the list of elements
     * \param equiv the range of elements to consider
     * \param result an association list between the list of equivalent elements and their value
     */
    void createAssociationList(std::vector<uint>& indexes, std::pair<uint,uint> equiv, std::vector<std::vector<uint>>& result);

    /*!
     * \brief filter supported points according to their bounds and order them and remove duplicates
     * \param left the left bound point
     * \param right the right bound point
     * \param points list of points computed
     * \param solutions list of associated solutions, if not null
     */
    void filterPoints(Score left, Score right, std::vector<Score>* points, std::vector<Body>* solutions = nullptr);

    /*!
     * \brief rank al the atoms for this instance, based on their relative distance
     * \param selectedAtoms list of atoms returned
     * \param atomScores scores of the selected atoms
     */
    void computeBestAtoms(std::vector<uint>& sortedAtoms, std::vector<double>& sortedScores);

  private:

    const Instance& _instance; /* instance of the classification problem */

};

/*!
 * \brief ostream operator << overload for std::pair<int,int> (Weight)
 * \param str a stream
 * \param tuple the tuple to display
 * \return the output stream
 */
std::ostream& operator<<(std::ostream& str, const std::pair<int,int>& tuple);

#endif
