/*!
 * \file Instance.hpp
 * \author S. Buchet
 * \brief definition of Instance class
 */

#pragma once

#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include <string>
#include <vector>

#include <iostream>

#include "DataFrame.hpp"

typedef std::vector<uint> Sample; /* alias to describe a line of the dataset */
typedef std::pair<uint, uint> Atom; /* alias for atoms (values of variables) */
typedef std::pair<uint, uint> Score; /* alias for a bi-objective score */
typedef std::vector<uint> Body; /* alias for the body of a rule: list of atom indexes */

/*!
 * \class Instance
 * \brief represent an classification instance
 */
class Instance {

  public:

    /*!
     * \brief constructor
     * \param dataset training dataset of the model
     * \param classVector vector of class [true,false] for each sample
     */
    Instance(const DataFrame<uint>& dataset, const std::vector<bool>& classVector);

    /*!
     * \brief alternative constructor
     * \param dataset training dataset of the model
     * \param classVector vector of class [true,false] for each sample
     * \param nval number of discrete values used per variables
     */
    Instance(const DataFrame<uint>& dataset, const std::vector<bool>& classVector, uint nval);

    /*!
     * \brief alternative constructor
     * \param dataset training dataset of the model
     * \param features list of dataset features to predict from (column index)
     * \param classVector vector of class [true,false] for each sample
     * \param nval number of discrete values used per variables
     */
    Instance(const DataFrame<uint>& dataset, const std::vector<uint>& features, const std::vector<bool>& classVector, uint nval);

    /*!
     * \brief alternative constructor
     * \param dataset training dataset of the model
     * \param features list of dataset features to predict from (column index)
     * \param classVector vector of class [true,false] for each sample
     */
    Instance(const DataFrame<uint>& dataset, const std::vector<uint>& features, const std::vector<bool>& classVector);

    /*!
     * \brief constructor: create an instance from an existing one
     * \param instance the existing instance
     * \param excluded: list of positive samples to exclude in the new instance
     * \param nval number of values per variables
     */
    Instance(const Instance& instance, const std::vector<uint>& excluded, uint nVal);


    /*!
     * \brief destructor
     */
    ~Instance();

    /*!
     * \brief get the index of an atom
     * \param atom the atom
     * \return the index of the atom
     */
    uint getAtomIndex(Atom atom) const;

    /*!
     * \brief get an atom given its index
     * \param index index of the atom
     * \return the corresponding atom
     */
    Atom getAtom(uint index) const;

    /*!
     * \brief get the atom with variable name, given its index
     * \param index index of the atom
     * \return the corresponding atom
     */
    std::pair<std::string, uint> getAtomLabel(uint index) const;

    /*!
     * \brief return the score of an atom as a pair
     * \param index index of the atom
     * \return positive score in result.first, negative score in result.second
     */
    Score getAtomScore(uint index) const;

    /*!
     * \brief return the score of an atom as a pair
     * \param atom the atom
     * \return positive score in result.first, negative score in result.second
     */
    Score getAtomScore(Atom atom) const;

    /*!
     * \brief return the positive score of an atom given its index
     * \param index index of the atom
     * \return positive score of the atom
     */
    uint getAtomPosScore(uint index) const;

    /*!
     * \brief return the positive score of an atom
     * \param atom the atom
     * \return positive score of the atom
     */
    uint getAtomPosScore(Atom& atom) const;

    /*!
     * \brief return the negative score of an atom given its index
     * \param index index of the atom
     * \return negative score of the atom
     */
    uint getAtomNegScore(uint index) const;

    /*!
     * \brief return the negative score of an atom
     * \param atom the atom
     * \return negative score of the atom
     */
     uint getAtomNegScore(Atom& atom) const;

     /*!
      * \brief compute the matching score of a body on a sample
      * \param body body of a rule
      * \param indSample index of the sample
      */
     uint matchingScore(const Body& body, uint indSample) const;

     /*!
      * \brief compute the matching score of a body on a posiive sample
      * \param body body of a rule
      * \param indPositive index of the positive sample
      */
     uint positiveMatchingScore(const Body& body, uint indPositive) const;

     /*!
      * \brief compute the matching score of a body on a negative sample
      * \param body body of a rule
      * \param indNegative index of the negative sample
      */
     uint negativeMatchingScore(const Body& body, uint indNegative) const;

    /*!
     * \brief get a sample from its index
     * \param ind index of the sample
     * \return the sample
     */
    const Sample& getSample(uint index) const;

    /*!
     * \brief get a positive sample from its index
     * \param ind index of the positive sample
     * \return the positive sample
     */
    const Sample& getPositive(uint index) const;

    /*!
     * \brief get a negative sample from its index
     * \param ind index of the negative sample
     * \return the negative sample
     */
    const Sample& getNegative(uint index) const;

    /*!
     * \brief return the index in the class vector of a positive sample
     * \param indPos the index in the set of positive samples
     * \return the index in the class vector
     */
    uint getPosClassIndex(uint indPos) const;

    /*!
     * \brief return the index in the class vector of a negative sample
     * \param indNeg the index in the set of negative samples
     * \return the index in the class vector
     */
    uint getNegClassIndex(uint indNeg) const;

    /*!
     * \brief get the number of samples
     * \return number of samples
     */
    uint nSamples() const;

    /*!
     * \brief get the number of positive samples
     * \return number of positive samples
     */
    uint nPositives() const;

    /*!
     * \brief get the number of negative samples
     * \return number of negative samples
     */
    uint nNegatives() const;

    /*!
     * \brief get the number of variable of the dataset
     * \return the number of variable
     */
    uint nVariables() const;

    /*!
     * \brief number of prediction features used
     * \return the number of prediction features
     */
    uint nFeatures() const;

    /*!
     * \brief get the features used for the prediction
     * \return a constant reference to the list of feature indexes
     */
    const std::vector<uint>& getFeatures() const;

    /*!
     * \brief get the number of discrete values for a given variable index
     * \param varInd the index of the variable
     * \return the number of value of the variable
     */
    uint nValues(uint varInd) const;

    /*!
     * \brief return the number of atoms in the instance: \sum v:variable nValues(v)
     * \return number of atoms
     */
    uint nAtoms() const;

    /*!
     * \brief return the dataset on which the instance is based
     * \return a constant reference toward the dataset
     */
    const DataFrame<uint>& getDataset() const;

    /*!
     * \brief put instance info in a string
     * \return a string
     */
    std::string toString() const;

  private: /* private methods */

    /*!
     * \brief initialize the instace: positive/negative samples, scores
     * \param classVector class [positive, negative] for each sample
     */
    void initialize(const std::vector<bool>& classVector);

    /* faster alternative ? */
    /* nvalue: number of values per variables */
    void initialize(const std::vector<bool>& classVector, uint nval);

    /*
     * \brief initialize the instance from a sub instance, by excluding some samples
     * \param instance the original instance
     * \param excluded a list of positive samples that are excluded
     * \param nValues number of values per variables
     */
    void initializeSubInstance(const Instance& instance, const std::vector<uint>& excluded, uint nValues);

  public: /* public static methods */

    /*!
     * \brief initialize a class vector representing positive/negative samples of a random classification instance
     * \param proportion proportion of positive samples in the dataset
     * \param classVector vector of the class computed: positive or negative
     * \return the number of positive samples of the class
     */
    static uint initRandomClass(double proportion, std::vector<bool>& classVector);

    /*!
     * \brief initialize a class vector representing a cell cluster label instance
     * \param dataset the single cell dataset
     * \param labels the set of labels for all cells
     * \param cellLabel the label to classify
     * \param classVector vector of the class computed: positive or negative
     * \return the number of positive samples in the class
     */
    static uint initClusterClass(DataFrame<uint>& dataset, DataFrame<std::string>& labels, std::string cellLabel, std::vector<bool>& classVector);

    /*!
     * \brief initialize a class vector to classify the co-expression of a gene
     * \param dataset the single cell dataset
     * \param geneIndex index of the gene to classify
     * \param value the value of the gene to classify
     * \param features the list of features for the instance (remove the predicted feature)
     * \return the number of positive samples in the instance
     * \param classVector vector of the class computed: positive or negative
     */
    static uint initCoexprClass(DataFrame<uint>& dataset, uint geneIndex, uint value, std::vector<uint>& features, std::vector<bool>& classVector);

    /*!
     * \brief initialize a class vector to classify the regulation of a gene
     * \param dataset the single cell dataset
     * \param graph the graph of successors
     * \param geneIndex index of the gene to classify
     * \param value the value of the gene to classify
     * \param proportion proportion of positive successors to be considered as positive samples, if <= 0.5: at least one
     * \param predNeq 1 will only select predecessor for which geneIndex value is not value, 2: geneIndex is value
     * \param classVector vector of the class computed: positive or negative
     * \return the number of positive samples in the instance
     */
    static uint initRegulationClass(DataFrame<uint>& dataset, NGraph& graph, uint geneIndex, uint value, double proportion, uint predNeq, std::vector<bool>& classVector);

  private:

    unsigned int _nvar;
    std::vector<uint> _nval; /* number of values for each feature of the dataset */

    std::vector<Atom> _atoms; /* list of atoms */
    std::vector<std::vector<uint>> _atomIndexes; /* index of each atom in the list */

    std::vector<uint> _positives; /* positive sample indexes in the dataset */
    std::vector<uint> _negatives; /* negative sample indexes in the dataset */

    std::vector<std::vector<uint>> _posScore; /* score for each atom in positive samples */
    std::vector<std::vector<uint>> _negScore; /* score for each atom in negative samples */

    const DataFrame<uint>& _dataset;

    std::vector<uint> _features; /* features used for the prediction in the dataset: column index corresponding to the instance variable index  */
};

/*!
 * \brief ostream operator << overload for Instance
 */
std::ostream& operator<<(std::ostream& str, const Instance& instance);

/*!
 * \brief ostream operator << overload for std::pair<uint,uint>
 * \param str a stream
 * \param tuple the tuple to display
 * \return the output stream
 */
std::ostream& operator<<(std::ostream& str, const std::pair<uint,uint>& tuple);

/*!
 * \brief ostream operator << overload for Body (std::vector<uint>)
 * \param str a stream
 * \param body the body to display
 * \return the output stream
 */
std::ostream& operator<<(std::ostream& str, const Body& body);

/*!
 * \brief return the body as a string
 * \param the body
 * \return the corresponding string
 */
std::string bodyToString(const Body& body);


#endif
