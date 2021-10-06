/*!
 * \file Histogram.hpp
 * \author S. Buchet
 * \brief definition of class Histogram
 */

#pragma once

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "Instance.hpp"

/*!
 * \class Histogram
 * \brief data structure representing the distribution of a matching score over positive and negative samples
 */
class Histogram {

  public: /* public enumeration */

    enum Type {POSITIVE, NEGATIVE, POSITIVE_NEGATIVE};


  public: /* pulic lethods */

    /*!
     * \brief default constructor
     * \param instance instance of the classification problem
     * \param solution solution of the histogram
     * \param type histogram type: pôsitive samples only, negative only or positive and negative
     */
    Histogram(const Instance& instance, const Body& solution, Type type);

    /*!
     * \brief return the histogram type
     */
    Type type() const;

    /*!
     * \brief compute the number of positive samples covered according to a threshold: i.e. error(solution,sample) <= threshold
     * \param theshold the covering threshold
     * \pre type() != Type::NEGATIVE
     */
    uint positiveCovered(uint threshold);


    /*!
     * \brief get a list of indexes of positive samples covered by the rule
     * \param threshold covering threshold
     * \param indexes indexes of positive samples covered
     */
    void getPositiveCovered(uint threshold, std::vector<uint>& indexes);

    /*!
     * \brief compute the number of negative samples covered according to a threshold: i.e. error(solution,sample) <= threshold
     * \param theshold the covering threshold
     * \pre type() != Type::POSITIVE
     */
    uint negativeCovered(uint threshold);

    /*!
     * \brief get the positive samples having a given score on the rule
     * \param score the score of the samples
     * \return a constant vector containing the sample indexes
     */
    const std::vector<uint>& getPosSampleScore(uint score);

    /*!
     * \brief get the negative samples having a given score on the rule
     * \param score the score of the samples
     * \return a constant vector containing the sample indexes
     */
    const std::vector<uint>& getNegSampleScore(uint score);

    /*!
     * \brief first function to display the histogram in the terminal
     * \param str the stream to display the histogram
     * \return the stream
     */
    std::ostream& display1(std::ostream& str) const;

    /*!
     * \brief first function to display the histogram in the terminal
     * \param str the stream to display the histogram
     * \return the stream
     */
    std::ostream& display2(std::ostream& str) const;

  private: /* private methods */

    /*!
     * \brief build the positive histogram
     */
    void buildPositive();

    /*!
     * \brief build the negative histogram
     */
    void buildNegative();

  private:

    friend std::ostream& operator<<(std::ostream& str, const Histogram& histogram);

  private: /* private attributes */

    const Instance& _instance;
    const Body& _solution;
    Type _type;

    std::vector<std::vector<uint>> _positives; /* positive samples corresponding to each possible score */
    std::vector<std::vector<uint>> _negatives; /* negatives samples corresponding to each possible score */

};

/*!
 * \brief ostream operator << overload for Histogram
 * \param str a stream
 * \param histogram the histogram to display
 * \return the output stream
 */
std::ostream& operator<<(std::ostream& str, const Histogram& tuple);

#endif
