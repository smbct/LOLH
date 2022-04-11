/*!
 * \file DataFrame.hpp
 * \author S. Buchet
 * \brief definition of class DataFrame
 */

#pragma once

#ifndef DATA_FRAME_HPP
#define DATA_FRAME_HPP

#include <string>

#include <sstream>
#include <vector>
#include <map>

typedef std::vector<std::vector<uint>> NGraph;

/*!
 * \class DataFrame
 * \brief data structure representing a dataset
 */
template <class T> class DataFrame {

  public: /* public methods */

    /*!
     * \brief constructor
     */
    DataFrame();

    /*!
     * \brief load a dataframe from a csv file
     * \param filename name of the file
     */
    void loadFromCSV(std::string filename);

    /*!
     * \brief display the DataFrame into a string
     * \return the string
     */
    std::string toString();

    /*!
     * \brief get the label of a column
     * \param index index of the column
     * \return label of the column
     */
    std::string getColLabel(uint index) const;

    /*!
     * \brief get the column labels as a list of string
     * \param labels list of labels returned
     */
    void getColumnLabels(std::vector<std::string>& labels);

    /*!
     * \brief return true iff the dataset contains label in the row indexes
     * \param label the index tested
     * \return true iif label is in the indexes
     */
    bool containsRowIndex(std::string label);

    /*!
     * \brief get the index of a column
     * \param label column label
     * \return index of the column
     */
    uint getColumnIndex(std::string label);

    /*!
     * \brief get the label of a row
     * \param index index of the row
     * \return label of the row
     */
    std::string getRowLabel(uint index);

    /*!
     * \brief get the index of a row
     * \param label row label
     * \return index of the row
     */
    uint getRowIndex(std::string label);

    /*!
     * \brief return the number of rows of the dataset
     * \return number of rows
     */
    uint nRows() const;

    /*!
     * \brief return the number of columns of the dataset
     * \return number of columns
     */
    uint nColumns() const;

    /*!
     * \brief get an element of the dataset at specific position
     * \param rowIndex row index of the element
     * \param colIndex columns index of the element
     * \return the element
     */
    T& getData(uint rowIndex, uint colIndex);

    /*!
     * \brief get a sample (row) in the dataset
     * \param rowIndex index of the row
     * \return constant reference to the row
     */
    const std::vector<T>& getRow(uint rowIndex) const;

    /*!
     * \brief get a row in the dataset from its label
     * \param rowLabel lael of the row
     * \return constant reference to the row
     */
    const std::vector<T>& getRow(std::string rowLabel) const;

    /*!
     * \brief number of unique values in a column
     * \param colInd index of the column
     * \return number of unique values
     */
    uint nUnique(uint colInd) const;

    /*!
     * \brief return a map between unique values and occurences for a given column
     * \param colInd index of the column
     * \return a map between unique values and counts
     */
    std::map<T, uint> uniqueCount(uint colInd) const;

    /*!
     * \brief return a map between unique values and occurences for a given column in specified rows
     * \param colInd index of the column
     * \param rows index of the rows
     * \return a map between unique values and counts in corresponding rows
     */
    std::map<T, uint> uniqueCount(uint colInd, const std::vector<uint>& rows) const;

    /*!
     * \brief count unique values and occurences for a given column in specified rows
     * \param colInd index of the column
     * \param rows index of the rows
     * \param list of values where the coutns are returned
     */
    void uniqueCount(uint colInd, const std::vector<uint>& rows, std::vector<uint>& counts) const;

    /*!
     * \brief count unique values and occurences for all columns, for two different classes
     * \param classVector the two classes corresponding to the two counts
     * \param posCounts the positive counts where classVector is true
     * \param negCounts the negative counts, where classVector is false
     */
    void uniqueCount(const std::vector<bool>& classVector, std::vector<std::vector<uint>>& posCounts, std::vector<std::vector<uint>>& negCounts) const;

    /*!
     * \brief count unique values and occurences for some columns, for two different classes
     * \param colIndexes indexes of the columns from which to count the occurences
     * \param classVector the two classes corresponding to the two counts
     * \param posCounts the positive counts where classVector is true
     * \param negCounts the negative counts, where classVector is false
     */
    void uniqueCount(const std::vector<uint>& colIndexes, const std::vector<bool>& classVector, std::vector<std::vector<uint>>& posCounts, std::vector<std::vector<uint>>& negCounts) const;

    /*!
     * \brief count value occurences in a subset of rows of the dataset
     * \param rows the rows from which the occurences are computed
     * \param counts the occurences for all variables
     */
    void uniqueCount(const std::vector<uint>& rows, std::vector<std::vector<uint>>& counts) const;

    void uniqueCount(std::vector<std::vector<uint>>& counts) const;

    /*!
     * \brief pre-compute the number of unique values for each variable
     * \pre the variables are integers or categorical
     */
    void computeUniqueVal();

  public: /* public static methods */

    /*!
     * \brief compute a neighborhood graph (list of successors) based on an observation dataframe and a transition dataframe
     * \param dataframe the dataframe containing all data
     * \param transitions the dataframe containing the transitions: pair of data indexes
     * \graph the graph created
     */
    static void createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<std::string>& transitions, NGraph& graph);

    /*!
     * \brief compute a neighborhood graph (list of successors) with a given delay, based on an observation dataframe and a transition dataframe
     * \param dataframe the dataframe containing all data
     * \param transitions the dataframe containing the transitions: pair of data indexes
     * \param delay the delay used in the original graph
     * \graph the graph created
     */
    static void createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<std::string>& transitions, uint delay, NGraph& graph);

    /*!
     * \brief compute a neighborhood graph (list of successors) with a given delay and pseudotime ordering, based on an observation dataframe and a transition dataframe
     * \param dataframe the dataframe containing all data
     * \param transitions the dataframe containing the transitions: pair of data indexes
     * \param ordering pseudotime ordering over all datapoint
     * \param delay the delay used in the original graph
     * \graph the graph created
     */
    static void createNeighbourhoodGraph(DataFrame<uint>& dataframe, DataFrame<std::string>& transitions, DataFrame<std::string>& ordering, uint delay, NGraph& graph);


  protected:

    /*!
     * \brief read a sample from the file
     * \param stream the stream to the data
     * \param the sample that stores the data
     */
    void readSample(std::istringstream& stream, std::vector<T>& sample);


  protected: /* protected attributes */

    std::vector<std::vector<T>> _data;

    /* list of columns and rows names */
    std::vector<std::string> _columns;
    std::vector<std::string> _rows;

    /* dictionary between rox/columns names and their indexes */
    std::map<std::string, uint> _colIndexes;
    std::map<std::string, uint> _rowIndexes;

    /* number of values for each variables, when var is integer or categorical */
    std::vector<uint> _nVal;
    uint _maxVal; // maximum possible value in the dataset


};

#endif /* DATA_FRAME_HPP */
