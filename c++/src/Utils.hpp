/*!
 * \file Utils.hpp
 * \author S. Buchet
 * \brief utilities functions
 */

#ifndef UTILS_HPP
#define UTILS_HPP


#include <vector>

#include <algorithm>
#include <random>

/*!
 * \class Utils
 * \brief utilities functions
 */
class Utils {

  public:

    /*!
     * \brief get the instance
     * \return the unique instance
     */
    static Utils& getInstance();

    /*!
     * \brief copy constructor: deleted
     */
    Utils(Utils const&) = delete;

    /*!
     * \brief affectation constructor: deleted
     */
    void operator=(Utils const&) = delete;


  public:


    /*!
     * \brief init the random generator
     */
    void initRand();

    /*!
     * \brief init the random generator with a seed
     * \pram seed seed of the random generator
     */
    void initRand(unsigned int seed);

    /*!
     * \brief return a random number between 0 and 1
     * \return the number
     */
    double rand01();

    /*!
     * \brief return a random number in [min, max]
     * \return the number
     */
    double randInt(double min, double max);

    /*!
     * \brief return a random integer number in [min, max]
     * \return the number
     */
    unsigned int randUInt(unsigned int min, unsigned int max);

    /*!
     * \brief return and epsilon (very small value)
     */
    double epsilon();

    /*!
     * \brief generate a random selection
     * \param size size of the random selection
     * \param selection the random selection
     */
    void generateRandomSelection(unsigned int size, std::vector<unsigned int>& selection);

    /*!
     * \brief shuffle a vector of integer values
     * \param values the list of values to shuffle
     */
    void shuffle(std::vector<unsigned int>& values);

    /*!
     * \brief split a string line per line
     * \param str the string to split
     * \param output the list of string returned
     */
    void splitlines(std::string& str, std::vector<std::string>& output);

    /*!
     * \brief split a string at each delimiter
     * \param str the string to split
     * \param del the delimiter string
     * \param output the list of string returned
     */
    void split(std::string& str, std::string del, std::vector<std::string>& output);

    /*!
     * \brief merge sort
     * \param indexes indexes to sort
     * \param values values used to compare the elements
     */
    void mergeSort(std::vector<uint>& indexes, std::vector<int>& values);

    /*!
     * \brief compute the Manhattan distance between two vectors (same size), result \in [0,size]
     * \param left the left vector
     * \param right the right vector
     * \return distance
     */
    uint manDistance(std::vector<uint>& left, std::vector<uint>& right);

    /*!
     * \brief compute the dot product between two vectors (same size)
     * \param left the left vector
     * \param right the right vector
     * \return the dot product
     */
    double dotProduct(std::vector<uint>& left, std::vector<uint>& right);

    /*!
     * \brief compute the euclidean norm of a vector: sqrt(sum val^2)
     * \param vec the vector
     * \return the norm
     */
    double euclideanNorm(std::vector<uint>& vec);

    /*!
     * \brief compute the norm of a 2d vector
     * \param vec the vector to normalize
     * \return the norm
     */
    double norm(std::pair<double, double>& vec);

    /*!
     * \brief normalize a 2d vector
     * \param vec the vector to normalize
     */
    void normalize(std::pair<double, double>& vec);

  private:

    std::default_random_engine _randomGenerator;

  private:

    /*!
     * \brief constructor
     */
    Utils();

};

#endif
