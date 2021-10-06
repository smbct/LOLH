/*!
 * \file Combinatorics.hpp
 * \author S. Buchet
 * \brief definition of class Combinatorics
 */

#pragma once

#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <vector>
#include <cstdlib>

/*!
 * \class Combinatorics
 * \brief utility functions for combinatorial computation
 */
class Combinatorics {

  public:

    /*!
     * \brief compute all subsets of k elements from a set of n elements
     * \param elements the set of elements
     * \param k the number of elements to pick
     * \param all possible subsets of size k
     */
    static void choose_kn(const std::vector<uint>& elements, uint k, std::vector<std::vector<uint>>& subsets);

    /*!
     * \brief compute all subsets of k elements from a set of element occurences
     * \param occurences occurence of each element
     * \param k the number of elements in the subsets
     * \param subsets the list of subsets
     */
    static void choose_kn_occurences(const std::vector<uint>& occurences, uint k, std::vector<std::vector<std::pair<uint,uint>>>& subsets);

    /*!
     * \brief generate all permutations of boolean vectors
     * \param size size of the vectors
     * \param permutations vector of the permutations, containing the result
     */
    static void generatePermutations(unsigned int size, std::vector<std::vector<bool>>& permutations);

    /*!
     * \brief generate a vector of indexes in specific range
     * \param length number of element: elements will be in 0..(length-1)
     * \param index the vector of index generated
     */
    static void generateRange(uint length, std::vector<uint>& index);

    /*!
     * \brief create a permutation of element given an index order
     * \param elements the elements to exchange
     * \param indexes the new order of the elements, given by indexes, can contain only a subset of indexes
     */
    template <class T>
    static void rearrenge(std::vector<T>& elements, const std::vector<uint>& indexes) {
      /* create a copy of elements for the exchange */
      std::vector<T> copy = elements;
      if(elements.size() > indexes.size()) {
        elements.resize(indexes.size());
      }
      for(uint ind = 0; ind < indexes.size(); ind ++) {
        elements[ind] = copy[indexes[ind]];
      }
    }

  private: /* private static methods */

    /*!
     * \brief recursive fuction to generate the permutations
     * \param index index in the current perm
     * \param perm the current permutation
     * \param permutations the list of all permutations
     */
    static void genereatePermutationsRec(unsigned int index, std::vector<bool>& perm, std::vector<std::vector<bool>>& permutations);

  private: /* private methods */

    /*!
     * \brief default constructor
     */
    Combinatorics();

};

#endif
