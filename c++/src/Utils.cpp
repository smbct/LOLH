/*!
 * \file Utils.cpp
 * \brief utilities functions
 */

#include "Utils.hpp"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;



/*----------------------------------------------------------------------------*/
Utils& Utils::getInstance() {

  static Utils instance;

  return instance;

}

/*----------------------------------------------------------------------------*/
Utils::Utils():
_randomGenerator(42)
{

}

/*----------------------------------------------------------------------------*/
void Utils::initRand() {
  srand(static_cast<unsigned int>(time(NULL)));
}

/*----------------------------------------------------------------------------*/
void Utils::initRand(unsigned int seed) {
  srand(seed);
}

/*----------------------------------------------------------------------------*/
double Utils::rand01() {
  return static_cast<double>(rand()%10001)/10000.;
}

/*----------------------------------------------------------------------------*/
double Utils::randInt(double min, double max) {
  return min + (max-min)*rand01();
}

/*----------------------------------------------------------------------------*/
unsigned int Utils::randUInt(unsigned int min, unsigned int max) {
  double randomVal = rand01();
  if(1.-randomVal < Utils::epsilon()) {
    randomVal -= Utils::epsilon();
  }
  double interval = 1./static_cast<double>(max-min+1); // cautious here if max-min+1 is very big
  return min + static_cast<unsigned int>(floor(randomVal/interval));
}

/*----------------------------------------------------------------------------*/
double Utils::epsilon() {
  return 10e-5;
}

/*----------------------------------------------------------------------------*/
void Utils::generateRandomSelection(unsigned int size, std::vector<unsigned int>& selection) {

  vector<unsigned int> indexes(size);
  for(unsigned int i = 0; i < size; i ++) {
    indexes.at(i) = i;
  }

  for(unsigned int i = 0; i < selection.size(); i ++) {

    unsigned int randInd = Utils::randUInt(0, static_cast<unsigned int>(indexes.size()));
    selection.at(i) = indexes.at(randInd);
    remove(indexes.begin(), indexes.end(), selection.at(i));

  }

}

/*----------------------------------------------------------------------------*/
void Utils::shuffle(std::vector<unsigned int>& values) {
  std::shuffle(values.begin(), values.end(), _randomGenerator);
}

/*----------------------------------------------------------------------------*/
void Utils::splitlines(string& str, vector<string>& output) {
  split(str, "\n", output);
}

/*----------------------------------------------------------------------------*/
void Utils::split(string& str, string del, vector<string>& output) {

  bool stop = false;

  std::size_t last = string::npos;

  uint del_size = static_cast<uint>(del.size());

  while(!stop) {

    std::size_t found;

    if(last == string::npos) {
      found = str.find(del);
    } else {
      found = str.find(del, last+del_size);
    }

    if (found != string::npos) {

      if(last == string::npos) {
        output.push_back(str.substr(0, found));
      } else {
        output.push_back(str.substr(last+del_size, found-last-del_size));
      }

      last = found;

    } else {

      /* one last sub-string */
      if(last == string::npos) {
        output.push_back(str.substr(0, str.size()));
      } else {
        output.push_back(str.substr(last+del_size, str.size()-last));
      }

      stop = true;
    }

  }

}

/*----------------------------------------------------------------------------*/
uint Utils::manDistance(vector<uint>& left, vector<uint>& right) {
  uint dist = 0;
  for(uint ind = 0; ind < left.size(); ind ++) {
    if(left[ind] != right[ind]) {
      dist += 1;
    }
  }
  return dist;
}

/*----------------------------------------------------------------------------*/
double Utils::dotProduct(vector<uint>& left, vector<uint>& right) {
  double res = 0.;
  for(uint ind = 0; ind < left.size(); ind ++) {
    res += left[ind]*right[ind];
  }
  return res;
}

/*----------------------------------------------------------------------------*/
double Utils::euclideanNorm(vector<uint>& vec) {
  double res = 0.;
  for(uint ind = 0; ind < vec.size(); ind ++) {
    res += vec[ind]*vec[ind];
  }
  res = sqrt(res);
  return res;
}

/*----------------------------------------------------------------------------*/
double Utils::norm(pair<double, double>& vec) {
  return sqrt(vec.first*vec.first+vec.second*vec.second);
}

/*----------------------------------------------------------------------------*/
void Utils::normalize(pair<double, double>& vec) {
  double norm = this->norm(vec);
  vec.first /= norm;
  vec.second /= norm;
}
