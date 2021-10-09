#include <iostream>

#include "Utils.hpp"
#include "ClassificationQuality.hpp"
#include "NetworkInduction.hpp"
#include "Parameters.hpp"

#include "Solver.hpp"



#include "LinearGraph.hpp"

using namespace std;


/*----------------------------------------------------------------------------*/
Parameters parametersExtraction(int argc, char* argv[]) {

  /* extraction of the commands */
  Parameters param;

  for(int index = 1; index < argc; index ++) {

    string argument = argv[index];

    if(argument == "-d") {
      param.debug = true;
    } else if(argument == "-v") {
      param.verbose = true;
    } else if(argument == "-r") {
        param.coexpression = false;
    } else if(argument == "-cq") {
      param.compute_network = false;
    } else if(argument == "-im") {
      if(index < argc-1) {
        param.input_matrix = argv[index+1];
        index ++;
      }
    } else if(argument == "-o") {
      if(index < argc-1) {
        param.output_file = argv[index+1];
        index ++;
      }
    } else if(argument == "-t") {
      if(index < argc-1) {
        param.threshold = stod(string(argv[index+1]));
        index ++;
      }
    } else if(argument == "-it") {
      if(index < argc-1) {
        param.input_transitions = argv[index+1];
        index ++;
      }
    } else if(argument == "-tr") {
      if(index < argc-1) {
        param.transition_rate = stod(string(argv[index+1]));
        index ++;
      }
    } else if(argument == "-pnq") {
      if(index < argc-1) {
        param.predecessor_neq = stoi(string(argv[index+1]));
        index ++;
      }
    } else if(argument == "-td") {
      if(index < argc-1) {
        param.transition_delay = stoi(string(argv[index+1]));
        index ++;
      }
    } else if(argument == "-h") {
      cout << "usage: " << endl;
      cout << "-d: execute the debug function" << endl;
      cout << "-v: display additional information on the computation" << endl;
      cout << "-r: compute the regulation network, if not, compute the coexpression network instead" << endl;
      cout << "-cq: only compute the instance quality, otherwise compute the network" << endl;
      cout << "-im: input matrix file (csv singl cell matrix)" << endl;
      cout << "-o: output file (gene network)" << endl;
      cout << "-t: selection threshold: 0. <= threshold <= 1." << endl;
      cout << "-it: input transitions file, for dynamic network only" << endl;
      cout << "-tr: transition rate between 0 and 1, or -1.: proportion of successors of a state required to verify the learned atom, to consider the state as positive, -1.: must be verified at least once" << endl;
      cout << "-pnq: 0 if the predecessor can be different or equal to the learned atom, 1 if it should be different, 2 if it should be equal" << endl;
      cout << "-td: transition delay (integer), > 1: length of the paths where the extremities are considered as a transition" << endl;
    } else {
      cout << "argument " << argument << " not recognised" << endl;
    }

  }

  return param;

}

/*----------------------------------------------------------------------------*/
void debug() {

  cout << "debugging function" << endl;

}

/*----------------------------------------------------------------------------*/
void computeRegulation(Parameters param) {

  /* load the discretized single cell matrix */
  DataFrame<uint> dataset;
  dataset.loadFromCSV(param.input_matrix);
  dataset.computeUniqueVal();

  /* load the transitions file */
  DataFrame<string> transitions;
  transitions.loadFromCSV(param.input_transitions);

  /* create the transition */
  NGraph successors;
  DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, param.transition_delay, successors);

  if(param.compute_network) {
    NetworkInduction::computeRegulationNetwork(dataset, successors, param.transition_rate, param.predecessor_neq, param.threshold, param.output_file);
  } else {
    ClassificationQuality::computeAtomRegulQuality(dataset, successors, param.transition_rate, param.predecessor_neq, param.output_file);
  }


}

/*----------------------------------------------------------------------------*/
void computeCoexpression(Parameters param) {

  DataFrame<uint> dataset;
  dataset.loadFromCSV(param.input_matrix);
  dataset.computeUniqueVal();

  if(param.compute_network) {
    NetworkInduction::computeNetwork(dataset, param.threshold, param.output_file);
  } else {
    ClassificationQuality::computeAtomCoexprQuality(dataset, param.output_file);
  }

}



/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[]) {

  cout << "*******************************************************************" << endl;
  cout << "        Computation of a gene network with LOLH algorithm          " << endl;
  cout << "*******************************************************************" << endl;

  Utils::getInstance().initRand(42);

  Parameters param = parametersExtraction(argc, argv);

  if(param.debug) {
    debug();
  } else {
    if(param.input_matrix.size() > 0 && param.output_file.size() > 0 && param.threshold >= 0) {
      if(param.coexpression) {
        computeCoexpression(param);
      } else {
        computeRegulation(param);
      }
    }
  }

  return 0;
}
