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
    } else if(argument == "-me") {
      if(index < argc-1) {
        param.max_edges = stoi(string(argv[index+1]));
        index ++;
      }
    } else if(argument == "-h") {
      cout << "usage: " << endl;
      cout << "-d: execute the debug function" << endl;
      cout << "-v: display additional information on the computation" << endl;
      cout << "-r: compute the regulation network, if not, compute the coexpression network instead" << endl;
      cout << "-cq: only compute the instance quality, otherwise compute the network" << endl;
      cout << "-im: input matrix file (csv single cell matrix)" << endl;
      cout << "-o: output file (gene network)" << endl;
      cout << "-t: selection threshold: 0. <= threshold <= 1." << endl;
      cout << "-it: input transitions file, for dynamic network only" << endl;
      cout << "-tr: transition rate between 0 and 1, or -1.: proportion of successors of a state required to verify the learned atom, to consider the state as positive, -1.: must be verified at least once" << endl;
      cout << "-pnq: 0 if the predecessor can be different or equal to the learned atom, 1 if it should be different, 2 if it should be equal" << endl;
      cout << "-td: transition delay (integer), > 1: length of the paths where the extremities are considered as a transition" << endl;
      cout << "-me: maximum number of edges for transition learning, default value = 0, meaning there is no limitation" << endl;
    } else {
      cout << "argument " << argument << " not recognised" << endl;
    }

  }

  return param;

}

/*----------------------------------------------------------------------------*/
void debug() {

  /* NK classification test */

  cout << "debugging function" << endl;

  string matrix_filename = "../dataset/Imagine/discrete_matrix.csv";
  string celltypes_filename = "../dataset/Imagine/cell_types.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(matrix_filename);
  dataset.computeUniqueVal();

  DataFrame<string> labels;
  labels.loadFromCSV(celltypes_filename);

  vector<bool> classVector(dataset.nRows(), false);
  vector<uint> features;

  uint nPos = Instance::initClusterClass(dataset, labels, "NK", classVector);

  Instance instance(dataset, classVector);

  Solver solver(instance);

  vector<uint> atomIndexes;
  vector<double> atomScores;

  double selectionThreshold = 0.55;
  solver.computeBestAtomsThreshold(selectionThreshold, atomIndexes, atomScores);

  for(uint ind = 0; ind < atomIndexes.size(); ind ++) {
    Atom atom = instance.getAtom(atomIndexes[ind]);
    cout << dataset.getColLabel(atom.first) << " " << atom.second << endl;
  }

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
    if(param.threshold >= 0) {
      NetworkInduction::computeRegulationNetwork(dataset, successors, param.transition_rate, param.predecessor_neq, param.threshold, param.output_file, param.max_edges);
    } else {
      /* error */
    }
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
    if(param.threshold >= 0) {
      NetworkInduction::computeNetwork(dataset, param.threshold, param.output_file, nullptr, param.max_edges);
    } else {
      /* error */
    }
  } else {
    ClassificationQuality::computeAtomCoexprQuality(dataset, param.output_file);
  }

}

/*----------------------------------------------------------------------------*/
void LIGER_expe_dynamics() {

  string matrix_file = "../dataset/Imagine/discrete_matrix.csv";
  string transitions_file = "../dataset/Imagine/transitions.csv";
  string output_file = "";

  // vector<float> tau_values = {-1, 0.2, 0.4, 0.6, 0.8};
  // vector<uint> rho_values = {0,1,2};
  // vector<uint> delta_values = {1,2};

  vector<float> tau_values = {0.4};
  vector<uint> rho_values = {0};
  vector<uint> delta_values = {2};

  Parameters param;
  param.compute_network = false;
  param.coexpression = false;
  param.input_matrix = matrix_file;
  param.input_transitions = transitions_file;


  /* load the discretized single cell matrix */
  DataFrame<uint> dataset;
  dataset.loadFromCSV(param.input_matrix);
  dataset.computeUniqueVal();

  /* load the transitions file */
  DataFrame<string> transitions;
  transitions.loadFromCSV(param.input_transitions);


  for(float tau: tau_values) {
    for(uint rho: rho_values) {
      for(uint delta: delta_values) {

        param.transition_rate = tau;
        param.predecessor_neq = rho;
        param.transition_delay = delta;
        param.output_file = "dynamics_quality_t_" + to_string(tau) + "_r_" + to_string(rho) + "_d_" + to_string(delta) + ".txt";

        // cout << tau << ", " << rho << ", " << delta << endl;
        // cout << param.output_file << endl;

        /* create the transition */
        NGraph successors;
        DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, param.transition_delay, successors);

        ClassificationQuality::computeAtomRegulQuality(dataset, successors, param.transition_rate, param.predecessor_neq, param.output_file);

      }
    }
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
    if(param.input_matrix.size() > 0 && param.output_file.size() > 0) {
      if(param.coexpression) {
        computeCoexpression(param);
      } else {
        computeRegulation(param);
      }
    }
  }

  // LIGER_expe_dynamics();

  return 0;
}
