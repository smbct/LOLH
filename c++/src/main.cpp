#include <iostream>

#include "Utils.hpp"
#include "ClassificationQuality.hpp"
#include "NetworkInduction.hpp"

#include "Solver.hpp"

#include "TransitionEmbedding.hpp"

#include "LinearGraph.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
void regulCorrelationsExperiment() {

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  // string datasetFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/IMAGINE_normalised_discrete_adaptive_Tcells.csv";


  // string transitionsFileName = "dataset/IMAGINE/transitions_bidir_095.csv";
  // string transitionsFileName = "dataset/IMAGINE/transitions_25.csv";

  // string transitionsFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/transitions_10_Tcells.csv";
  // string transitionsFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/transitions_25_new_Tcells.csv";

  // string transitionsFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/transitions_25_new_new_Tcells.csv";

  string transitionsFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/transitions_dist4_Tcells.csv";



  // string result_path = "../../Learning/IMAGINE_dataset/CMSB/Tcells/";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  DataFrame<string> transitions;
  transitions.loadFromCSV(transitionsFileName);

  /* parameters */
  uint bodyLength = 20;
  double trRate = -1;
  uint predNeq = 1;
  uint delay_bis = 1; /* delay in neighborhood graph */

  /* compute a list of successors (indexes) from the transitions */
  NGraph successors;
  // DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay_bis, successors);

  // string output_filename = "regul_quality_25_pn1.txt";
  // output_filename = result_path+output_filename;

  // ClassificationQuality::computeRegulQuality(dataset, successors, bodyLength, trRate, predNeq, output_filename);
  // ClassificationQuality::computeCoexprQuality(dataset, bodyLength, "res/coexpr_quality_body.txt");

  ClassificationQuality::computeAtomCoexprQuality(dataset, "res/coexpr_quality_atom.txt");

  // NetworkInduction::computeNetwork(dataset, 0.3, "res/coexpr_network.txt");
  // NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_dist4_pn0_-1.txt");

  // uint bodyLength = 20;
  //
  // /* different trRate to tests */
  // vector<uint> delays = {1,2,3,4};
  //
  // // vector<double> trRateValues = {-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  // double trRate = -1;
  //
  // vector<uint> predNeqValues = {0, 1/*, 2*/};
  // // uint predNeq = 2;
  //
  // for(uint delay: delays) {
  //   NGraph successors;
  //   DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);
  //
  //   for(uint predNeq: predNeqValues) {
  //
  //     string output = "area_regul_IMAGINE_tr" + to_string(trRate) + "_b" + to_string(bodyLength) + "_d" + to_string(delay) + "_pred" + to_string(predNeq) + ".txt";
  //     output = "res/" + output;
  //
  //     ClassificationQuality::computeRegulQuality(dataset, successors, bodyLength, trRate, predNeq, result_path);
  //   }
  // }


  // for(double trRate: trRateValues) {
  //   for(uint predNeq: predNeqValues) {
  //
  //     // cout << "parameters: d=" << delay << " ; trRate= " << trRate << " ; predNeq=" << predNeq << endl;
  //
  //     ClassificationQuality::computeRegulQuality(dataset, successors, bodyLength, delay, trRate, predNeq);
  //
  //   }
  // }


}

/*----------------------------------------------------------------------------*/
void coexpressionNetworkExperiment() {

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  string result_path = "res/";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  /* parameters */
  double trRate = -1;
  uint predNeq = 1;
  uint delay = 1; /* delay in neighborhood graph */

  /* neighborhgood graph: 0.4 */
  {
    string transitionsFileName = "dataset/IMAGINE/transitions_bidir_04.csv";
    DataFrame<string> transitions;
    transitions.loadFromCSV(transitionsFileName);

    delay = 1;
    NGraph successors;
    DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);

    /* predNeq == 0 */
    predNeq = 0;
    trRate = 0.6;
    NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng04_d1_tr06_d1_pn0.txt");

    /* predNeq == 1 */
    predNeq = 1;
    trRate = -1.;
    NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng04_d1_tr-1_d1_pn1.txt");

    /* predNeq == 2 */
    predNeq = 2;
    NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng04_d1_tr-1_d1_pn2.txt");

  }

  /* neighborhood graph 0.95 */
  {

    string transitionsFileName = "dataset/IMAGINE/transitions_bidir_095.csv";
    DataFrame<string> transitions;
    transitions.loadFromCSV(transitionsFileName);

    /* delay == 2 */
    {

      delay = 2;
      NGraph successors;
      DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);

      /* predNeq == 0 */
      predNeq = 0;
      trRate = 0.6;
      NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng095_d2_tr06_d2_pn0.txt");

    }

    /* delay == 1 */
    {

      delay = 1;
      NGraph successors;
      DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);

      /* predNeq == 1 */
      predNeq = 1;
      trRate = -1.;
      NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng095_d1_tr-1_d1_pn1.txt");


      /* predNeq == 2 */
      predNeq = 2;
      trRate = 0.8;
      NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng095_d1_tr08_d1_pn2.txt");

    }

    /* delay == 4 */
    {

      delay = 4;
      NGraph successors;
      DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);

      /* predNeq == 2 comparison (red) */
      predNeq = 2;
      trRate = 0.7;
      NetworkInduction::computeRegulationNetwork(dataset, successors, trRate, predNeq, 0.3, result_path + "regul_network_ng095_d4_tr07_d1_pn2.txt");

    }

  }
}

/*----------------------------------------------------------------------------*/
void debugging() {

  /* load the dataset */
  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  /* load the transitions */
  string transitionsFileName = "dataset/IMAGINE/transitions_bidir_095.csv";
  DataFrame<string> transitions;
  transitions.loadFromCSV(transitionsFileName);


  /* parameters */
  double trRate = -1;
  uint predNeq = 1;
  uint delay = 1; /* delay in neighborhood graph */

  uint bodyLength = 20;

  /* create the graph of successors */
  NGraph successors;
  DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, delay, successors);

  uint nN = 0;
  for(auto& elt: successors) {
    nN += static_cast<uint>(elt.size());
  }
  // cout << "mean n neighbours: " << static_cast<double>(nN)/static_cast<double>(successors.size()) << endl;
  cout << "n transitions: " << nN << endl;

  /* target: CST3 0 */
  cout << "CST3 index: " << dataset.getColumnIndex("CST3") << endl;

  /* create the class vector for the classification instance */
  vector<bool> classVector(dataset.nRows(), false);

  Atom elt = Atom(12044, 0);

  uint nPos =  Instance::initRegulationClass(dataset, successors, elt.first, elt.second, trRate, predNeq, classVector);


  // vector<uint> features;
  // uint nPos =  Instance::initCoexprClass(dataset, elt.first, elt.second, features, classVector);

  cout << "Atom learned: " << dataset.getColLabel(12044) << endl;
  cout << "nPos: " << nPos << endl;

  // Instance instance(dataset, features, classVector);
  Instance instance(dataset, classVector); // adaptive discretization
  cout << instance.nPositives() << " vs " << instance.nNegatives() << endl;

  Solver solver(instance);

  vector<Body> solutions;
  vector<Score> points;
  vector<Weight> weights;

  solver.computeSupported(bodyLength, &points, &weights, &solutions);

  double area = Solver::relativeArea(points);

  cout << "relative area: " << area << endl;

  vector<uint> atomsNetwork;
  vector<double> atomsNetworkScores;

  cout << "Computation of the best atoms:" << endl;

  // solver.computekBestAtoms(nGenes, atomsNetwork, atomsNetworkScores);
  // selectionThreshold = 0.3;
  solver.computeBestAtomsThreshold(0.3, atomsNetwork, atomsNetworkScores);

  cout << "n atom selected: " + to_string(atomsNetwork.size()) << endl;
  cout << "n scores: " + to_string(atomsNetworkScores.size()) << endl;

}

/*----------------------------------------------------------------------------*/
void transitionEmbeddingTest() {

  // string matrixFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  string matrixFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/IMAGINE_normalised_discrete_adaptive_Tcells.csv";


  string embeddingFileName = "../../Learning/IMAGINE_dataset/dataset/embedding_coord.csv";

  /* test with a coexpression network */
  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/coexpression_network_0_14205_t03.txt";

  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/regul_network_25_pn1.txt";
  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/regul_network_10_pn1.txt";

  // string graphFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/regul_network_10_pn0_05.txt";
  // string graphFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/regul_network_25_new_pn0_05.txt";

  string graphFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/regul_network_dist4_pn0_-1_reduced.txt";


  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/pn1_random_labels.txt";
  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/regul_network_ng095_d1_tr-1_d1_pn1.txt";
  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/regul_network_ng095_d2_tr06_d2_pn0.txt";
  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/expe_regulation_network/regul_network_ng095_d4_tr07_d1_pn2.txt";

  // string graphFileName = "../../Learning/IMAGINE_dataset/Analysis/data/regul_network_25_pn1.txt";


  // string neighborhoodGraphFileName = "../../Learning/IMAGINE_dataset/dataset/neighborhood_graph_base.csv";
  // string neighborhoodGraphFileName = "../../Learning/IMAGINE_dataset/dataset/transitions_25.csv";
  // string neighborhoodGraphFileName = "../../Learning/IMAGINE_dataset/dataset/transitions_150.csv";

  string neighborhoodGraphFileName = "../../Learning/IMAGINE_dataset/CMSB/Tcells/neighbours_150_Tcells.csv";


  TransitionEmbedding trEmbedding(matrixFileName, embeddingFileName, graphFileName, neighborhoodGraphFileName);

}

/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[]) {

  Utils::getInstance().initRand(42);

  /* rapid example: use the gene clustering method to compute different logic rules corresponding to the gene clusters */

  DataFrame<uint> dataset;
  // dataset.loadFromCSV("../../Learning/IMAGINE_dataset/dataset/IMAGINE_normalised_discrete_adaptive_NK.csv");
  dataset.loadFromCSV("../../Learning/IMAGINE_dataset/dataset/IMAGINE_normalised_discrete_adaptive.csv");
  dataset.computeUniqueVal();

  cout << dataset.getRowLabel(0) << endl;

  /* load the cell types, to know which is NK */
  DataFrame<string> cellTypes;
  cellTypes.loadFromCSV("../../Learning/IMAGINE_dataset/dataset/cell_types.csv");

  /* list of cells that are not NK: true == not NK */
  vector<bool> negativeCells(dataset.nRows(), true);
  for(uint indRow = 0; indRow < cellTypes.nRows(); indRow ++) {
    if(cellTypes.getData(indRow, 0) == "NK") {
      negativeCells[dataset.getRowIndex(cellTypes.getRowLabel(indRow))] = false;
    }
  }


  NetworkInduction::computeNetwork(dataset, 0.3, "res/coexpr_network_NK_exclude.txt", &negativeCells);

  //NetworkInduction::computeNetwork(dataset, 0.5, "res/coexpr_network_global.txt");



  // LinearGraph graph(dataset);
  // graph.computeGraph(0.5, "res/linear_graph.txt");


  // regulCorrelationsExperiment();
  // transitionEmbeddingTest();


  // coexpressionNetworkExperiment();

  // debugging();

  // if(argc > 1) {
  //
  // cout << atoi(argv[1]) << endl;
  // } else {
  //
  // }


  return 0;

}
