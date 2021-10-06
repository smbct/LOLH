/*!
 * \file Tests.cpp
 * \author S. Buchet
 * \brief implementatoin of tests functions
 */

#include "Tests.hpp"

/* time */
#include <chrono>
#include <unistd.h>

#include "Solver.hpp"
#include "LogicProgram.hpp"
#include "TransitionGraph.hpp"
#include "Induction.hpp"

#include "ClassificationQuality.hpp"
#include "NetworkInduction.hpp"

using namespace std;

/*----------------------------------------------------------------------------*/
void test1() {

  DataFrame<uint> dataframe;
  dataframe.loadFromCSV("dataset/mammalian/states.csv");

  // DataFrame<uint> df2;
  // df2.loadFromCSV("dataset/binary_logcpm.csv");

  // for(unsigned int indCol = 0; indCol < df2.nColumns(); indCol ++) {
  //   // cout << indCol << " : " << df2.nUnique(indCol) << endl;
  //   df2.nUnique(indCol);
  // }

  // cout << dataframe.toString() << endl;

  DataFrame<string> transitions_df;
  transitions_df.loadFromCSV("dataset/mammalian/transitions.csv");

  // vector<bool> a_1class;
  // Instance::createClassVector(dataframe, transitions_df, "c", 1, a_1class);

  /* create an instance */
  // Instance instance(dataframe, a_1class);

}


/*----------------------------------------------------------------------------*/
void test2() {

  auto start = chrono::steady_clock::now();
  DataFrame<uint> df;
  df.loadFromCSV("dataset/binary_logcpm.csv");
  auto stop = chrono::steady_clock::now();
  cout << "dataset loading execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;


  DataFrame<string> labels;
  labels.loadFromCSV("dataset/cells_labels_corrected_2.csv");

  uint df_nRows = df.nRows();
  vector<bool> classVector(df_nRows, false);

  uint labels_nRows = labels.nRows();

  /* check that the row index is indeed in the dataframe */
  // string feature = "MAIT1";
  string feature = "MAIT1";

  for(uint indRow = 0; indRow < labels_nRows; indRow ++) {
    string label = labels.getRowLabel(indRow);

    if(labels.getData(indRow, 0) == feature && df.containsRowIndex(label)) {
      classVector[df.getRowIndex(label)] = true;
    }
  }

  cout << "* creation of the instance" << endl;
  start = chrono::steady_clock::now();
  Instance instance(df, classVector, 2);
  stop = chrono::steady_clock::now();
  cout << "instance creation execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

  cout << instance.toString() << endl;

  // uint bodyLength = 10;
  vector<uint> solution(10, 0);

  cout << "* solving the instance" << endl;
  Solver solver(instance);

  /*--------------------------------------------------------------------------*/
  /* computation of all supported points */
  /*--------------------------------------------------------------------------*/

  vector<Body> solutions;
  vector<Score> points;
  vector<Weight> weights;

  uint bodyLength = 10;

  // solver.solveScalarization(bodyLength, Weight(5, -10), &points, &solutions);
  // solver.solveScalarization(bodyLength, Weight(10, -1), solution);
  // solver.solveScalarization(bodyLength, Weight(lbda1, lbda2), solution);

  cout << "Computation of the supported solutions: " << endl;

  start = chrono::steady_clock::now();
  solver.computeSupported(bodyLength, &points, &weights, &solutions);
  stop = chrono::steady_clock::now();
  cout << "computation supported time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

  cout << "Points: " << endl;
  for(auto& score : points) {
    cout << score << "  ";
  }
  cout << endl;

  cout << "Weights: " << endl;
  for(auto& weight : weights) {
    cout << weight << "  ";
  }
  cout << endl;

  cout << "Bodies: " << endl;
  for(auto& body : solutions) {
    cout << body << "  ";
  }
  cout << endl;

  cout << "Relative area: " << Solver::relativeArea(points) << endl;

  /* create and display an histogram of one rule */
  // Body& middle = solutions[solutions.size()/2];
  // Histogram histogram(instance, middle, Histogram::Type::POSITIVE_NEGATIVE);
  // cout << "histogram: " << endl;
  // cout << histogram << endl;
  //
  // cout << "positive covering 03: " << histogram.positiveCovered(3) << endl;
  // cout << "negative covering 03: " << histogram.negativeCovered(3) << endl;

  /*--------------------------------------------------------------------------*/
  /* computation of a target rule */
  /*--------------------------------------------------------------------------*/

  // Body solution2;
  // Score point;
  // Weight weight;
  //
  // uint bodyLength = 10;
  // uint targetCover = 20;
  // uint threshold = 2;
  //
  // cout << "body length: " << bodyLength << endl;
  // cout << "covering target: " << targetCover << " over " << instance.nPositives() << endl;
  // cout << "threshold: " << threshold << endl;
  //
  // start = chrono::steady_clock::now();
  // solver.computeTargetSolution(bodyLength, targetCover, threshold, point, weight, solution2);
  // stop = chrono::steady_clock::now();
  //
  // cout << "1st instance solving time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;
  //
  // cout << "Point result: " << point << endl;
  // cout << "Weight result: " << weight << endl;
  // cout << "Solution result: " << solution2 << endl;
  //
  // Histogram histogram(instance, solution2, Histogram::Type::  POSITIVE_NEGATIVE);
  // cout << "covering result: " << histogram.positiveCovered(threshold) << endl;
  // cout << "negative covering: " << histogram.negativeCovered(threshold) << endl;
  // // cout << histogram << endl;
  //
  // vector<uint> covered;
  // histogram.getPositiveCovered(threshold, covered);
  // cout << "positive samples covered: " << covered << endl;
  //
  // /* create a sub-instance */
  // start = chrono::steady_clock::now();
  // Instance subInstance(instance, covered, 2);
  // stop = chrono::steady_clock::now();
  //
  // cout << "sub instance creation time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;
  //
  //
  // Solver solver2(subInstance);
  // start = chrono::steady_clock::now();
  // solver2.computeTargetSolution(bodyLength, targetCover, threshold, point, weight, solution2);
  // stop = chrono::steady_clock::now();
  //
  // cout << "sub instance solving time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;
  //
  // cout << "sub solution: " << endl;
  // cout << solution2 << endl;
  // cout << point << endl;
  // cout << weight << endl;

  // Histogram histogram2(instance, solution2, Histogram::Type::POSITIVE_NEGATIVE);
  // cout << histogram2 << endl;

}

/*----------------------------------------------------------------------------*/
void test3() {

  LogicProgram program;

  program.loadFromFile("dataset/mammalian/mammalian.lp");
  // program.loadFromFile("dataset/mammalian/model_learned.lp");

  // program.loadFromFile("dataset/ERBB_G1-S/model_learned.lp");
  // program.loadFromFile("dataset/ERBB_G1-S/ERBB_G1-S.lp");

  LogicProgram subProgram;
  subProgram.createSubProgram(program, 10);

  cout << subProgram.toString() << endl;


  TransitionGraph graph(subProgram);
  graph.generateCompleteGraph(Semantics::ASYNCHRONOUS);
  // graph.exportCSV("dataset/ERBB_G1-S/original_sub10_complete");
  graph.exportCSV("dataset/mammalian/original_asynchronous_complete");

}

/*----------------------------------------------------------------------------*/
void test4() {

  // vector<uint> elements;
  // for(uint ind = 0; ind < 5; ind ++) {
  //   elements.push_back(ind);
  // }
  // vector<vector<uint>> subsets;
  // Combinatorics::choose_kn(elements, 4, subsets);
  // for(auto& subset: subsets) {
  //   for(uint ind = 0; ind < subset.size(); ind ++) {
  //     cout << subset[ind] << " ";
  //   }
  //   cout << endl;
  // }

  vector<uint> occurences = {10,5,4,6};
  // vector<uint> occurences = {5,3,5};
  vector<vector<pair<uint,uint>>> subsets2;

  Combinatorics::choose_kn_occurences(occurences, 10, subsets2);
  for(auto& subset : subsets2) {
    for(auto& elt : subset) {
      cout << "(" << elt.first << "," << elt.second << ")   ";
    }
    cout << endl;
  }

}

/*----------------------------------------------------------------------------*/
void test5() {

  auto start = chrono::steady_clock::now();
  DataFrame<uint> df;
  df.loadFromCSV("dataset/binary_logcpm.csv");
  auto stop = chrono::steady_clock::now();
  cout << "dataset loading execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;


  DataFrame<string> labels;
  labels.loadFromCSV("dataset/cells_labels_corrected_2.csv");

  uint df_nRows = df.nRows();
  vector<bool> classVector(df_nRows, false);

  uint labels_nRows = labels.nRows();

  /* check that the row index is indeed in the dataframe */
  // string feature = "MAIT1";
  string feature = "MAIT1";

  for(uint indRow = 0; indRow < labels_nRows; indRow ++) {
    string label = labels.getRowLabel(indRow);

    if(labels.getData(indRow, 0) == feature && df.containsRowIndex(label)) {
      classVector[df.getRowIndex(label)] = true;
    }
  }

  cout << "* creation of the instance" << endl;
  start = chrono::steady_clock::now();
  Instance instance(df, classVector, 2);
  stop = chrono::steady_clock::now();
  cout << "instance creation execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

  cout << instance.toString() << endl;

  cout << "* solving the instance" << endl;

  /*--------------------------------------------------------------------------*/
  /* inference of a variable's program */
  /*--------------------------------------------------------------------------*/

  uint bodyLength = 10;
  uint nbRules = 10;
  uint threshold = 2;

  cout << "body length: " << bodyLength << endl;
  cout << "nb rules goal: " << nbRules << endl;
  cout << "threshold: " << threshold << endl;

  // start = chrono::steady_clock::now();
  // Inference::inferRules(instance, bodyLength, nbRules, threshold);
  // stop = chrono::steady_clock::now();
  //
  // cout << "program inference time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;


}

/*----------------------------------------------------------------------------*/
void test6() {

  auto start = chrono::steady_clock::now();
  DataFrame<uint> df;
  df.loadFromCSV("dataset/binary_logcpm.csv");
  auto stop = chrono::steady_clock::now();
  cout << "dataset loading execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;


  DataFrame<string> labels;
  labels.loadFromCSV("dataset/cells_labels_corrected_2.csv");

  uint df_nRows = df.nRows();


  uint labels_nRows = labels.nRows();

  /* check that the row index is indeed in the dataframe */
  // string feature = "MAIT1";

  vector<string> features = {"MAIT0", "MAIT1", "cluster 7", "MAIT17a", "MAIT17b", "cyclingS", "cyclingG2M"};

  uint bodyLength = 10;

  for(auto feature : features) {

    cout << "current feature: " << feature << endl;

    /* classification column in the dataset */
    vector<bool> classVector(df_nRows, false);
    for(uint indRow = 0; indRow < labels_nRows; indRow ++) {
      string label = labels.getRowLabel(indRow);

      if(labels.getData(indRow, 0) == feature && df.containsRowIndex(label)) {
        classVector[df.getRowIndex(label)] = true;
      }
    }

    cout << "* creation of the instance" << endl;
    start = chrono::steady_clock::now();
    Instance instance(df, classVector, 2);
    stop = chrono::steady_clock::now();
    cout << "instance creation execution time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

    cout << instance.toString() << endl;

    cout << "* solving the instance" << endl;
    Solver solver(instance);

    vector<Body> solutions;
    vector<Score> points;
    vector<Weight> weights;

    cout << "Computation of the supported solutions: " << endl;

    start = chrono::steady_clock::now();
    solver.computeSupported(bodyLength, &points, &weights, &solutions);
    stop = chrono::steady_clock::now();
    cout << "computation supported time: " << chrono::duration_cast<chrono::milliseconds>(stop-start).count() << endl;

    cout << "Relative area of " << feature << ": " << Solver::relativeArea(points) << endl;
    cout << endl << endl << endl;
  }

}

/*----------------------------------------------------------------------------*/
void testParetoQuality() {

  DataFrame<uint> df;
  df.loadFromCSV("dataset/binary_logcpm_reduced.csv");

  DataFrame<string> labels;
  labels.loadFromCSV("dataset/cells_labels_corrected_2.csv");

  uint bodyLength = 10;

  string targetLabel = "MAIT1";

  /* classification column in the dataset */
  vector<bool> classVector(df.nRows(), false);
  Instance::initClusterClass(df, labels, targetLabel, classVector);

  vector<Body> solutions;
  vector<Score> points;
  vector<Weight> weights;

  Instance instance(df, classVector, 2);

  Solver solver(instance);

  solver.computeSupported(bodyLength, &points, &weights, &solutions);

  cout << "supported points: " << endl;
  for(auto& point: points) {
    cout << point << " ";
  }
  cout << endl << endl;

  cout << "Relative area of " << targetLabel << ": " << Solver::relativeArea(points) << endl;


}

/*----------------------------------------------------------------------------*/
void testLearnProgram() {

  string datasetFileName = "dataset/MAIT/binary_logcpm_reduced.csv";
  string transitionsFileName = "dataset/MAIT/transitions_n10_bidir.csv";

  Parameters param;
  param.bodyLength = 10; /* length of the rules body */
  param.graphDelay = 3; /* distance between cells of a transition in the graph */
  param.nRules = 7; /* number of rules to learn */
  param.posRate = 0.2; /* percentage of positive successor to consider the transition as positive */
  param.matchingThreshold = 2; /* min error to accept a sample on a rule */

  // Induction::createLFITProgram(datasetFileName, transitionsFileName, param, 0, 5);

}

/*----------------------------------------------------------------------------*/
void testQuality(uint indGene, uint val) {

  string datasetFileName = "dataset/MAIT/cells/binary_scanpy.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  cout << "dataset n columns: " << dataset.nColumns() << endl;
  cout << "atom target: " << dataset.getColLabel(indGene) << "_" << val << endl;

  /* create the class vector for the classification instance */
  vector<bool> classVector(dataset.nRows(), false);
  vector<uint> features;
  uint nPos =  Instance::initCoexprClass(dataset, indGene, val, features, classVector);

  cout << "features size: " << features.size() << endl;

  Instance instance(dataset, features, classVector, 2);

  uint bodyLength = 10;

  /* make sure the instance is well conditions -> non constant feature */
  if(nPos > 0 && nPos < classVector.size()) {

    Solver solver(instance);

    vector<Body> solutions;
    vector<Score> points;
    vector<Weight> weights;

    solver.computeSupported(bodyLength, &points, &weights, &solutions);

    cout << "n points: " << points.size() << endl;
    cout << "n solutions: " << solutions.size() << endl;

    double area = Solver::relativeArea(points);
    double maxDist, relativePos;
    Solver::maxMinDist(points, maxDist, relativePos);

    pair<double,double> posVarMean, negVarMean;
    Solver::meanVar(points, solutions, instance, posVarMean, negVarMean);

    cout << "area: " << area << endl;
    cout << "max dist: " << maxDist << endl;
    cout << "relative pos: " << relativePos << endl;

    cout << "pos mean mean: " << posVarMean.first << endl;
    cout << "pos var mean: " << posVarMean.second << endl;

    cout << "neg mean mean: " << negVarMean.first << endl;
    cout << "neg var mean: " << negVarMean.second << endl;

    Solver::meanVarMed(points, solutions, instance, posVarMean, negVarMean);
    cout << "with medians: " << endl;
    cout << "pos mean med: " << posVarMean.first << endl;
    cout << "pos var med: " << posVarMean.second << endl;

    cout << "neg mean med: " << negVarMean.first << endl;
    cout << "neg var med: " << negVarMean.second << endl;

    double posMean, negMean;
    Solver::medianMean(points, instance, posMean, negMean);
    cout << "median: of the means: " << posMean << " ; " << negMean << endl;

  }

}


/*----------------------------------------------------------------------------*/
void tuneTest(uint indGene, uint val) {

  string datasetFileName = "dataset/MAIT/cells/binary_scanpy.csv";
  // string transitionsFileName = "dataset/MAIT/transitions/transitions_bidir_scanpy.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  /* create the class vector for the classification instance */
  vector<bool> classVector(dataset.nRows(), false);
  vector<uint> features;
  uint nPos =  Instance::initCoexprClass(dataset, indGene, val, features, classVector);

  Instance instance(dataset, classVector);

  uint bodyLength = 10;

  /* make sure the instance is well conditions -> non constant feature */
  if(nPos > 0 && nPos < classVector.size()) {

    vector<bool> classVectorModified;
    Solver::tuneInstance(instance, bodyLength, 0.17, classVector, classVectorModified);

    /* create the modified instance */
    Instance instanceModified(dataset, classVectorModified);
    Solver solver(instanceModified);

    vector<Body> solutions;
    vector<Score> points;
    vector<Weight> weights;

    solver.computeSupported(bodyLength, &points, &weights, &solutions);

    cout << "Relative area (after tune): " << Solver::relativeArea(points) << endl;

  }

}



/*----------------------------------------------------------------------------*/
void computeCoexprArea() {
  // string datasetFileName = "dataset/MAIT/cells/binary_logcpm_reduced.csv";
  // string datasetFileName = "dataset/MAIT/cells/binary_scanpy_median.csv";
  // string datasetFileName = "dataset/MAIT/cells/adaptive_scanpy.csv";
  // string datasetFileName = "dataset/IMAGINE/IMAGINE_binary_mean.csv";
  // string datasetFileName = "dataset/IMAGINE/IMAGINE_binary_median.csv";
  // string datasetFileName = "dataset/IMAGINE/IMAGINE_scaled_3_discrete.csv";
  // string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_binary_0_reduced.csv";
  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  // string datasetFileName = "dataset/IMAGINE/IMAGINE_4_discrete.csv";
  // string datasetFileName = "dataset/IMAGINE/IMAGINE_log1p_4_discrete.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  uint bodyLength = 10;
  // ClassificationQuality::computeFastCoexprQuality(datasetFileName, bodyLength);
  ClassificationQuality::computeCoexprQuality(dataset, bodyLength, "res/corexpr_quality.txt");

}

/*----------------------------------------------------------------------------*/
void computeCoexprNetwork() {

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  // uint bodyLength = 50;

  // ClassificationQuality::computeNetworkOld(datasetFileName, bodyLength, 50, 0.4, min, max);
  NetworkInduction::computeNetwork(dataset, 0.3, "res/coexpr_network.txt");

}

/*----------------------------------------------------------------------------*/
void computeRegulationArea() {

  // string datasetFileName = "dataset/MAIT/cells/binary_logcpm_reduced.csv";
  // string datasetFileName = "dataset/MAIT/cells/binary_scanpy.csv";
  // string datasetFileName = "dataset/MAIT/cells/binary_scanpy_shuffled.csv";

  // string transitionsFileName = "dataset/MAIT/transitions/transitions_n2_invpt.csv";
  // string transitionsFileName = "dataset/MAIT/transitions/transitions_bidir_scanpy.csv";

  // string transitionsFileName = "dataset/MAIT/transitions/transitions_allpairs.csv"; /* default 15 neighbours */

  // string transitionsFileName = "dataset/MAIT/transitions/transitions_n10_bidir.csv";


  string orderingFileName = "dataset/MAIT/pseudotime.csv";

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";
  string transitionsFileName = "dataset/IMAGINE/transitions_bidir_04.csv";

  // string transitionsFileName = "dataset/MAIT/transitions_linear_rndpt.csv"; /* 5 neighbours */

  // string transitionsFileName = "dataset/MAIT/transitions_rndpt.csv"; /* random pseudotime */
  // string transitionsFileName = "dataset/MAIT/transitions_bidir_pt.csv"; /* bidirectional */
  // string transitionsFileName = "dataset/MAIT/transitions_cells_ng.csv";
  // string transitionsFileName = "dataset/MAIT/transitions_linear_tree.csv";

  // uint bodyLength = 20;
  // double trRate = 0.5;
  // bool predNeq = false;
  // ClassificationQuality::computeRegulQuality(datasetFileName, transitionsFileName, orderingFileName, bodyLength, trRate, predNeq);
}





/*----------------------------------------------------------------------------*/
void computeLFITProgram() {

  // string datasetFileName = "dataset/MAIT/cells/binary_logcpm_reduced.csv";
  // string transitionsFileName = "dataset/MAIT/transitions/transitions_n2_bidir.csv"; /* 5 neighbours */

  string datasetFileName = "dataset/MAIT/cells/binary_scanpy_selected.csv";
  string transitionsFileName = "dataset/MAIT/transitions/transitions_bidir_scanpy.csv"; /* 5 neighbours */

  string cellLabelFileName = "dataset/MAIT/cells_labels_corrected_2.csv";

  Parameters param;

  param.coExpression = true;
  param.bodyLength = /*4*/10; /* length of the rules body */
  param.graphDelay = 3; /* distance between cells of a transition in the graph */
  param.nRules = 20; /* number of rules to learn */
  param.posRate = 0.1; /* percentage of positive successor to consider the transition as positive */
  param.matchingThreshold = /*2*/0; /* min error to accept a sample on a rule */

  param.falsePositiveRate = /*0.1*/0.7; /* max false positive rate for a rule */

  param.clusterBodyLength = 10;

  Induction::createLFITProgram(datasetFileName, transitionsFileName, cellLabelFileName, param);

}

void test() {

  string datasetFileName = "dataset/MAIT/cells/binary_scanpy_selected.csv";
  string cellLabelFileName = "dataset/MAIT/cells_labels_corrected_2.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  cout << "gene: " << dataset.getColLabel(366) << endl;
  Parameters param;

  param.coExpression = true;
  param.bodyLength = /*4*/50; /* length of the rules body */
  param.graphDelay = 3; /* distance between cells of a transition in the graph */
  param.nRules = 10; /* number of rules to learn */
  param.posRate = 0.2; /* percentage of positive successor to consider the transition as positive */
  param.matchingThreshold = /*2*/0; /* min error to accept a sample on a rule */

  param.falsePositiveRate = 0.1; /* max false positive rate for a rule */

  vector<uint> features;

  Atom target(366,0);
  vector<bool> classVector(dataset.nRows(), false);

  uint nPos = Instance::initCoexprClass(dataset, target.first, target.second, features, classVector);

  /* make sure the instance is not ill conditioned */
  if(nPos > 0 && nPos < classVector.size()) {

    Instance instance(dataset, features, classVector, 2);
    vector<Body> solutions;
    Induction::ruleInduction(instance, param, solutions, false);
  }

}

/*----------------------------------------------------------------------------*/
uint getGeneIndex(string gene_name) {
  string datasetFileName = "dataset/MAIT/cells/binary_scanpy.csv";
  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  return dataset.getColumnIndex(gene_name);
}

/*----------------------------------------------------------------------------*/
void computeGeneInfluence(uint indGene, uint val) {

  string datasetFileName = "dataset/MAIT/cells/binary_scanpy.csv";
  string transitionsFileName = "dataset/MAIT/transitions/transitions_bidir_scanpy.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  cout << "dataset n columns: " << dataset.nColumns() << endl;

  cout << "atom target: " << dataset.getColLabel(indGene) << "_" << val << endl;

  /* create the class vector for the classification instance */
  vector<bool> classVector(dataset.nRows(), false);
  vector<uint> features;
  uint nPos =  Instance::initCoexprClass(dataset, indGene, val, features, classVector);

  cout << "features size: " << features.size() << endl;

  cout << nPos << endl;

  Instance instance(dataset, features, classVector, 2);

  // uint bodyLength = 10;

  // vector<uint> geneInfluence(100, 0);

  /* make sure the instance is well conditions -> non constant feature */
  // if(nPos > 0 && nPos < classVector.size()) {
  //
  //   Solver solver(instance);
  //   // should be replaced by compute best k atoms
  //   solver.computeInfluentialVar(bodyLength, geneInfluence);
  //
  // }

}

/*----------------------------------------------------------------------------*/
void regulCorrelationsExperiment2(uint delay) {

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";
  string transitionsFileName = "dataset/IMAGINE/transitions_bidir_095.csv";

  /* old experiment with an ordering from pseudotime */

  string orderingFileName("");

  cout << delay << endl;

  // DataFrame<string> ordering;
  // ordering.loadFromCSV(orderingFileName);

  // DataFrame<string> labelsDataframe;
  // labelsDataframe.loadFromCSV("dataset/MAIT/cells_labels_corrected_2.csv");
  //
  // set<string> labels;
  // map<string, string> cellLabel;
  // for(uint ind = 0; ind < labelsDataframe.nRows(); ind ++) {
  //   cellLabel.insert(pair<string, string>(labelsDataframe.getRowLabel(ind), labelsDataframe.getData(ind, 0)));
  //   labels.insert(labelsDataframe.getData(ind, 0));
  // }

  /*for(auto& elt: cellLabel) {
    cout << elt.first << " : " << elt.second << endl;
  }*/

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);

  DataFrame<string> transitions;
  transitions.loadFromCSV(transitionsFileName);

  /* parameters */
  uint bodyLength = 20;
  double trRate = 0.4;
  uint predNeq = 0;
  uint delay_bis = 4; /* delay in neighborhood graph */

  /* compute a list of successors (indexes) from the transitions */
  NGraph successors;
  // DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, successors);
  DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, /*ordering,*/ delay_bis, successors);

  string output_filename = "res/regul_quality.txt";

  ClassificationQuality::computeRegulQuality(dataset, successors, bodyLength, trRate, predNeq, output_filename);

  // NGraph successors;
  // DataFrame<uint>::createNeighbourhoodGraph(dataset, transitions, /*ordering,*/ delay, successors);
  //
  // uint bodyLength = 20;
  //
  // /* different trRate to tests */
  // vector<double> trRateValues = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
  // vector<bool> predNeqValues = {false, true};
  //
  // for(double trRate: trRateValues) {
  //   for(bool predNeq: predNeqValues) {
  //
  //     // cout << "parameters: d=" << delay << " ; trRate= " << trRate << " ; predNeq=" << predNeq << endl;
  //
  //     ClassificationQuality::computeRegulQuality(dataset, successors, bodyLength, delay, trRate, predNeq);
  //
  //   }
  // }


}

/*----------------------------------------------------------------------------*/
void testRelativeAreaAtoms() {

  string datasetFileName = "dataset/IMAGINE/IMAGINE_normalised_discrete_adaptive.csv";

  DataFrame<uint> dataset;
  dataset.loadFromCSV(datasetFileName);
  dataset.computeUniqueVal();

  uint indVar = dataset.getColumnIndex("POU2F2");
  uint val = 1;

  vector<bool> classVector(dataset.nRows(), false);
  vector<uint> features;
  uint nPos = Instance::initCoexprClass(dataset, indVar, val, features, classVector);

  if(nPos > 0 && nPos < classVector.size()) {

    Instance instance(dataset, features, classVector);

    Solver solver(instance);

    vector<uint> atoms;
    Combinatorics::generateRange(instance.nAtoms(), atoms);

    vector<uint> nondominatedAtoms;

    solver.computeNonDominatedAtoms(atoms, nondominatedAtoms);

    double relativeArea = solver.computeRelativeAreaAtoms(nondominatedAtoms);

    cout << "relative area: " << relativeArea << endl;

  }

}
