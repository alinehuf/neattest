//
//  main.c
//  NeatTest_0.1
//
//  Created by AlineHUF on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//
//  largely based on Mat Buckland 2002

/* TODO : 
 * genome.c / getCompatibilityScore() : 
 *                   Stanley propose to have longest=1 if both genomes are small
 *                   (fewer than 20 genes) => compatibility threshold must be
 *                   altered accordingly ...
 *
 * NOT IMPLEMENTED :
 * Stanley : « There was a 75% chance that an inherited gene was disabled if it 
 * was disabled in either parent. In each generation, 25% of offspring resulted 
 * from mutation without crossover. The interspecies mating rate was 0.001. »
 *
 * DONE  
 * population / epoch() :
 *                   Stanley recommends to clone species leader only if the 
 *                   species has more than five networks, in Neat1.1 he used the
 *                   futur number of individual instead of the current one...
 *                   => does it improve something ? 
 * phenotype.c / sigmoid() :
 *                   Stanley used a constant activation to have a sigmoid slope
 *                   near to linear : 1/(1+(exp(-4.924273 * activesum)))
 *                   => is that the activation variable, subject to mutations,
 *                      causes a disturbance ?
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>   // init srand()
#include <unistd.h> // getpid()
#ifdef GRAPHVIZ
  #include <gvc.h>  // graphViz
#endif
#include "data.h"
#include "global.h"

/* GLOBALES */
int Verbose = 1;
int UltraVerbose = 0;
int GlobalLog = 1;
int Dump = 1;
int GraphPs = 0;

FILE * GlobalLogFile = NULL;

//const char * Log_dir = "log/";
const char * Log_dir = "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/log/";

/*******************************************************************************
 * prototypes
 ******************************************************************************/

void testOneGame(char * gamePath, char * paramsPath);
char * getName(char * path);
void createGlobalLogFile(const char * dir);
void createDumpFile(const char * dir, const char * gameName, sPopulation * pop,
                    testBase * data, sParams * params);

void dumpInputsFactsAndGenome(FILE * out, testBase * data, sGenome * maxGen);
void populationTest(); // XOR test (population, speciation...)
void genomeToGraphVizDotFile(sGenome * gen, testBase * data, sPopulation * pop,
                             const char * dir, char * gameName);
void genomeToGraphViz(sGenome * gen, testBase * data, sPopulation * pop,
                      const char * dir, char * gameName);

/*******************************************************************************
 * MAIN
 ******************************************************************************/

int main(int argc, char * argv[]) {
  srand((unsigned int) time(NULL)); // inits random

//  if (argc < 2) {
//    fprintf(stderr, "usage: %s XOR | data_file + parameter_file\n",
//            argv[0]);
//    return EXIT_FAILURE;
//  }
//  else if (strcasecmp("XOR", argv[1]) == 0) populationTest(); // XOR test
//  else if (argc < 3) {
//    fprintf(stderr, "usage: %s data_file parameter_file\n", argv[0]);
//    return EXIT_FAILURE;
//  }
//  else testOneGame(argv[1], argv[2]);

//  testOneGame("/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/data/tictactoe_P1.txt",
//              "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/params/params_1.ini");

  populationTest();

  return 0;
}


/*******************************************************************************
 * computes one game
 ******************************************************************************/

/*
 * computes one game final states file with NEAT
 */
void testOneGame(char * gamePath, char * paramsPath) {
  sParams * params = loadConf(paramsPath);
  testBase * data = loadData(gamePath);

  // cut file extension
  char * gameName = getName(gamePath);

  // create global log file
  if (GlobalLog) createGlobalLogFile(Log_dir);
  if (Verbose) dumpParams(stdout, params);
  if (GlobalLog) dumpParams(GlobalLogFile, params);

  // print game data, list of differents facts and all finals states
  if (UltraVerbose) dumpFacts(data);
  if (UltraVerbose) dumpFinalStates(data);

  // convert data into a simple form to feed NEAT algorithm
  NeatDataBase * simpleData = testToNeatDataBase(data);
  if (UltraVerbose) dumpSimpleData(simpleData);

  // Create a new population
  sPopulation * pop = createPopulation( simpleData->iInputSize,
                                        simpleData->iOutputSize, params );
  
  int iter, indiv, entry, player, spec;
  double error, totalError;  

  // creates first phenotype for each individual
  for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
    int depth = calculateNetDepth(pop->vGenomes[indiv]);
    createPhenotype(pop->vGenomes[indiv], depth);
    if (UltraVerbose)
      dumpPhenotype( pop->vGenomes[indiv]->pPhenotype,
                     pop->vGenomes[indiv]->iId );
  }

  // iterations : compute fitnesses and generates next epoch
  for (iter = 0; iter < params->iNumEpoch; iter++) {
    printf("\n--------------------- epoch %d\n", pop->iGeneration);
    if (Verbose) printf("fitnesses :\n");
    if (GlobalLog)
      fprintf(GlobalLogFile, "\n--------------------- epoch %d\n fitnesses :\n",
              pop->iGeneration);
    for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
      totalError = 0;
      for (entry = 0; entry < simpleData->iNumEntries; entry++) {
        double * outputs = updateANNresponse(pop->vGenomes[indiv]->pPhenotype,
                                             simpleData->mInputs[entry],
                                             simpleData->iInputSize,
                                             simpleData->iOutputSize, SNAPSHOT);
        for (player = 0; player < simpleData->iOutputSize; player++) {
          error = simpleData->mOutputs[entry][player] - outputs[player];
          if (error < 0) error = -error;
          totalError += error;
        }
        free(outputs);
      }
      totalError /= simpleData->iOutputSize * simpleData->iNumEntries;
      pop->vGenomes[indiv]->dFitness = (1 - totalError) * 100;
      if (Verbose) printf("%f\t", pop->vGenomes[indiv]->dFitness);
      if (GlobalLog)
        fprintf(GlobalLogFile, "%f\t", pop->vGenomes[indiv]->dFitness);
    }
    if (Verbose) printf("\n");
    if (GlobalLog) fprintf(GlobalLogFile, "\n");

    // dump every iDumpEvery generations to have regular repports
    if (Dump && pop->iGeneration % params->iDumpEvery == 0)
      createDumpFile(Log_dir, gameName, pop, data, params);

    if (iter < params->iNumEpoch - 1) epoch(pop); // new generation
  }

  // show best genome and fitness
  double max = 0;
  sGenome * maxGen = NULL;
  for (spec = 0; spec < pop->iNumSpecies; spec++) {
    if (pop->vSpecies[spec]->dBestFitness > max) {
      max = pop->vSpecies[spec]->dBestFitness;
      maxGen = pop->vSpecies[spec]->sLeader;
    }
  }
  // final report
  if (Verbose)
    dumpInputsFactsAndGenome(stdout, data, maxGen);
  if (GlobalLog) {
    dumpInputsFactsAndGenome(GlobalLogFile, data, maxGen);
    dumpInnovTable(GlobalLogFile, pop->sInnovTable);
    if (GraphPs) genomeToGraphViz(maxGen, data, pop, Log_dir, gameName);
    else genomeToGraphVizDotFile(maxGen, data, pop, Log_dir, gameName);
  }

  // close files
  fclose(GlobalLogFile);
  // memory free
  free(gameName);
  freePopulation(pop);
  // free game data
  freeSimpleData(simpleData);
  freeData(data);
  free(params);
}

/*******************************************************************************
 * convenient - get filename without extension
 ******************************************************************************/

char * getName(char * path) {
  char name[BUFSIZ];
  int i = 0, j;
  while (*path != '\0') {
    if (*path == '/') i = 0;
    else name[i++] = *path;
    path++;
  }
  for (j = 0; j < i; j++)
  if (name[j] == '.') name[j] = '\0';
  return strdup(name);
}

/*******************************************************************************
 * convert genome graph into graphiz DOT format / or DOT and PS format
 ******************************************************************************/

void genomeToGraphVizDotFile(sGenome * gen, testBase * data, sPopulation * pop,
                             const char * dir, char * gameName) {
  char filePath[512];
  sprintf(filePath, "%s%d_%s_gen%d_fit%f.gv", dir, getpid(), gameName,
          gen->iId, gen->dFitness);
  FILE * GVFile;
  if ((GVFile = fopen(filePath, "w")) == NULL) {
    fprintf(stderr, "error in main() : can't create graphiz DOT file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(GVFile, "/*--------------- %s - fitness %f ---------------*/\n",
          gameName, gen->dFitness);
  int i, j;
  // directed graph
  fprintf(GVFile, "digraph Genome_%d {\n\trankdir=BT;\n\tsize=\"10,8\";\n"
          "\trotate=90;\n\tranksep=\"3 equally\";\n", // /tsplines=polyline;\n
          gen->iId);
  // labels
  fprintf(GVFile, "\t{\n\t\tnode [shape=plaintext, fontsize=16];\n"
          "\t\t\"inputs\" -> \"outputs\";\n\t}\n");
  // rank setting for inputs
  fprintf(GVFile, "\t{\n\t\trank = same; \"inputs\";\n\t\t");
  for (i = 0; i <= data->iNumFacts; i++)
    fprintf(GVFile, "%d; ", i);
  fprintf(GVFile, "\n\t}\n");
  // rank setting for ouputs
  fprintf(GVFile, "\t{\n\t\trank = same; \"outputs\";\n\t\t");
  for (i = data->iNumFacts + 1; i <= data->iNumFacts + data->iNumPlayers; i++)
    fprintf(GVFile, "%d; ", i);
	fprintf(GVFile, "\n\t}\n");
  // rank setting for hidden layers
  for (i = 0; i < pop->iNumDepth; i++) {
    fprintf(GVFile, "\t{\n\t\trank = same;\n\t\t");
    if (i == 0) fprintf(GVFile, "\"inputs\";\n\t\t");
    if (i == 1) fprintf(GVFile, "\"outputs\";\n\t\t");
    for (j = 0; j < gen->iNumNeurons; j++) {
      if (gen->vNeurons[j]->dSplitY == pop->vDepths[i])
        fprintf(GVFile, "%d; ", gen->vNeurons[j]->iId);
    }
    fprintf(GVFile, "\n\t}\n");
  }

  // al neurons
  for (i = 0; i < gen->iNumNeurons; i++) {
    fprintf(GVFile, "\t%d [", gen->vNeurons[i]->iId);
    if (i < data->iNumFacts)
      fprintf(GVFile, "label=\"(%s)\\n", data->vAllFacts[i]);
    else if (i == data->iNumFacts)
      fprintf(GVFile, "label=\"bias\\n");
    else if (i <= data->iNumFacts + data->iNumPlayers)
      fprintf(GVFile, "label=\"score\\nplayer %d\\n", i - data->iNumFacts);
    else
      fprintf(GVFile, "label=\"");
    fprintf(GVFile, "%3.1f\",style=filled,fillcolor=\"0 0 %f\"];\n",
            gen->vNeurons[i]->dSigmoidCurvature,
            gen->vNeurons[i]->dSigmoidCurvature);
  }
  // all links
  for (i = 0; i < gen->iNumLinks; i++) {
    if (gen->vLinks[i]->bEnabled == E_TRUE)
      fprintf(GVFile, "\t%d -> %d [label=\"%+3.1f\",penwidth=%f,color=%s];\n",
              gen->vLinks[i]->iFromNeuron, gen->vLinks[i]->iToNeuron,
              gen->vLinks[i]->dWeight,
              (gen->vLinks[i]->dWeight < 0)? -5 * gen->vLinks[i]->dWeight :
              5 * gen->vLinks[i]->dWeight,
              (gen->vLinks[i]->dWeight < 0)? "blue":"red");
  }
  fprintf(GVFile, "}\n");
  fclose(GVFile);

  //exemple : dot -Tps Desktop/graph3.gv -o Desktop/graph3.ps
  //char psFilePath[512];
  //sprintf(psFilePath, "%s%d_%s_gen%d.ps", dir, getpid(), gameName, gen->iId);
  //execl(Graphiz_dot, "dot", "-Tps", filePath, "-o", psFilePath, NULL);
}

void genomeToGraphViz(sGenome * gen, testBase * data, sPopulation * pop,
                      const char * dir, char * gameName) {
#ifdef GRAPHVIZ
  char string_buffer1[512];
  char string_buffer2[512];
  GVC_t * gvc = gvContext();
  Agnode_t * node1 = NULL;
  Agnode_t * node2 = NULL;
  Agedge_t * edge = NULL;
  int i, j;

  // new graph
  sprintf(string_buffer1, "Genome_%d", gen->iId);
  Agraph_t * g = agopen(string_buffer1, Agstrictdirected, NULL);
  // set attribute of the graph
  agattr(g, AGRAPH, "rankdir", "BT");
  agattr(g, AGRAPH, "size", "10,8");
  agattr(g, AGRAPH, "rotate", "90");
  agattr(g, AGRAPH, "splines", "polyline");
  agattr(g, AGRAPH, "ranksep", "3 equally");
  agattr(g, AGNODE, "style", "filled");
  agattr(g, AGNODE, "shape", "");

  // subgraph legend
  Agraph_t * legend = agsubg(g, NULL, 1);
  agattr(legend, AGNODE, "style", "");
  agattr(legend, AGNODE, "shape", "plaintext");
  agattr(legend, AGNODE, "fontsize", "16");
  Agnode_t * inputs = agnode(legend, "inputs", 1);
  Agnode_t * outputs = agnode(legend, "outputs", 1);
  agedge(g, inputs, outputs, NULL, 1);
  // subgraph layers
  Agraph_t ** subg = calloc(pop->iNumDepth, sizeof(*subg));
  Agraph_t * inputs_subgraph = NULL;
  Agraph_t * outputs_subgraph = NULL;
  for (i = 0; i < pop->iNumDepth; i++) {
    sprintf(string_buffer1, "Subg_%d", i);
    subg[i] = agsubg(g, string_buffer1, 1);
    agsafeset(subg[i], "rank", "same", "");
    if (i == 0) {
      inputs_subgraph = subg[i];
      agnode(inputs_subgraph, "inputs", 1); // inputs
    }
    else if (i == 1) {
      outputs_subgraph = subg[i];
      agnode(outputs_subgraph, "outputs", 1); // outputs
    }
  }
  // all neurons
  for (i = 0; i < gen->iNumNeurons; i++) {
    sprintf(string_buffer1, "%d", gen->vNeurons[i]->iId);
    if (i < data->iNumFacts) {                             // input node
      node1 = agnode(inputs_subgraph, string_buffer1, 1);
      sprintf(string_buffer2, "(%s)\\n%3.1f",
              data->vAllFacts[i], gen->vNeurons[i]->dSigmoidCurvature);
      agset(node1, "label", string_buffer2);
    } else if (i == data->iNumFacts) {                     // bias
      node1 = agnode(inputs_subgraph, string_buffer1, 1);
      sprintf(string_buffer2, "bias\\n%3.1f",
                             gen->vNeurons[i]->dSigmoidCurvature);
      agset(node1, "label", string_buffer2);
    } else if (i <= data->iNumFacts + data->iNumPlayers) { // output node
      node1 = agnode(outputs_subgraph, string_buffer1, 1);
      sprintf(string_buffer2, "score\\nplayer %d\\n%3.1f",
                             i - data->iNumFacts,
                             gen->vNeurons[i]->dSigmoidCurvature);
      agset(node1, "label", string_buffer2);
    } else {                                               // hidden node
      for (j = 0; j < pop->iNumDepth; j++) {
        if (gen->vNeurons[i]->dSplitY == pop->vDepths[j]) {
          node1 = agnode(subg[j], string_buffer1, 1);
          break;
        }
      }
      sprintf(string_buffer2, "%3.1f", gen->vNeurons[i]->dSigmoidCurvature);
      agset(node1, "label", string_buffer2);
    }
    sprintf(string_buffer1, "0 0 %f", gen->vNeurons[i]->dSigmoidCurvature);
    agset(node1, "fillcolor", string_buffer1);
  }
  // all links
  agattr(g, AGEDGE, "label", "");
  agattr(g, AGEDGE, "penwidth", "");
  agattr(g, AGEDGE, "color", "");
  for (i = 0; i < gen->iNumLinks; i++) {
    if (gen->vLinks[i]->bEnabled == E_TRUE) {
      sprintf(string_buffer1, "%d", gen->vLinks[i]->iFromNeuron);
      sprintf(string_buffer2, "%d", gen->vLinks[i]->iToNeuron);
      node1 = agnode(g, string_buffer1, 1);
      node2 = agnode(g, string_buffer2, 1);
      edge = agedge(g, node1, node2, NULL, 1);
      sprintf(string_buffer1, "%+3.1f", gen->vLinks[i]->dWeight);
      agset(edge, "label", string_buffer1);
      sprintf(string_buffer1, "%f",
              (gen->vLinks[i]->dWeight < 0) ? -5 * gen->vLinks[i]->dWeight :
                                               5 * gen->vLinks[i]->dWeight);
      agset(edge, "penwidth", string_buffer1);
      agset(edge, "color", (gen->vLinks[i]->dWeight < 0)? "blue":"red");
    }
  }
  // free subgraph layers vector
  free(subg);

  gvLayout (gvc, g, "dot");
  // create DOT file
  sprintf(string_buffer1, "%s%d_%s_gen%d_fit%f.gv", dir, getpid(), gameName,
                                                    gen->iId, gen->dFitness);
  gvRenderFilename(gvc, g, "gv", string_buffer1);
  // create postscript file
  sprintf(string_buffer2, "%s%d_%s_gen%d_fit%f.ps", dir, getpid(), gameName,
                                                    gen->iId, gen->dFitness);
  gvRenderFilename(gvc, g, "ps", string_buffer2);

  gvFreeLayout(gvc, g);

  agclose(g);
  gvFreeContext(gvc);
#endif
}


/*******************************************************************************
 * log & dump files
 ******************************************************************************/

void createGlobalLogFile(const char * dir) {
  char filePath[512];
  sprintf(filePath, "%s%d__GLOBAL-LOG.txt", dir, getpid());
  if ((GlobalLogFile = fopen(filePath, "w")) == NULL) {
    fprintf(stderr, "error in main() : can't open global log file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(GlobalLogFile, "------------ test bench %d ------------\n", getpid());
}

void createDumpFile(const char * dir, const char * gameName, sPopulation * pop,
                    testBase * data, sParams * params){
  FILE * DumpFile;
  char filePath[512];
  sprintf(filePath, "%s%d_%s_epoch%05d_dump.txt",
          dir, getpid(), gameName, pop->iGeneration);
  if ((DumpFile = fopen(filePath, "w")) == NULL) {
    fprintf(stderr, "error in main() : can't open dump file\n");
    exit(EXIT_FAILURE);
  }
  fprintf(DumpFile, "%s\n\n--------------- epoch %d ---------------\n",
                    gameName, pop->iGeneration);
  dumpSpecies(DumpFile, pop);
  int spec;
  char fileName[512];
  for (spec = 0; spec < pop->iNumSpecies; spec++) {
    dumpGenome(DumpFile, pop->vSpecies[spec]->sLeader);
    sprintf(fileName, "%s_epoch%05d", gameName, pop->iGeneration);
    if (GraphPs)
      genomeToGraphViz(pop->vSpecies[spec]->sLeader, data, pop, Log_dir,
                      fileName);
    else
      genomeToGraphVizDotFile(pop->vSpecies[spec]->sLeader, data, pop, Log_dir,
                             fileName);
    fprintf(DumpFile, "\n");
  }
  fclose(DumpFile);
}

// dump final Genome and facts
void dumpInputsFactsAndGenome(FILE * out, testBase * data, sGenome * maxGen) {
  fprintf(out, "--------------FACTS----------------\n");
  int i;
  for (i = 0; i < data->iNumFacts; i++)
    fprintf(out, "%d - (%s)\n", i, data->vAllFacts[i]);
  fprintf(out, "--------------BEST FITNESS = %f\n", maxGen->dFitness);
  dumpGenome(out, maxGen);
}


/*******************************************************************************
 * test - convenient for debug - XOR
 ******************************************************************************/

/*
 * test with XOR : population, speciation, reproduction
 */
void populationTest() {
  // set globales
  Verbose = 0;
  UltraVerbose = 0;
  GlobalLog = 0;
  Dump = 0;
  GraphPs = 0;
  
  sParams * params = loadConf(NULL);
  dumpParams(stdout, params);

  double data_in[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  double data_out[4] = {0, 1, 1, 0};

  int num_tests = 100;
  int test;
  int winned = 0;
  double average_time = 0;
  for (test = 0; test < num_tests; test++) {
    printf("test %d\n", test);
    // Create a new population
    sPopulation * pop = createPopulation(2, 1, params);

    int indiv;
    // creates first phenotype for each individual
    for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
      //int depth = calculateNetDepth(pop->vGenomes[indiv]);
      createPhenotype(pop->vGenomes[indiv],
                      pop->iNumDepth + pop->vGenomes[indiv]->iNumRecur);
      if (Verbose) dumpGenome(stdout, pop->vGenomes[indiv]);
      if (Verbose) dumpPhenotype(pop->vGenomes[indiv]->pPhenotype, pop->vGenomes[indiv]->iId);
    }
    if (Verbose) dumpInnovTable(stdout, pop->sInnovTable);

    int iter, entry, spec;
    double error, totalError;
    double collected_outputs[4] = {0};

    // iterations : compute fitnesses and generates next epoch
    for (iter = 0; iter < params->iNumEpoch; iter++) {
      if (Verbose) printf("------- epoch %d\n", pop->iGeneration);
      for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
        //printf("indiv %-3d répond = ", pop->vGenomes[indiv]->iId);
        totalError = 0;
        for (entry = 0; entry < 4; entry++) {
          double * outputs = updateANNresponse(pop->vGenomes[indiv]->pPhenotype,
                                               data_in[entry], 3, 1, SNAPSHOT);
          //printf("%6.2f ", outputs[0]);
          error = data_out[entry] - outputs[0];
          if (error < 0) error = -error;
          totalError += error;
          collected_outputs[entry] = outputs[0];
          free(outputs);
        }
        totalError /= 4;
        pop->vGenomes[indiv]->dFitness = (1 - totalError) * 100;
        
        if (collected_outputs[0] < 0.5 && collected_outputs[1] > 0.5 &&
            collected_outputs[2] > 0.5 && collected_outputs[3] < 0.5) {
          printf("\n0 0 => %f - 0 1 => %f - 1 0 => %f - 1 1 => %f - "
                 "erreur=%f - fitness=%f\n", collected_outputs[0],
                 collected_outputs[1], collected_outputs[2], collected_outputs[3],
                 totalError, pop->vGenomes[indiv]->dFitness);
          winned++;
          average_time += iter;
          goto END;
        }
        if (Verbose) printf("%f ", pop->vGenomes[indiv]->dFitness);
      }
      if (Verbose) printf("\n");

      // create next epoch
      if (iter < params->iNumEpoch-1) epoch(pop);

      // find best genome
      double max = 0;
      sGenome * maxGen = NULL;
      for (spec = 0; spec < pop->iNumSpecies; spec++) {
        if (pop->vSpecies[spec]->dBestFitness > max) {
          max = pop->vSpecies[spec]->dBestFitness;
          maxGen = pop->vSpecies[spec]->sLeader;
        }
      }
      // show best genome
      if (Verbose) dumpGenome(stdout, maxGen);
    }
    //dumpInnovTable(stdout, pop->sInnovTable);

  END:

    freePopulation(pop);
  }

  printf("%d tests - solution found %d times - win : %d %% - it tooks an "
         "average of %d generation", num_tests, winned,
         (int) ((double) winned / (double) num_tests * 100),
         (int) (average_time / winned));
  free(params);
}


