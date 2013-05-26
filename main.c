//
//  main.c
//  NeatTest_0.1
//
//  Created by AlineHUF on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//
//  largely based on Mat Buckland 2002

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gvc.h> // graphiz
#include "data.h"
#include "global.h"

/* GLOBALES */
int Verbose = 1;
int UltraVerbose = 0;
int GlobalLog = 1;
int Dump = 1;

FILE * GlobalLogFile = NULL;

//const char * Games_dir = "data/";
//const char * Log_dir = "log/";
const char * Games_dir = "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/data/";
const char * Log_dir = "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/log/";
const char * Graphiz_dot = "/usr/local/bin/dot";


void testOneGame(char * gamePath, char * paramsPath);
char * getName(char * path);
void createGlobalLogFile(const char * dir);
void createDumpFile(const char * dir, const char * gameName, sPopulation * pop,
                    testBase * data, sParams * params);
void genomeToGraphizDotFile(sGenome * gen, testBase * data, sPopulation * pop,
                            const char * dir, char * gameName);
void dumpInputsFactsAndGenome(FILE * out, testBase * data, sGenome * maxGen);
void populationTest(); // XOR test (population, speciation...)

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

  testOneGame("/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/data/tictactoe_P1.txt",
              "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/params/params_1.ini");

  //populationTest();

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
    genomeToGraphizDotFile(maxGen, data, pop, Log_dir, gameName);
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
    genomeToGraphizDotFile(pop->vSpecies[spec]->sLeader, data, pop, Log_dir,
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
 * convert genome graph into graphiz DOT format
 ******************************************************************************/

//void genomeToGraphizDotFile(sGenome * gen, testBase * data, sPopulation * pop,
//                            const char * dir, char * gameName) {
//  char filePath[512];
//  sprintf(filePath, "%s%d_%s_gen%d.gv", dir, getpid(), gameName, gen->iId);
//  FILE * GVFile;
//  if ((GVFile = fopen(filePath, "w")) == NULL) {
//    fprintf(stderr, "error in main() : can't create graphiz DOT file\n");
//    exit(EXIT_FAILURE);
//  }
//  fprintf(GVFile, "/*--------------- %s - fitness %f ---------------*/\n",
//          gameName, gen->dFitness);
//  int i;
//  // directed graph 
//  fprintf(GVFile, "digraph Genome_%d {\n\trankdir=BT;\n\tsize=\"10,8\";\n"
//          "\trotate=90;\nsplines=polyline;\nranksep=\"3 equally\";\n",
//          gen->iId);
//  // labels
//  fprintf(GVFile, "\t{\n\t\tnode [shape=plaintext, fontsize=16];\n"
//                  "\t\t\"inputs\" -> \"outputs\";\n\t}\n");
//  // rank setting for inputs
//  fprintf(GVFile, "\t{\n\t\trank = same; \"inputs\";\n\t\t");
//  for (i = 0; i <= data->iNumFacts; i++)
//    fprintf(GVFile, "%d; ", i);
//  fprintf(GVFile, "\n\t}\n");
//  // rank setting for ouputs
//  fprintf(GVFile, "\t{\n\t\trank = same; \"outputs\";\n\t\t");
//  for (i = data->iNumFacts + 1; i <= data->iNumFacts + data->iNumPlayers; i++)
//    fprintf(GVFile, "%d; ", i);
//	fprintf(GVFile, "\n\t}\n");
//  // al neurons
//  for (i = 0; i < gen->iNumNeurons; i++) {
//    fprintf(GVFile, "\t%d [", gen->vNeurons[i]->iId);
//    if (i < data->iNumFacts)
//      fprintf(GVFile, "label=\"(%s)\\n", data->vAllFacts[i]);
//    else if (i == data->iNumFacts)
//      fprintf(GVFile, "label=\"bias\\n");
//    else if (i <= data->iNumFacts + data->iNumPlayers)
//      fprintf(GVFile, "label=\"score\\nplayer %d\\n", i - data->iNumFacts);
//    else
//      fprintf(GVFile, "label=\"");
//    fprintf(GVFile, "%3.1f\",style=filled,fillcolor=\"0 0 %f\"];\n",
//           gen->vNeurons[i]->dSigmoidCurvature,
//           gen->vNeurons[i]->dSigmoidCurvature);
//  }
//  // all links
//  for (i = 0; i < gen->iNumLinks; i++) {
//    if (gen->vLinks[i]->bEnabled == E_TRUE)
//      fprintf(GVFile, "\t%d -> %d [label=\"%+3.1f\",penwidth=%f,color=%s];\n",
//           gen->vLinks[i]->iFromNeuron, gen->vLinks[i]->iToNeuron,
//           gen->vLinks[i]->dWeight,
//           (gen->vLinks[i]->dWeight < 0)? -5 * gen->vLinks[i]->dWeight :
//                                          5 * gen->vLinks[i]->dWeight,
//           (gen->vLinks[i]->dWeight < 0)? "blue":"red");
//  }
//  fprintf(GVFile, "}\n");
//  fclose(GVFile);
//
//  //exemple : dot -Tps Desktop/graph3.gv -o Desktop/graph3.ps
//  //char psFilePath[512];
//  //sprintf(psFilePath, "%s%d_%s_gen%d.ps", dir, getpid(), gameName, gen->iId);
//  //execl(Graphiz_dot, "dot", "-Tps", filePath, "-o", psFilePath, NULL);
//}

void genomeToGraphizDotFile(sGenome * gen, testBase * data, sPopulation * pop,
                            const char * dir, char * gameName) {
  FILE * fp;
  char string_buffer1[512];
  char string_buffer2[512];
  GVC_t * gvc = gvContext();
  Agnode_t * node1;
  Agnode_t * node2;
  Agedge_t * edge;
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
  Agraph_t * inputs_subgraph;
  Agraph_t * outputs_subgraph;
  for (i = 0; i < pop->iNumDepth; i++) {
    subg[i] = agsubg(g, NULL, 1);
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

  gvLayout(gvc, g, "dot");
  // create DOT file
  sprintf(string_buffer1, "%s%d_%s_gen%d_fit%f.gv", dir, getpid(), gameName,
                                                    gen->iId, gen->dFitness);
  fp = fopen(string_buffer1, "w");
  gvRender(gvc, g, "gv", fp);
  fclose(fp);
  // create postscript file
  sprintf(string_buffer2, "%s%d_%s_gen%d_fit%f.ps", dir, getpid(), gameName,
                                                    gen->iId, gen->dFitness);
  fp = fopen(string_buffer2, "w");
  gvRender(gvc, g, "ps", fp);
  fclose(fp);

  gvFreeLayout(gvc, g);
  agclose(g);

  gvFreeContext(gvc);
}

/*******************************************************************************
 * tests - convenient for debug
 ******************************************************************************/

/*
 * test with XOR : population, speciation, reproduction
 */
void populationTest() {
  sParams * params = loadConf(NULL);
  dumpParams(stdout, params);

  double data_in[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  double data_out[4] = {0, 1, 1, 0};

  // Create a new population
  sPopulation * pop = createPopulation(2, 1, params);

  int indiv;
  // creates first phenotype for each individual
  for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
    //int depth = calculateNetDepth(pop->vGenomes[indiv]);
    createPhenotype(pop->vGenomes[indiv],
                    pop->iNumDepth + pop->vGenomes[indiv]->iNumRecur);
    dumpGenome(stdout, pop->vGenomes[indiv]);
    dumpPhenotype(pop->vGenomes[indiv]->pPhenotype, pop->vGenomes[indiv]->iId);
  }

  int iter, entry, spec;
  double error, totalError;
  
  // iterations : compute fitnesses and generates next epoch
  for (iter = 0; iter < params->iNumEpoch; iter++) {
    printf("------- epoch %d\n", pop->iGeneration);
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
        free(outputs);
      }
      totalError /= 4;
      //printf("- erreur=%6.2f\n", totalError);
      pop->vGenomes[indiv]->dFitness = (1 - totalError) * 100;
      printf("%f ", pop->vGenomes[indiv]->dFitness);
    }
    printf("\n");

    // create next epoch
    epoch(pop);

    // show best genome and outputs
    double max = 0;
    sGenome * maxGen = NULL;
    for (spec = 0; spec < pop->iNumSpecies; spec++) {
      if (pop->vSpecies[spec]->dBestFitness > max) {
        max = pop->vSpecies[spec]->dBestFitness;
        maxGen = pop->vSpecies[spec]->sLeader;
      }
    }
    dumpGenome(stdout, maxGen);
    int depth = calculateNetDepth(maxGen);
    createPhenotype(maxGen, depth);
    printf("outputs = ");
    totalError = 0;
    for (entry = 0; entry < 4; entry++) {
      double * outputs = updateANNresponse(maxGen->pPhenotype,
                                           data_in[entry], 3, 1, SNAPSHOT);
      error = data_out[entry] - outputs[0];
      if (error < 0) error = -error;
      totalError += error;
      printf("%d %d => %8.6f\t",
             (int) data_in[entry][0], (int) data_in[entry][1], outputs[0]);
      free(outputs);
    }
    totalError /= 4;
    printf("- erreur=%6.2f - fitness=%6.2f\n", totalError,
                                               (1 - totalError) * 100);
  }

  freePopulation(pop);
  free(params);
}

