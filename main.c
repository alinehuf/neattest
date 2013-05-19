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
#include "data.h"
#include "genome.h"
#include "population.h"

/* GLOBALES */
int UltraVerbose = 1;
//const char * games_dir = "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/data/";
//const char * log_dir = "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/log/";
const char * games_dir = "data/";
const char * log_dir = "log/";


void testOneGame(char * gamePath, char * paramsPath);
void genomeTest(); // first test of genomes
void populationTest(); // second test (population, speciation...)

/*******************************************************************************
 * MAIN
 ******************************************************************************/

int main(int argc, char * argv[]) {
  srand((unsigned int) time(NULL)); // inits random
  
  //genomeTest(); // first test
  //populationTest(); // second test
  //testOneGame("/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/data/tictactoe_P1.txt", "params_1game.ini");
  testOneGame("data/tictactoe_P1.txt", "params_1game.ini");

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
  // print game data, list of differents facts and all finals states
  if (UltraVerbose) dumpFacts(data);
  if (UltraVerbose) dumpFinalStates(data);

  // convert data into a simple form to feed NEAT algorithm
  NeatDataBase * simpleData = testToNeatDataBase(data);
  if (UltraVerbose) dumpSimpleData(simpleData);

  // Create a new population
  sPopulation * pop = createPopulation( params->iNumIndividuals,
                                        simpleData->iInputSize,
                                        simpleData->iOutputSize, params );
  
  // allocate memory for the fitnesses of each individual
  double * fitnesses = calloc(params->iNumIndividuals, sizeof(*fitnesses));
  int iter, indiv, entry, player;
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
  for (iter = 0; iter < 1; iter++) { // params->iNumEpoch
    printf("------- epoch %d\n", pop->iGeneration);
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
      fitnesses[indiv] = (1 - totalError) * 100;
    }

    for (indiv = 0; indiv < params->iNumIndividuals; indiv++) {
      printf("%f ", fitnesses[indiv]);
    }
    printf("\n");

    epoch(pop, fitnesses, params->iNumIndividuals);
  }

  // memory free
  free(fitnesses);
  freePopulation(pop);
  // free game data
  freeSimpleData(simpleData);
  freeData(data);
  //free(params);
}


/*******************************************************************************
 * tests - convenient for debug
 ******************************************************************************/

/*
 * first test of files : params, utils, genes, genome, innovation, phenotype
 */
void genomeTest() {
  sParams * params = loadConf(NULL);
  
  sGenome * gen1 = createInitialGenome(0, 3, 2);
  sGenome * gen2 = createInitialGenome(1, 3, 2);
  sInnovTable * innovTable = createNewInnovTable(gen1->vNeurons,
                              gen1->iNumNeurons, gen1->vLinks, gen1->iNumLinks);
  dumpGenome(gen1);
  dumpGenome(gen2);
  int i;
  for (i = 0; i < 100; i ++) {
    addNeuron(gen1, params->dChanceAddNode, innovTable,
              params->iNumTrysToFindOldLink);
    addNeuron(gen2, params->dChanceAddNode, innovTable,
              params->iNumTrysToFindOldLink);
  }
  for (i = 0; i < 100; i ++) {
    addLink(gen1, params->dChanceAddLink, params->dChanceAddRecurrentLink,
            innovTable, params->iNumTrysToFindLoop, params->iNumTrysToAddLink);
    addLink(gen2, params->dChanceAddLink, params->dChanceAddRecurrentLink,
            innovTable, params->iNumTrysToFindLoop, params->iNumTrysToAddLink);
  }
  puts("------------MUTATION--------");
  dumpGenome(gen1);
  dumpGenome(gen2);
  printf("compatibility : %f \n",
         getCompatibilityScore(gen1, gen2, params));
  puts("------------CROSSOVER--------");
  sGenome * gen3 = crossover(gen1, gen2);
  gen3->iId = 2;
  dumpGenome(gen3);

  puts("------------INNOVATIONS--------");
  dumpInnovTable(innovTable);

  // free memory
  freeGenome(gen1);
  freeGenome(gen2);
  
  puts("------------PHENOTYPE--------"); // CREER LE PHENOTYPE - CALCULER LA SORTIE
  createPhenotype(gen3, 3);
  dumpPhenotype(gen3->pPhenotype, gen3->iId);
  double inputs[] = {0.3, 0.5, 0.7};
  double * outputs = updateANNresponse(gen3->pPhenotype, inputs,3, 2, SNAPSHOT);
  printf("outputs : ");
  for(i = 0; i < 2; i++) printf("%6.2f ", outputs[i]);
  printf("\n");

  // free memory
  free(outputs);
  freeGenome(gen3);
  freeInnovTable(innovTable);
  free(params);
}

/*
 * second test : population, speciation, reproduction
 */
void populationTest() {
  sParams * params = loadConf(NULL);

  double data_in[4][3] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
  double data_out[4] = {0, 1, 1, 0};

  // Create a new population
  sPopulation * pop = createPopulation(params->iNumIndividuals, 2, 1, params);

  // allocate memory for the fitnesses of each individual
  double * fitnesses = calloc(pop->iNumGenomes, sizeof(*fitnesses));

  int indiv;
  // creates first phenotype for each individual
  for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
    int depth = calculateNetDepth(pop->vGenomes[indiv]);
    createPhenotype(pop->vGenomes[indiv], depth);
    dumpGenome(pop->vGenomes[indiv]);
    dumpPhenotype(pop->vGenomes[indiv]->pPhenotype, pop->vGenomes[indiv]->iId);
  }

  int iter, entry;
  double error, totalError;
  
  // iterations : compute fitnesses and generates next epoch
  for (iter = 0; iter < 2; iter++) {
    printf("------- epoch %d\n", pop->iGeneration);
    for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
      totalError = 0;
      for (entry = 0; entry < 4; entry++) {
        double * outputs = updateANNresponse(pop->vGenomes[indiv]->pPhenotype,
                                             data_in[entry], 3, 1, SNAPSHOT);
        error = data_out[entry] - outputs[0];
        if (error < 0) error = -error;
        totalError += error;
        free(outputs);
      }
      totalError /= 4;
      fitnesses[indiv] = (1 - totalError) * 100;
    }

    for (indiv = 0; indiv < pop->iNumGenomes; indiv++) {
      printf("%f ", fitnesses[indiv]);
    }
    printf("\n");

//    puts("------------TEST SPECIES--------");
//    pop->listSpecies = addOneSpecies( pop->listSpecies, pop->vGenomes[0], 0,
//                                      pop->iNumGenomes );
//    freeSpecies(pop);

    epoch(pop, fitnesses, pop->iNumGenomes);
  }

  free(fitnesses);
  freePopulation(pop);
  free(params);
}

