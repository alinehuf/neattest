//
//  population.h
//  NeatTest_0.1
//
//  Created by dexter on 25/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_population_h
#define NeatTest_0_1_population_h

#include "species.h"

// this structure is used in the creation of a network depth lookup table.
typedef struct splitDepth {
  double val;
  int    depth;
  struct splitDepth * cdr;
} sSplitDepth;

typedef struct {
  sGenome ** vGenomes;      // current population (fixed size = iNumIndividuals)
  int iNumGenomes;
  listSpecies * listSpecies; // all the species
  int iNumSpecies;
  sInnovTable * sInnovTable; // to keep track of innovations
  // current generation
  int iGeneration;
	int iNextGenomeId;
  // adjusted fitness scores
  double dTotFitAdj;
  double dAvFitAdj;
  // index into the genomes for the fittest genome
  int iFittestGenome;
  double dBestEverFitness;
  sParams * sParams;
  // this holds the precalculated split depths. They are used
  // to calculate a neurons x/y position for rendering and also
  // for calculating the flush depth of the network when a
  // phenotype is working in 'snapshot' mode.
  sSplitDepth * listSplits;
} sPopulation;

// population
sPopulation * createPopulation(int popSize, int nbIn, int nbOut, sParams * p);
void freePopulation(sPopulation * pop);
void freeSpecies(sPopulation * pop);
// epoch
void epoch(sPopulation * pop, double * fitnessScores, int size);
void resetAndKill(sPopulation * pop, int iNumGensAllowedNoImprovement);
void sortGenomes(sPopulation * pop);
sGenome * tournamentSelection(sPopulation * pop, const int numComparisons);
// speciation
void speciateAndCalculateSpawnLevels(sPopulation * pop);
// crossover
sGenome * crossover(sGenome * mum, sGenome * dad);
// convenient
bool idNotIntoVect(int id, int * vect, int size);
int cmpLinksByInnovIds(const sLinkGene ** a, const sLinkGene ** b);
int cmpNeuronsByIds(const sNeuronGene ** a, const sNeuronGene ** b);
int cmpGenomesByFitness(const sGenome ** a, const sGenome ** b);
// depth lookup table
sSplitDepth * split(sSplitDepth * vecSplits, double low, double high,
                    int depth, int total);
//int calculateNetDepth(sGenome * gen, sSplitDepth * listSplit);
int calculateNetDepth(sGenome * gen);
#endif
