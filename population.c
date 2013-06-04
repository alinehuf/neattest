//
//  population.c
//  NeatTest_0.1
//
//  Created by dexter on 25/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

extern int UltraVerbose;
extern int Verbose;
extern int GlobalLog;
extern FILE * GlobalLogFile;

sPopulation * createPopulation(int nbIn, int nbOut, sParams * params) {
  int i;
  sPopulation * pop = malloc(sizeof(*pop));
  // genomes
  pop->iNumGenomes = params->iNumIndividuals;
  pop->iNextGenomeId = 0;
  pop->vGenomes = calloc(pop->iNumGenomes, sizeof(*pop->vGenomes));
  //create the population of genomes
  for (i = 0; i < params->iNumIndividuals; i++) {
    pop->vGenomes[i] = createInitialGenome(pop->iNextGenomeId++, nbIn, nbOut);
    // DUMP GENOMES
    if (UltraVerbose) dumpGenome(stdout, pop->vGenomes[i]);
  }
  // the species
  pop->iTotalSpecies = INIT_VECT_SIZE;
  pop->vSpecies = calloc(pop->iTotalSpecies, sizeof(*pop->vSpecies));
  pop->iNumSpecies = 0;
  pop->iNextSpeciesId = 0;
  // to keep track of innovations (take the first genome as model)
  pop->sInnovTable = createNewInnovTable(pop->vGenomes[0]->vNeurons,
                                        pop->vGenomes[0]->iNumNeurons,
                                        pop->vGenomes[0]->vLinks,
                                        pop->vGenomes[0]->iNumLinks);
  // DUMP INNOVATION TABLE
  if (UltraVerbose) dumpInnovTable(stdout, pop->sInnovTable);
  //current generation
  pop->iGeneration = 1;
  //adjusted fitness scores
  pop->dTotFitAdj = 0;
  pop->dAvFitAdj = 0;
  //index into the genomes for the fittest genome
  pop->iFittestGenome = 0;
  pop->dBestEverFitness = 0;
  pop->sParams = params;
  // depth lookup table
  pop->iTotalDepths = INIT_VECT_SIZE;
  pop->vDepths = calloc(pop->iTotalDepths, sizeof(*pop->vDepths));
  pop->iNumDepth = 0;
  addDepth(pop, 0);
  addDepth(pop, 1);
  return pop;
}

void freePopulation(sPopulation * pop) {
  int i;
  for (i = 0; i < pop->iNumGenomes; i++) {
    freeGenome(pop->vGenomes[i]);
  }
  free(pop->vGenomes);
  free(pop->vDepths);
  freeInnovTable(pop->sInnovTable);
  freeSpecies(pop);
  free(pop);
}

/*******************************************************************************
 * epoch
 ******************************************************************************/

/*  This function performs one epoch of the genetic algorithm and returns
 *  a vector of pointers to the new phenotypes
 */
void epoch(sPopulation * pop) {
  int spec, gen;
  // reset appropriate values and kill off the existing phenotypes and
  // any poorly performing species
  resetAndKill(pop);

  // sort genomes
  sortGenomes(pop);

  // separate the population into species of similar topology, adjust
  // fitnesses and calculate spawn levels
  speciateAndCalculateSpawnLevels(pop);
  if (Verbose) dumpSpecies(stdout, pop);
  if (GlobalLog) dumpSpecies(GlobalLogFile, pop);
  
  // this will hold the new population of genomes
  sGenome ** newGenomes = calloc( pop->sParams->iNumIndividuals,
                                  sizeof(*newGenomes) );
  // request the offspring from each species. The number of children to
  // spawn is a double which we need to convert to an int.
  int numSpawnedSoFar = 0;
  sGenome * baby;
  // now to iterate through each species selecting offspring to be mated and
  // mutated
  for (spec = 0; spec < pop->iNumSpecies; spec++) {
    // because of the number to spawn from each species is a double
    // rounded up or down to an integer it is possible to get an overflow
    // of genomes spawned. This statement just makes sure that doesn't happen
    if (numSpawnedSoFar < pop->sParams->iNumIndividuals) {
      //this is the amount of offspring this species is required to
      // spawn. Round simply rounds the double up or down.
      int numToSpawn = round(pop->vSpecies[spec]->dSpawnsRqd);

      bool bChosenBestYet = E_FALSE;

      while (numToSpawn--) {
        // first grab the best performing genome from this species and transfer
        // to the new population without mutation. This provides per species
        // elitism
        // => Stanley recommends to do this only if the species has more than
        //    five networks, in Neat1.1 he used the futur number of individual
        //    instead of the current one... 
        if (!bChosenBestYet) { //  && numToSpawn > 5
          baby = copyGenome(pop->vSpecies[spec]->sLeader);
          bChosenBestYet = E_TRUE;
        } else {
          // if the number of individuals in this species is only one
          // then we can only perform mutation
          if (pop->vSpecies[spec]->iNumMembers == 1) {
            // spawn a child
            baby = copyGenome(randomAmongBest( pop->vSpecies[spec],
                                               pop->sParams->dSurvivalRate ));
          }
          //if greater than one we can use the crossover operator
          else {
            // spawn1
            sGenome * g1 = randomAmongBest( pop->vSpecies[spec],
                                            pop->sParams->dSurvivalRate );
            if (randFloat() < pop->sParams->dCrossoverRate) {
              // spawn2, make sure it's not the same as g1
              sGenome * g2 = randomAmongBest( pop->vSpecies[spec],
                                              pop->sParams->dSurvivalRate );
              // number of attempts at finding a different genome
              int numAttempts = pop->sParams->iNumTrysToFindMate;
              while (g1->iId == g2->iId && numAttempts--)
                g2 = randomAmongBest( pop->vSpecies[spec],
                                      pop->sParams->dSurvivalRate);
              if (g1->iId != g2->iId)
                baby = crossover(g1, g2);
              else
                baby = copyGenome(g1); // crossover fail
            } else {
              baby = copyGenome(g1);   // no crossover
            }
          }

          // adjust new genome id
          baby->iId = pop->iNextGenomeId++;

          // now we have a spawned child lets mutate it! First there is the
          // chance a neuron may be added
          if (randFloat() < pop->sParams->dChanceAddNode)
            addNeuron(baby, pop);
          // now there's the chance a link may be added
          else if (randFloat() < pop->sParams->dChanceAddLink)
            addLink(baby, pop);
          else {
            // mutate the weights
            mutateWeigth(baby, pop->sParams);
            // mutate the activation response
            mutateActivationResponse(baby, pop->sParams);
          }
       } // end choice of a baby

        //sort the babies genes by their innovation numbers
        qsort(baby->vLinks, baby->iNumLinks, sizeof(sLinkGene *),
              (int (*) (const void *, const void *)) cmpLinksByInnovIds);

        //add to new pop
        newGenomes[numSpawnedSoFar++] = baby;
        
        if (numSpawnedSoFar == pop->sParams->iNumIndividuals)
          goto newGenerationReady;

      } // end while not enougth babies
    } // end if to much babies
  } // next species

  // if there is an underflow due to the rounding error and the amount
  // of offspring falls short of the population size additional children
  // need to be created and added to the new population. This is achieved
  // simply, by using tournament selection over the entire population.
  if (numSpawnedSoFar < pop->sParams->iNumIndividuals) {
    //calculate amount of additional children required
    int rqd = pop->sParams->iNumIndividuals - numSpawnedSoFar;
    //grab them
    while (rqd--) {
      sGenome * chosenOne = tournamentSelection(pop, pop->iNumGenomes / 5);
      newGenomes[numSpawnedSoFar++] = copyGenome(chosenOne);
    }
  }

newGenerationReady:

  // free the old vector of genomes
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    freeGenome(pop->vGenomes[gen]); // free the genome itself
  free(pop->vGenomes);              // free the vector containing the genomes
  
  //replace the current population with the new one
  pop->vGenomes = newGenomes;

  //create the new phenotypes
  for (gen = 0; gen < pop->iNumGenomes; gen++) {
    if (UltraVerbose) dumpGenome(stdout, pop->vGenomes[gen]);
    //calculate max network depth
    //int depth = calculateNetDepth(pop->vGenomes[gen]);
    createPhenotype(pop->vGenomes[gen],
                    pop->iNumDepth + pop->vGenomes[gen]->iNumRecur);
    if (UltraVerbose)
      dumpPhenotype(pop->vGenomes[gen]->pPhenotype, pop->vGenomes[gen]->iId);
  }

  //increase generation counter
  pop->iGeneration++;
}

/* resets some values ready for the next epoch, kills off all the phenotypes
 * and any poorly performing species.
 */
void resetAndKill(sPopulation * pop) {
  pop->dTotFitAdj = 0;
  pop->dAvFitAdj  = 0;

  // control the number of species
  if (pop->iNumSpecies > pop->sParams->iMaxSpecies)
    pop->sParams->dCompatibilityThreshold += 0.005;

  // purge the species
  int i = 0;
  while (i < pop->iNumSpecies) {
    purgeSpecies(pop->vSpecies[i]);
    // kill off species if not improving and if not the species which contains
    // the best genome found so far
    if ( pop->vSpecies[i]->iGensNoImprovement >
         pop->sParams->iNumGensAllowedNoImprov
         && pop->vSpecies[i]->dBestFitness < pop->dBestEverFitness )
      removeOneSpecies(pop, i);
    else i++;
  }
  // we can also delete the phenotypes
  int gen;
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    freePhenotype(pop->vGenomes[gen]);
}

/* sorts the population into descending fitness, 
 * and updates any fitness statistics accordingly
 */
void sortGenomes(sPopulation * pop) {
  // sort the genomes according to their unadjusted (no fitness sharing)
  // fitnesses
  qsort(pop->vGenomes, pop->iNumGenomes, sizeof(sGenome *),
        (int (*) (const void *, const void *)) cmpGenomesByFitness);
  // is the best genome this generation the best ever?
  if (pop->vGenomes[0]->dFitness > pop->dBestEverFitness)
    pop->dBestEverFitness = pop->vGenomes[0]->dFitness;
}

// realises a tournament selection by comparing numComparisons times the firness
// of random genomes in the population and selecting the best
sGenome * tournamentSelection(sPopulation * pop, const int numComparisons) {
  double bestFitnessSoFar = 0;
  int chosenOne = 0;
  // select numComparisons members from the population at random 
  // testing against the best found so far
  int i;
  for (i = 0; i < numComparisons; i++) {
    int thisTry = randInt(0, pop->iNumGenomes - 1);

    if (pop->vGenomes[thisTry]->dFitness > bestFitnessSoFar) {
      chosenOne = thisTry;
      bestFitnessSoFar = pop->vGenomes[thisTry]->dFitness;
    }
  }
  //return the champion
  return pop->vGenomes[chosenOne];
}

/*******************************************************************************
 * convenient
 ******************************************************************************/

//bool idNotIntoVect(int id, int * vect, int size) {
//  int i;
//  for (i = 0; i < size; i++) {
//    if (id == vect[i]) return E_FALSE;
//  }
//  return E_TRUE;
//}

/* to sort link genes by their innovation id
 */
int cmpLinksByInnovIds(const sLinkGene ** a, const sLinkGene ** b) {
  if ((*a)->iInnovId == (*b)->iInnovId) return 0;
  else if ((*a)->iInnovId < (*b)->iInnovId) return -1;
  else return 1;
}

/* to sort neuron genes by their id
 */
int cmpNeuronsByIds(const sNeuronGene ** a, const sNeuronGene ** b) {
  if ((*a)->iId == (*b)->iId) return 0;
  else if ((*a)->iId < (*b)->iId) return -1;
  else return 1;
}

/* to sort genome by their fitness score 
 * in descending order */
int cmpGenomesByFitness(const sGenome ** a, const sGenome ** b) {
  if ((*a)->dFitness == (*b)->dFitness) return 0;
  else if ((*a)->dFitness > (*b)->dFitness) return -1;
  else return 1;
}

/*******************************************************************************
 * depth lookup table
 ******************************************************************************/

/* add a Y level in the global list for depths of the population
 * usefull to update ANN response
 */
void addDepth(sPopulation * pop, double splitY) {
  int i;
  for (i = 0; i < pop->iNumDepth; i++)
    if (pop->vDepths[i] == splitY) return;
  if (pop->iNumDepth == pop->iTotalDepths) {
    pop->iTotalDepths += INIT_VECT_SIZE;
    pop->vDepths = realloc(pop->vDepths,
                           pop->iTotalDepths * sizeof(*pop->vDepths));
  }  
  pop->vDepths[pop->iNumDepth++] = splitY;
}

/* get all the different level on Y axis of the nodes and
 * returns the depth of the network based on this figure
 */
int calculateNetDepth(sGenome * gen) {
  int i, j;
  double depths[BUFSIZ];
  int numDepths = 0;
  int numRecurrent = 0;
  for (i = 0; i < gen->iNumNeurons; i++) {
    if (gen->vNeurons[i]->bRecurrent == E_TRUE) numRecurrent++;
    for (j = 0; j < numDepths; j++)
      if (gen->vNeurons[i]->dSplitY == depths[j]) goto end;
    depths[numDepths++] = gen->vNeurons[i]->dSplitY;
    end:
    //printf("depth = %d, level = %6.2f\n", numDepths-1, depths[numDepths-1])
    ;
  }
  return numDepths + numRecurrent;
}
