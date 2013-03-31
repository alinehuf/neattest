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
#include "genome.h"
#include "phenotype.h"
#include "population.h"

/*******************************************************************************
 * epoch
 ******************************************************************************/

/*  This function performs one epoch of the genetic algorithm and returns
 *  a vector of pointers to the new phenotypes
 */
sPhenotype * epoch(sPopulation * pop, double * fitnessScores, int size) {
  //first check to make sure we have the correct amount of fitness scores
  if (size != pop->iNumGenomes) {
    fprintf(stderr, "epoch() : error scores/genomes mismatch!");
    exit(EXIT_FAILURE);
  }

  // reset appropriate values and kill off the existing phenotypes and
  // any poorly performing species
  resetAndKill(pop, pop->sParams->iNumGensAllowedNoImprov);
  
  // update the genomes with the fitnesses scored in the last run
  int gen;
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    pop->vGenomes[gen].dFitness = fitnessScores[gen];

  // sort genomes
  sortGenomes(pop);

  // separate the population into species of similar topology, adjust
  // fitnesses and calculate spawn levels
  speciateAndCalculateSpawnLevels(pop);

  // this will hold the new population of genomes
  sGenome * newGenomes;

  // request the offspring from each species. The number of children to
  // spawn is a double which we need to convert to an int.
  int numSpawnedSoFar = 0;
  sGenome baby;
  // now to iterate through each species selecting offspring to be mated and
  // mutated
  listSpecies * curSpc = pop->listSpecies;
  while (curSpc != NULL) {
    // because of the number to spawn from each species is a double
    // rounded up or down to an integer it is possible to get an overflow
    // of genomes spawned. This statement just makes sure that doesn't happen
    if (numSpawnedSoFar < pop->sParams->iNumIndividuals) {
      //this is the amount of offspring this species is required to
      // spawn. Rounded simply rounds the double up or down.
      int numToSpawn = round(curSpc->sSpecies->dSpawnsRqd);

      bool bChosenBestYet = FALSE;

      while (numToSpawn--) {
        // first grab the best performing genome from this species and transfer
        // to the new population without mutation. This provides per species
        // elitism
        if (!bChosenBestYet) {
          baby = curSpc->sSpecies->sLeader;
          bChosenBestYet = TRUE;
        } else {
          // if the number of individuals in this species is only one
          // then we can only perform mutation
          if (curSpc->sSpecies->iNumMembers == 1) {
            // spawn a child
            baby = randomSpawn(curSpc->sSpecies, pop->sParams->dSurvivalRate);
          }
          //if greater than one we can use the crossover operator
          else {
            // spawn1
            sGenome g1 = 
                     randomSpawn(curSpc->sSpecies, pop->sParams->dSurvivalRate);
            if (randFloat() < pop->sParams->dCrossoverRate) {
              // spawn2, make sure it's not the same as g1
              sGenome g2 =
                     randomSpawn(curSpc->sSpecies, pop->sParams->dSurvivalRate);

              // number of attempts at finding a different genome
              int numAttempts = pop->sParams->iNumTrysToFindMate;
              while (g1.iId == g2.iId && numAttempts--)
                g2 = randomSpawn(curSpc->sSpecies, pop->sParams->dSurvivalRate);
              if (g1.iId != g2.iId)
                baby = crossover(g1, g2);
            } else {
              baby = g1;
            }
          }

          // adjust new genome id
          baby.iId = pop->iNextGenomeId++;
//
//          //now we have a spawned child lets mutate it! First there is the
//          //chance a neuron may be added
//          if (baby.NumNeurons() < CParams::iMaxPermittedNeurons)
//          {
//            baby.AddNeuron(CParams::dChanceAddNode,
//                           *m_pInnovation,
//                           CParams::iNumTrysToFindOldLink);
//          }
//
//          //now there's the chance a link may be added
//          baby.AddLink(CParams::dChanceAddLink,
//                       CParams::dChanceAddRecurrentLink,
//                       *m_pInnovation,
//                       CParams::iNumTrysToFindLoopedLink,
//                       CParams::iNumAddLinkAttempts);
//
//          //mutate the weights
//          baby.MutateWeights(CParams::dMutationRate,
//                             CParams::dProbabilityWeightReplaced,
//                             CParams::dMaxWeightPerturbation);
//
//          baby.MutateActivationResponse(CParams::dActivationMutationRate,
//                                        CParams::dMaxActivationPerturbation);
        }
//
//        //sort the babies genes by their innovation numbers
//        baby.SortGenes();
//
//        //add to new pop
//        NewPop.push_back(baby);
//
//        ++NumSpawnedSoFar;
//
//        if (NumSpawnedSoFar == CParams::iNumSweepers)
//        {
//          NumToSpawn = 0;
//        }
//
      } //end while
    } //end if
    curSpc = curSpc->next;
  }


  // if there is an underflow due to the rounding error and the amount
  // of offspring falls short of the population size additional children
  // need to be created and added to the new population. This is achieved
  // simply, by using tournament selection over the entire population.
//  if (NumSpawnedSoFar < CParams::iNumSweepers)
//  {
//
//    //calculate amount of additional children required
//    int Rqd = CParams::iNumSweepers - NumSpawnedSoFar;
//
//    //grab them
//    while (Rqd--)
//    {
//      NewPop.push_back(TournamentSelection(m_iPopSize/5));
//    }
//  }
//
//  //replace the current population with the new one
//  m_vecGenomes = NewPop;
//
//  //create the new phenotypes
//  vector<CNeuralNet*> new_phenotypes;
//
//  for (gen=0; gen<m_vecGenomes.size(); ++gen)
//  {
//    //calculate max network depth
//    int depth = CalculateNetDepth(m_vecGenomes[gen]);
//    
//    CNeuralNet* phenotype = m_vecGenomes[gen].CreatePhenotype(depth);
//    
//    new_phenotypes.push_back(phenotype);
//  }
//  
//  //increase generation counter
//  ++m_iGeneration;
//  
//  return new_phenotypes;

}


/* resets some values ready for the next epoch, kills off all the phenotypes
 * and any poorly performing species.
 */
void resetAndKill(sPopulation * pop, int iNumGensAllowedNoImprovement) {
  pop->dTotFitAdj = 0;
  pop->dAvFitAdj  = 0;

  //purge the species
  listSpecies * curSpec = pop->listSpecies;
  while (curSpec != NULL) {
    purgeSpecies(curSpec->sSpecies);

    // kill off species if not improving and if not the species which contains
    // the best genome found so far
    if ( curSpec->sSpecies->iGensNoImprovement > iNumGensAllowedNoImprovement
         && curSpec->sSpecies->dBestFitness < pop->dBestEverFitness ) {
      curSpec = removeOneSpecies(curSpec);
    } else {
       curSpec = curSpec->next;
    }
  }
  //we can also delete the phenotypes
  int gen;
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    deletePhenotype(&pop->vGenomes[gen]);
}

/* sorts the population into descending fitness, 
 * and updates any fitness statistics accordingly
 */
void sortGenomes(sPopulation * pop) {
  // sort the genomes according to their unadjusted (no fitness sharing)
  // fitnesses
  qsort(pop->vGenomes, pop->iNumGenomes, sizeof(sGenome),
        (int (*) (const void *, const void *)) cmpGenomesByFitness);
  // is the best genome this generation the best ever?
  if (pop->vGenomes[0].dFitness > pop->dBestEverFitness)
    pop->dBestEverFitness = pop->vGenomes[0].dFitness;
}

/*******************************************************************************
 * speciation and spawn level
 ******************************************************************************/

/* separates each individual into its respective species by calculating
 * a compatibility score with every other member of the population and
 * niching accordingly. The function then adjusts the fitness scores of
 * each individual by species age and by sharing and also determines
 * how many offspring each individual should spawn.
 */
void speciateAndCalculateSpawnLevels(sPopulation * pop) {
  int gen;
  listSpecies * curSpc;
  bool bAdded;
  
  // iterate through each genome and speciate
  for (gen = 0; gen < pop->iNumGenomes; gen++) {
    bAdded = FALSE;
    // calculate its compatibility score with each species leader
    // if compatible add to species. If not, create a new species
    curSpc = pop->listSpecies;
    while (curSpc != NULL) {
      double compatibility = getCompatibilityScore(pop->vGenomes[gen],
                                       curSpc->sSpecies->sLeader, pop->sParams);
      // if this individual is similar to this species add to species
      if (compatibility <= pop->sParams->dCompatibilityThreshold) {
        addMember(curSpc->sSpecies, &pop->vGenomes[gen]);
        bAdded = TRUE;
        break;
      }
      if (curSpc->next == NULL) break;
      
    }
    // if we have not found a compatible species, let's create a new one
    if (!bAdded)
      addOneSpecies(curSpc, &pop->vGenomes[gen], pop->iNumSpecies++,
                                                 pop->sParams->iNumIndividuals);
  }

  // now all the genomes have been assigned a species the fitness scores
  // need to be adjusted to take into account sharing and species age
  curSpc = pop->listSpecies;
  while (curSpc != NULL) {
    adjustFitnesses(curSpc->sSpecies, pop->sParams);
    curSpc = curSpc->next;
  }

  // calculate new adjusted total & average fitness for the population
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    pop->dTotFitAdj += pop->vGenomes[gen].dAjustedFitness;

  pop->dAvFitAdj = pop->dTotFitAdj / pop->iNumGenomes;

  // calculate how many offspring each member of the population should spawn
  for (gen = 0; gen < pop->iNumGenomes; gen++) {
    double toSpawn = pop->vGenomes[gen].dAjustedFitness / pop->dAvFitAdj;
    pop->vGenomes[gen].dAmountToSpawn = toSpawn;
  }

  // iterate through all the species and calculate how many offspring
  // each species should spawn
  curSpc = pop->listSpecies;
  while (curSpc != NULL) {
    speciesSpawnAmount(curSpc->sSpecies);
    curSpc = curSpc->next;
  }
}


/*******************************************************************************
 * crossover
 ******************************************************************************/

sGenome crossover(sGenome mum, sGenome dad) {
  typedef enum {MUM, DAD, NONE} parent_type;
  // first, calculate the genome we will using the disjoint/excess genes from.
  // This is the fittest genome.
  parent_type best;

  // if they are of equal fitness use the shorter (because we want to keep
  // the networks as small as possible)
  if (mum.dFitness == dad.dFitness) {
    // if they are of equal fitness and length just choose one at random
    if (mum.iNumLinks == dad.iNumLinks)
      best = (parent_type) randInt(0, 1);
    else if (mum.iNumLinks < dad.iNumLinks)
      best = MUM;
    else
      best = DAD;
  } else if (mum.dFitness > dad.dFitness) {
    best = MUM;
  } else {
    best = DAD;
  }

  // create en empty genome for the baby (id will be adjusted later)
  sGenome baby = createEmptyGenome(-1, mum.iNumInputs, mum.iNumOuputs);
  // this will hold a copy of the gene we wish to add at each step
  sLinkGene selectedGene;
  selectedGene.iInnovId = -1;     // to control that the choice of the gene
  parent_type selectedFrom = NONE;  // has been done

  int curMum = 0;
  int curDad = 0;
  // step through each parents genes until we reach the end of both
  while (curMum < mum.iNumLinks || curDad < dad.iNumLinks) {
    
    // the end of mum's genes have been reached
    if (curMum >= mum.iNumLinks && curDad < dad.iNumLinks) {
      // if dad is fittest
      if (best == DAD) { //add dads genes
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto dad's next gene
      curDad++;
    }
    // the end of dads's genes have been reached
    else if ( curDad >= dad.iNumLinks && curMum < mum.iNumLinks) {
      // if mum is fittest
      if (best == MUM) { //add mums genes
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
       }
      // move onto mum's next gene
      curMum++;
    }
    // if mums innovation number is less than dads
    else if (mum.vLinks[curMum].iInnovId < dad.vLinks[curDad].iInnovId) {
      // if mum is fittest add gene
      if (best == MUM) {
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
      }
      // move onto mum's next gene
      curMum++;
    }
    // if dads innovation number is less than mums
    else if (dad.vLinks[curDad].iInnovId < mum.vLinks[curMum].iInnovId) {
      // if dad is fittest add gene
      if (best == DAD) {
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto dad's next gene
      curDad++;
    }
    // if innovation numbers are the same
    else if (mum.vLinks[curMum].iInnovId == dad.vLinks[curDad].iInnovId)
    {
      // grab a gene from either parent
      if (randFloat() < 0.5f) {
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
      } else {
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto next gene of each parent
      curMum++;
      curDad++;
    }

    // this should never happen because the gene *should* have already been
    // selected
    if (selectedFrom == NONE || selectedGene.iInnovId == -1) {
      fprintf(stderr, "crossover() : error selectedGene and selectedFrom " \
                      "should have a value\n");
      exit(EXIT_FAILURE);
    }

    // add the selected gene if not already added
    if ( baby.iNumLinks == 0 ||
         baby.vLinks[baby.iNumLinks-1].iInnovId != selectedGene.iInnovId )
      genomeAddLink(&baby, &selectedGene);

    // Check if we already have the nodes referred to in SelectedGene.
    // If not, they need to be added.
    if (!alreadyHaveThisNeuronId(baby, selectedGene.iFromNeuron)) {
      if (selectedFrom == MUM)
        genomeAddNeuron(&baby,
                   &mum.vNeurons[getNeuronPos(mum, selectedGene.iFromNeuron)] );
      else
        genomeAddNeuron(&baby,
                   &dad.vNeurons[getNeuronPos(dad, selectedGene.iFromNeuron)] );
    }
    if (!alreadyHaveThisNeuronId(baby, selectedGene.iToNeuron)) {
      if (selectedFrom == MUM)
        genomeAddNeuron(&baby,
                     &mum.vNeurons[getNeuronPos(mum, selectedGene.iToNeuron)] );
      else
        genomeAddNeuron(&baby,
                     &dad.vNeurons[getNeuronPos(dad, selectedGene.iToNeuron)] );
    }
  } // end while

  //now sort the neurons into order
  qsort(baby.vNeurons, baby.iNumNeurons, sizeof(sNeuronGene),
        (int (*) (const void *, const void *)) cmpNeuronsByIds);

  return baby;
}

/*******************************************************************************
 * convenient
 ******************************************************************************/

bool idNotIntoVect(int id, int * vect, int size) {
  int i;
  for (i = 0; i < size; i++) {
    if (id == vect[i]) return FALSE;
  }
  return TRUE;
}

int cmpNeuronsByIds(const sNeuronGene * a, const sNeuronGene * b) {
  if (a->iId == b->iId) return 0;
  else if (a->iId < b->iId) return -1;
  else return 1;
}

int cmpGenomesByFitness(const sGenome * a, const sGenome * b) {
  if (a->dFitness == b->dFitness) return 0;
  else if (a->dFitness > b->dFitness) return -1;
  else return 1;
}
