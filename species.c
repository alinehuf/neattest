//
//  species.c
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/*******************************************************************************
 * gestion of species
 ******************************************************************************/

/* creates an instance of a new species. A local copy of the initializing genome
 * is kept in sLeader and the first element of vMembers is a pointer to that 
 * genome.
 */
sSpecies * createSpecies(sGenome * firstOrg, int specId, int iNumIndividuals) {
  sSpecies * spec = (sSpecies *) malloc(sizeof(*spec));
  spec->sLeader = copyGenome(firstOrg);
  spec->vMembers = (sGenome* *) calloc(iNumIndividuals,sizeof(*spec->vMembers));
  spec->vMembers[0] = firstOrg;
  spec->iNumMembers = 1;
  spec->iTotalMembers = iNumIndividuals;
  spec->iSpeciesId = specId;
  spec->dBestFitness = firstOrg->dFitness;
  spec->iGensNoImprovement = 0;
  spec->iAge = 0;
  spec->dSpawnsRqd = 0;
  return spec;
}

/* free all species */
void freeSpecies(sPopulation * pop) {
  int i;
  for (i = 0; i < pop->iNumSpecies; i++) {
    free(pop->vSpecies[i]->vMembers);
    freeGenome(pop->vSpecies[i]->sLeader);
    free(pop->vSpecies[i]);
  }
  free(pop->vSpecies);
}

/* adds a new member to this species and updates the member variables
 * accordingly
 */
void addMember(sSpecies * spec, sGenome * newMember) {
  //is the new member's fitness better than the best fitness?
  if (newMember->dFitness > spec->dBestFitness) {
    freeGenome(spec->sLeader);
    spec->dBestFitness = newMember->dFitness;
    spec->iGensNoImprovement = 0;
    spec->sLeader = copyGenome(newMember);
  }
  spec->vMembers[spec->iNumMembers++] = newMember;
}

/*  clears out all the members from the last generation, updates the age 
 * and gens no improvement.
 */
void purgeSpecies(sSpecies * spec) {
  spec->iNumMembers = 0;
  //update age etc
  spec->iAge++;
  spec->iGensNoImprovement++;
  spec->dSpawnsRqd = 0;
}

/*******************************************************************************
 * dump species - convenient for debug
 ******************************************************************************/

void dumpSpecies(FILE * out, sPopulation * pop) {
  int i;
  fprintf(out, "list of species (%d - CompatibilityThreshold=%f) :\n",
          pop->iNumSpecies, pop->sParams->dCompatibilityThreshold);
  for (i = 0; i < pop->iNumSpecies; i++)
    fprintf(out, "species %-2d - age %-3d - no improvement since %-2d "
            "generations - %-2d members - best fitness : %f - spawn : %6.2f - "
            "leader : genome %d\n",
            pop->vSpecies[i]->iSpeciesId, pop->vSpecies[i]->iAge,
            pop->vSpecies[i]->iGensNoImprovement,
            pop->vSpecies[i]->iNumMembers, pop->vSpecies[i]->dBestFitness,
            pop->vSpecies[i]->dSpawnsRqd, pop->vSpecies[i]->sLeader->iId);
}

/*******************************************************************************
 * speciation and spawn level
 ******************************************************************************/

/* add one species
 * chesks if the vector is big enougth, and extend it if necessary
 */
void addOneSpecies(sPopulation * pop, sGenome * firstOrg) {
  sSpecies * newSpec = createSpecies( firstOrg, pop->iNextSpeciesId++,
                                     pop->sParams->iNumIndividuals );
  if (pop->iNumSpecies == pop->iTotalSpecies) {
    pop->iTotalSpecies += INIT_VECT_SIZE;
    pop->vSpecies = realloc(pop->vSpecies,
                            pop->iTotalSpecies * sizeof(*pop->vSpecies));
  }
  pop->vSpecies[pop->iNumSpecies++] = newSpec;
}

/* remove one species
 * replace the empty slot by the last species
 */
void removeOneSpecies(sPopulation * pop, int id) {
  free(pop->vSpecies[id]->vMembers);
  freeGenome(pop->vSpecies[id]->sLeader);
  free(pop->vSpecies[id]);
  pop->iNumSpecies--;
  pop->vSpecies[id] = pop->vSpecies[pop->iNumSpecies];
}

/* separates each individual into its respective species by calculating
 * a compatibility score with every other member of the population and
 * niching accordingly. The function then adjusts the fitness scores of
 * each individual by species age and by sharing and also determines
 * how many offspring each individual should spawn.
 */
void speciateAndCalculateSpawnLevels(sPopulation * pop) {
  int gen, spec;
  bool bAdded;

  // iterate through each genome and speciate
  for (gen = 0; gen < pop->iNumGenomes; gen++) {
    bAdded = E_FALSE;
    // calculate its compatibility score with each species leader
    // if compatible add to species. If not, create a new species
    for (spec = 0; spec < pop->iNumSpecies; spec++) {
      double compatibility = getCompatibilityScore(pop->vGenomes[gen],
                                                   pop->vSpecies[spec]->sLeader, pop->sParams);
      // if this individual is similar to this species add to species
      if (compatibility <= pop->sParams->dCompatibilityThreshold) {
        addMember(pop->vSpecies[spec], pop->vGenomes[gen]);
        bAdded = E_TRUE;
        break;
      }
    }
    // if we have not found a compatible species, let's create a new one
    if (!bAdded)
      addOneSpecies(pop, pop->vGenomes[gen]);
  }

  // now all the genomes have been assigned a species the fitness scores
  // need to be adjusted to take into account sharing and species age
  for (spec = 0; spec < pop->iNumSpecies; spec++)
    adjustFitnesses(pop->vSpecies[spec], pop->sParams);

  // calculate new adjusted total & average fitness for the population
  for (gen = 0; gen < pop->iNumGenomes; gen++)
    pop->dTotFitAdj += pop->vGenomes[gen]->dAjustedFitness;

  pop->dAvFitAdj = pop->dTotFitAdj / pop->iNumGenomes;

  // calculate how many offspring each member of the population should spawn
  for (gen = 0; gen < pop->iNumGenomes; gen++) {
    double toSpawn = pop->vGenomes[gen]->dAjustedFitness / pop->dAvFitAdj;
    pop->vGenomes[gen]->dAmountToSpawn = toSpawn;
  }

  // iterate through all the species and calculate how many offspring
  // each species should spawn
  for (spec = 0; spec < pop->iNumSpecies; spec++)
    speciesSpawnAmount(pop->vSpecies[spec]);
}

/*******************************************************************************
 * tools to manipulate the species
 ******************************************************************************/

/* adjusts the fitness of each individual by first examining the species age 
 * and penalising if old, boosting if young.
 * performs fitness sharing by dividing the fitness by the number of individuals
 * in the species. This ensures a species does not grow too large
 */
void adjustFitnesses(sSpecies * spec, sParams * param) {
  int gen;
  for (gen = 0; gen < spec->iNumMembers; ++gen) {
    double fitness = spec->vMembers[gen]->dFitness;

    // boost the fitness scores if the species is young
    if (spec->iAge < param->iYoungBonusAgeThreshhold)
      fitness *= param->dYoungFitnessBonus;

    // punish older species
    if (spec->iAge > param->iOldAgeThreshold)
      fitness *= param->dOldAgePenalty;

    // apply fitness sharing to adjusted fitnesses
    spec->vMembers[gen]->dAjustedFitness = fitness / spec->iNumMembers;
  }
}

/* adds up the expected spawn amount for each individual in the species to 
 * calculate the amount of offspring this species should spawn
 */
void speciesSpawnAmount(sSpecies * spec) {
  int i;
  for (i = 0; i < spec->iNumMembers; i++)
    spec->dSpawnsRqd += spec->vMembers[i]->dAmountToSpawn;
}

/* returns a random genome selected from the best individuals of the species
 */
sGenome * randomAmongBest(sSpecies * spec, double dSurvivalRate) {
  sGenome * baby;
  if (spec->iNumMembers == 1) baby = spec->vMembers[0];
  else {
    int maxIndexSize = (int) (dSurvivalRate * spec->iNumMembers) + 1;
    int theOne = randInt(0, maxIndexSize);
    baby = spec->vMembers[theOne];
  }
  return baby;
}
