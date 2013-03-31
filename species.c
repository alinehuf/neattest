//
//  species.c
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "species.h"
#include "utils.h"

/*******************************************************************************
 * gestion of species
 ******************************************************************************/

/* creates an instance of a new species. A local copy of the initializing genome
 * is kept in sLeader and the first element of vMembers is a pointer to that 
 * genome.
 */
sSpecies * createSpecies(sGenome * firstOrg,int speciesId, int iNumIndividuals){
  sSpecies * spec = (sSpecies *) malloc(sizeof(*spec));
  spec->sLeader = *firstOrg;
  spec->vMembers = (sGenome* *) calloc(iNumIndividuals,sizeof(*spec->vMembers));
  spec->vMembers[0] = firstOrg;
  spec->iNumMembers = 1;
  spec->iTotalMembers = iNumIndividuals;
  spec->iSpeciesId = speciesId;
  spec->dBestFitness = firstOrg->dFitness;
  spec->iGensNoImprovement = 0;
  spec->iAge = 0;
  spec->dSpawnsRqd = 0;
  return spec;
}

/* adds a new member to this species and updates the member variables
 * accordingly
 */
void addMember(sSpecies * spec, sGenome * newMember) {
  //is the new member's fitness better than the best fitness?
  if (newMember->dFitness > spec->dBestFitness) {
    spec->dBestFitness = newMember->dFitness;
    spec->iGensNoImprovement = 0;
    spec->sLeader = *newMember;
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

/* adjusts the fitness of each individual by first examining the species age 
 * and penalising if old, boosting if young.
 * performs fitness sharing by dividing the fitness by the number of individuals
 * in the species. This ensures a species does not grow too large
 */
void adjustFitnesses(sSpecies * spec, sParams * param) {
  int gen;
  for (gen = 0; gen < spec->iNumMembers; ++gen) {
    double fitness = spec->vMembers[gen]->dFitness;

    //boost the fitness scores if the species is young
    if (spec->iAge < param->iYoungBonusAgeThreshhold)
      fitness *= param->dYoungFitnessBonus;

    //punish older species
    if (spec->iAge > param->iOldAgeThreshold)
      fitness *= param->dOldAgePenalty;

    //apply fitness sharing to adjusted fitnesses
    spec->vMembers[gen]->dAjustedFitness = fitness / spec->iNumMembers;;
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

/* returns a random genome selected from the best individuals
 */
sGenome randomSpawn(sSpecies * spec, double dSurvivalRate) {
  sGenome baby;
  if (spec->iNumMembers == 1) baby = *spec->vMembers[0];
  else {
    int maxIndexSize = (int) (dSurvivalRate * spec->iNumMembers) + 1;
    int theOne = randInt(0, maxIndexSize);
    baby = *spec->vMembers[theOne];
  }
  return baby;
}

/*******************************************************************************
 * manipulation of doubly linked list of species
 ******************************************************************************/

listSpecies * addOneSpecies(listSpecies * lastSpec, sGenome * firstOrg,
                                           int speciesId, int iNumIndividuals) {
  sSpecies * newSpec = createSpecies(firstOrg, speciesId, iNumIndividuals);
  listSpecies * newSlot = (listSpecies *) malloc(sizeof(*newSlot));
  newSlot->sSpecies = newSpec;
  newSlot->prev = lastSpec;
  newSlot->next = NULL;
  lastSpec->next = newSlot;
  return newSlot;
}

listSpecies * removeOneSpecies(listSpecies * curSpec) {
  curSpec->prev->next = curSpec->next;
  curSpec->next ->prev= curSpec->prev;
  listSpecies * next = curSpec->next;
  free(curSpec->sSpecies->vMembers);
  free(curSpec->sSpecies);
  free(curSpec);
  return next;
}
