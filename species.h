//
//  species.h
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_species_h
#define NeatTest_0_1_species_h

#include "genome.h"

typedef struct {
  sGenome * sLeader;   // keep a local copy of the first member of this species
  sGenome ** vMembers; // pointers to all the genomes within this species
  int iNumMembers;
  int iTotalMembers;
  int iSpeciesId;      // the species needs an identification number
  double dBestFitness;    // best fitness found so far by this species
  int iGensNoImprovement; // generations since fitness has improved
  int iAge;
  double dSpawnsRqd; // how many of this species should be spawned for next pop
} sSpecies;

// doubly linked list of species
typedef struct listSpecies {
  sSpecies * sSpecies;
  struct listSpecies * cdr;
} listSpecies;

// species
sSpecies * createSpecies(sGenome * firstOrg, int specId, int iNumIndividuals);
void addMember(sSpecies * spec, sGenome * newMember);
void purgeSpecies(sSpecies * spec);
void adjustFitnesses(sSpecies * spec, sParams * param);
void speciesSpawnAmount(sSpecies * spec);
sGenome * randomAmongBest(sSpecies * spec, double dSurvivalRate);
// doubly linked list of species
listSpecies * addOneSpecies(listSpecies * lastSpec, sGenome * firstOrg,
                                            int speciesId, int iNumIndividuals);
listSpecies * removeOneSpecies(listSpecies ** listSpec, int specId);
// convenient for debug
void dumpSpecies(listSpecies * listSpec);

#endif
