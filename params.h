//
//  params.h
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_params_h
#define NeatTest_0_1_params_h

typedef struct {
  // weight mutation
  double dWeightMutationRate;
  double dMaxWeightPerturbation;
  double dProbabilityWeightReplaced;
  // activation (sigmoid curvature) mutation
  double dActivationMutationRate;
  double dMaxActivationPerturbation;
  // mutation by adding a link
  double dChanceAddLink;
  double dChanceAddRecurrentLink;
  int iNumTrysToFindLoop;
  int iNumTrysToAddLink;
  // mutation by adding a node
  double dChanceAddNode;
  int iNumTrysToFindOldLink;
  // speciation compatibility distance
  double dCompatibilityThreshold;
  double dExcessGenesCoef;
  double dDisjointGenesCoef;
  double dWeightDiffCoef;
  // speciation distribution
  int iMaxSpecies;
  int iYoungBonusAgeThreshhold;
  double dYoungFitnessBonus;
  int iOldAgeThreshold;
  double dOldAgePenalty;
  // spawn
  double dSurvivalRate;
  double dCrossoverRate;
  int iNumTrysToFindMate;
  // population
  int iNumGensAllowedNoImprov;
  //global
  int iNumIndividuals;
  int iNumEpoch;
} sParams;

/*******************************************************************************
 * prototypes
 ******************************************************************************/

sParams * loadConf(char * file);
void error(char * format, char * message);
char * loadString(char * prefix, FILE * fileref);
void loadNumber(char * prefix, void * var, int type, FILE * fileref);

#endif
