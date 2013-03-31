//
//  params.c
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"

#define DEFAULT_CONF_FILE "params.ini"
enum {DOUBLE_TYPE = 0, FLOAT_TYPE, INT_TYPE};

/*******************************************************************************
 * LOADING THE CONFIGURATION FILE
 ******************************************************************************/

/* load configuration from a file
 * allocates , memory must be freed after with freeData()
 */
sParams loadConf(char * file) {
  if (! file) file = DEFAULT_CONF_FILE;
  // open file
  FILE * fileref = fopen(file, "r");
  if (! fileref)
    error("Error in loadConf(): Unable to open file %s.\n", file);

  sParams params;
  loadNumber("dWeightMutationRate", &params.dWeightMutationRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dMaxWeightPerturbation", &params.dMaxWeightPerturbation,
             DOUBLE_TYPE, fileref);  
  loadNumber("dProbabilityWeightReplaced", &params.dProbabilityWeightReplaced,
             DOUBLE_TYPE, fileref);
  loadNumber("dActivationMutationRate", &params.dActivationMutationRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dMaxActivationPerturbation", &params.dMaxActivationPerturbation,
             DOUBLE_TYPE, fileref);
  loadNumber("dChanceAddLink", &params.dChanceAddLink,
             DOUBLE_TYPE, fileref);
  loadNumber("dChanceAddRecurrentLink", &params.dChanceAddRecurrentLink,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindLoop", &params.iNumTrysToFindLoop,
             INT_TYPE, fileref);
  loadNumber("iNumTrysToAddLink", &params.iNumTrysToAddLink,
             INT_TYPE, fileref);
  loadNumber("dChanceAddNode", &params.dChanceAddNode,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindOldLink", &params.iNumTrysToFindOldLink,
             INT_TYPE, fileref);
  loadNumber("dCompatibilityThreshold", &params.dCompatibilityThreshold,
             DOUBLE_TYPE, fileref);
  loadNumber("dExcessGenesCoef", &params.dExcessGenesCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("dDisjointGenesCoef", &params.dDisjointGenesCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("dWeightDiffCoef", &params.dWeightDiffCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("iMaxSpecies", &params.iMaxSpecies,
             INT_TYPE, fileref);
  loadNumber("iYoungBonusAgeThreshhold", &params.iYoungBonusAgeThreshhold,
             INT_TYPE, fileref);
  loadNumber("dYoungFitnessBonus", &params.dYoungFitnessBonus,
             DOUBLE_TYPE, fileref);
  loadNumber("iOldAgeThreshold", &params.iOldAgeThreshold,
             INT_TYPE, fileref);
  loadNumber("dOldAgePenalty", &params.dOldAgePenalty,
             DOUBLE_TYPE, fileref);
  loadNumber("dSurvivalRate", &params.dSurvivalRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dCrossoverRate", &params.dCrossoverRate,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindMate", &params.iNumTrysToFindMate,
             INT_TYPE, fileref);
  loadNumber("iNumGensAllowedNoImprov", &params.iNumGensAllowedNoImprov,
             INT_TYPE, fileref);
  loadNumber("iNumIndividuals", &params.iNumIndividuals,
             INT_TYPE, fileref);
  loadNumber("iNumEpoch", &params.iNumEpoch,
             INT_TYPE, fileref);
  
  fclose(fileref);

  // dump
  puts("------------ parameters :");
  printf("dWeightMutationRate :        %f\n", params.dWeightMutationRate);
  printf("dMaxWeightPerturbation :     %f\n", params.dMaxWeightPerturbation);
  printf("dProbabilityWeightReplaced : %f\n",params.dProbabilityWeightReplaced);
  printf("dActivationMutationRate :    %f\n", params.dActivationMutationRate);
  printf("dMaxActivationPerturbation : %f\n",params.dMaxActivationPerturbation);
  printf("dChanceAddLink :             %f\n", params.dChanceAddLink);
  printf("dChanceAddRecurrentLink :    %f\n", params.dChanceAddRecurrentLink);
  printf("iNumTrysToFindLoop :         %d\n", params.iNumTrysToFindLoop);
  printf("iNumTrysToAddLink :          %d\n", params.iNumTrysToAddLink);
  printf("dChanceAddNode :             %f\n", params.dChanceAddNode);
  printf("iNumTrysToFindOldLink :      %d\n", params.iNumTrysToFindOldLink);
  printf("dCompatibilityThreshold :    %f\n", params.dCompatibilityThreshold);
  printf("dExcessGenesCoef :           %f\n", params.dExcessGenesCoef);
  printf("dDisjointGenesCoef :         %f\n", params.dDisjointGenesCoef);
  printf("dWeightDiffCoef :            %f\n", params.dWeightDiffCoef);
  printf("iMaxSpecies :                %d\n", params.iMaxSpecies);
  printf("iYoungBonusAgeThreshhold     %d\n", params.iYoungBonusAgeThreshhold);
  printf("dYoungFitnessBonus :         %f\n", params.dYoungFitnessBonus);
  printf("iOldAgeThreshold :           %d\n", params.iOldAgeThreshold);
  printf("dOldAgePenalty :             %f\n", params.dOldAgePenalty);
  printf("dSurvivalRate :              %f\n", params.dSurvivalRate);
  printf("dCrossoverRate :             %f\n", params.dCrossoverRate);
  printf("iNumTrysToFindMate :         %d\n", params.iNumTrysToFindMate);
  printf("iNumGensAllowedNoImprov :    %d\n", params.iNumGensAllowedNoImprov);
  printf("iNumIndividuals :            %d\n", params.iNumIndividuals);
  printf("iNumEpoch :                  %d\n", params.iNumEpoch);
  puts("-------------------------");

  return params;
}

/* convenient for error messages
 */
void error(char * format, char * message) {
  if(!format) format = "%s";
  fprintf(stderr, format, message); fprintf(stderr, "\n");
  exit(1);
}

/* load a string from the file and throw away the comments
 */
char * loadString(char * prefix, FILE * fileref) {
  char buffer[BUFSIZ] = {0};       // to get the comments
  char format[BUFSIZ] = {0};       // for the format string
  char * var;                      // for the data
  if ((var = calloc(BUFSIZ, sizeof(char))) == NULL) // allocated memory
    error(0, "Error : calloc faillure");
  strcat(format, prefix);          // prefix to be sure to take the good value
  strcat(format, " %s");
  if(fscanf(fileref, format, var) < 1)
    error("Error : unable to find %s in the configuration file", prefix);
  if (fgets(buffer, BUFSIZ, fileref) == 0) // to throw away the comments
    error(0, "Error : while throwing away comments");
  return var;
}

/* load an integer from the file and throw away the comments
 */
void loadNumber(char * prefix, void * var, int type, FILE * fileref) {
  char buffer[BUFSIZ] = {0};       // to get the comments
  char format[BUFSIZ] = {0};       // for the format string
  strcat(format, prefix);
  if (type == DOUBLE_TYPE) strcat(format, " %lf");
  else if (type == FLOAT_TYPE) strcat(format, " %f");
  else strcat(format, " %d");
  if(fscanf(fileref, format, var) < 1)
    error("Error : unable to find %s in the configuration file", prefix);
  if (fgets(buffer, BUFSIZ, fileref) == 0) // to throw away the comments
    error(0, "Error : while throwing away comments");
}