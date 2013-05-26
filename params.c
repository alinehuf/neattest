//
//  params->c
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"

#define DEFAULT_CONF_FILE "/Users/dex/FAC/M1_MEMOIRE/NeatTest_0.1/NeatTest_0.1/params/params.ini"
enum {DOUBLE_TYPE = 0, FLOAT_TYPE, INT_TYPE};

/*******************************************************************************
 * LOADING THE CONFIGURATION FILE
 ******************************************************************************/

/* load configuration from a file
 */
sParams * loadConf(char * file) {
  if (! file) file = DEFAULT_CONF_FILE;
  // open file
  FILE * fileref = fopen(file, "r");
  if (! fileref)
    error("Error in loadConf(): Unable to open file %s.\n", file);

  sParams * params = calloc(1, sizeof(*params));
  
  loadNumber("dWeightMutationRate", &params->dWeightMutationRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dMaxWeightPerturbation", &params->dMaxWeightPerturbation,
             DOUBLE_TYPE, fileref);  
  loadNumber("dProbabilityWeightReplaced", &params->dProbabilityWeightReplaced,
             DOUBLE_TYPE, fileref);
  loadNumber("dActivationMutationRate", &params->dActivationMutationRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dMaxActivationPerturbation", &params->dMaxActivationPerturbation,
             DOUBLE_TYPE, fileref);
  loadNumber("dChanceAddLink", &params->dChanceAddLink,
             DOUBLE_TYPE, fileref);
  loadNumber("dChanceAddRecurrentLink", &params->dChanceAddRecurrentLink,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindLoop", &params->iNumTrysToFindLoop,
             INT_TYPE, fileref);
  loadNumber("iNumTrysToAddLink", &params->iNumTrysToAddLink,
             INT_TYPE, fileref);
  loadNumber("dChanceAddNode", &params->dChanceAddNode,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindOldLink", &params->iNumTrysToFindOldLink,
             INT_TYPE, fileref);
  loadNumber("dCompatibilityThreshold", &params->dCompatibilityThreshold,
             DOUBLE_TYPE, fileref);
  loadNumber("dExcessGenesCoef", &params->dExcessGenesCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("dDisjointGenesCoef", &params->dDisjointGenesCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("dWeightDiffCoef", &params->dWeightDiffCoef,
             DOUBLE_TYPE, fileref);
  loadNumber("iMaxSpecies", &params->iMaxSpecies,
             INT_TYPE, fileref);
  loadNumber("iYoungBonusAgeThreshhold", &params->iYoungBonusAgeThreshhold,
             INT_TYPE, fileref);
  loadNumber("dYoungFitnessBonus", &params->dYoungFitnessBonus,
             DOUBLE_TYPE, fileref);
  loadNumber("iOldAgeThreshold", &params->iOldAgeThreshold,
             INT_TYPE, fileref);
  loadNumber("dOldAgePenalty", &params->dOldAgePenalty,
             DOUBLE_TYPE, fileref);
  loadNumber("dSurvivalRate", &params->dSurvivalRate,
             DOUBLE_TYPE, fileref);
  loadNumber("dCrossoverRate", &params->dCrossoverRate,
             DOUBLE_TYPE, fileref);
  loadNumber("iNumTrysToFindMate", &params->iNumTrysToFindMate,
             INT_TYPE, fileref);
  loadNumber("iNumGensAllowedNoImprov", &params->iNumGensAllowedNoImprov,
             INT_TYPE, fileref);
  loadNumber("iNumIndividuals", &params->iNumIndividuals,
             INT_TYPE, fileref);
  loadNumber("iNumEpoch", &params->iNumEpoch,
             INT_TYPE, fileref);
  loadNumber("iDumpEvery", &params->iDumpEvery,
             INT_TYPE, fileref);
  fclose(fileref);

  return params;
}

// dump
void dumpParams(FILE * out, sParams * params) {
  fprintf(out, "------------ parameters :\n");
  fprintf(out, "dWeightMutationRate :        %f\n", params->dWeightMutationRate);
  fprintf(out, "dMaxWeightPerturbation :     %f\n", params->dMaxWeightPerturbation);
  fprintf(out, "dProbabilityWeightReplaced : %f\n", params->dProbabilityWeightReplaced);
  fprintf(out, "dActivationMutationRate :    %f\n", params->dActivationMutationRate);
  fprintf(out, "dMaxActivationPerturbation : %f\n", params->dMaxActivationPerturbation);
  fprintf(out, "dChanceAddLink :             %f\n", params->dChanceAddLink);
  fprintf(out, "dChanceAddRecurrentLink :    %f\n", params->dChanceAddRecurrentLink);
  fprintf(out, "iNumTrysToFindLoop :         %d\n", params->iNumTrysToFindLoop);
  fprintf(out, "iNumTrysToAddLink :          %d\n", params->iNumTrysToAddLink);
  fprintf(out, "dChanceAddNode :             %f\n", params->dChanceAddNode);
  fprintf(out, "iNumTrysToFindOldLink :      %d\n", params->iNumTrysToFindOldLink);
  fprintf(out, "dCompatibilityThreshold :    %f\n", params->dCompatibilityThreshold);
  fprintf(out, "dExcessGenesCoef :           %f\n", params->dExcessGenesCoef);
  fprintf(out, "dDisjointGenesCoef :         %f\n", params->dDisjointGenesCoef);
  fprintf(out, "dWeightDiffCoef :            %f\n", params->dWeightDiffCoef);
  fprintf(out, "iMaxSpecies :                %d\n", params->iMaxSpecies);
  fprintf(out, "iYoungBonusAgeThreshhold     %d\n", params->iYoungBonusAgeThreshhold);
  fprintf(out, "dYoungFitnessBonus :         %f\n", params->dYoungFitnessBonus);
  fprintf(out, "iOldAgeThreshold :           %d\n", params->iOldAgeThreshold);
  fprintf(out, "dOldAgePenalty :             %f\n", params->dOldAgePenalty);
  fprintf(out, "dSurvivalRate :              %f\n", params->dSurvivalRate);
  fprintf(out, "dCrossoverRate :             %f\n", params->dCrossoverRate);
  fprintf(out, "iNumTrysToFindMate :         %d\n", params->iNumTrysToFindMate);
  fprintf(out, "iNumGensAllowedNoImprov :    %d\n", params->iNumGensAllowedNoImprov);
  fprintf(out, "iNumIndividuals :            %d\n", params->iNumIndividuals);
  fprintf(out, "iNumEpoch :                  %d\n", params->iNumEpoch);
  fprintf(out, "iDumpEvery :                 %d\n", params->iDumpEvery);
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