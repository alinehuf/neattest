//
//  data.h
//  gaGGPfinalState
//
//  Created by dexter on 05/04/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef gaGGPfinalState_loaddata_h
#define gaGGPfinalState_loaddata_h

#define MAX_LINE_SIZE 2048
#define INIT_DATA_VECT_SIZE 1000
#define MAX_PLAYERS_NUMBER 10

typedef struct {
  int * vFactIdxs;
  int * vScore;
  int iNumFacts;
} finalState;

typedef struct {
  int iNumPlayers;
  char ** vAllFacts;
  int iNumFacts;
  int iMaxFacts;
  finalState ** vAllFinalStates;
  int iNumFinalStates;
  int iMaxFinalStates;
} testBase;

typedef struct {
  double ** mInputs;
  double ** mOutputs;
  int iNumEntries;
  int iInputSize;
  int iOutputSize;
} NeatDataBase;

/* prototypes */

testBase * loadData(const char * filename);
void freeData(testBase * data);
int readGroup(char * text, int textIdx, char * buffer);
int * readScore(char * scoreBuffer, int * nbPlayers);
int getNumberFacts(char * buffer);
int getFactIdx(char * buffer, testBase * data);
int addFact(char * buffer, testBase * data);
void addFinalState(finalState * fstate, testBase * data);
void dumpFacts(testBase * data);
void dumpFinalStates(testBase * data);

NeatDataBase * testToNeatDataBase(testBase * dataBase);
void freeSimpleData(NeatDataBase * simpleData);
void dumpSimpleData(NeatDataBase * simpleData);

#endif
