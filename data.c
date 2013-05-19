//
//  data.c
//  gaGGPfinalState
//
//  Created by dexter on 05/04/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

/* load the data, 
 * on each line a list of facts describing the final state of a game
 * the first group (...) gives the scores
 */
testBase * loadData(const char * filename) {
  FILE * fileref;
  if ((fileref = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "error : unable to open dump file %s", filename);
    exit(EXIT_FAILURE);
  }
  // test database structure
  // calloc: to initialize all other values ​​to 0
  testBase * data = (testBase *) calloc(1, sizeof(*data));
  data->iMaxFacts = INIT_DATA_VECT_SIZE;
  data->vAllFacts = (char* *)
                       malloc(INIT_DATA_VECT_SIZE * sizeof(*(data->vAllFacts)));
  data->iMaxFinalStates = INIT_DATA_VECT_SIZE;
  data->vAllFinalStates = (finalState* *)
                 malloc(INIT_DATA_VECT_SIZE * sizeof(*(data->vAllFinalStates)));
  int i;
  int idxFact = -1;
  char lineBuffer[BUFSIZ];            // a line of the data file
  int charIdx;                        // id of the current char in the line
  char scoreBuffer[BUFSIZ];           // buffer for the string of scores
  char buffer[BUFSIZ];                // general buffer (for differents facts)

  while(fgets(lineBuffer, BUFSIZ, fileref) != NULL) {
    // new final state
    finalState * fstate = malloc(sizeof(*fstate));

    // read score and player strings
    charIdx = readGroup(lineBuffer, 0, scoreBuffer);      // score

    // get the score
    fstate->vScore = readScore(scoreBuffer, &data->iNumPlayers);

    // get the facts => final state of the game
    // get the number of facts
    fstate->iNumFacts = getNumberFacts(&lineBuffer[charIdx]);
    // allocate memory
    fstate->vFactIdxs =
               (int *) malloc(fstate->iNumFacts * sizeof(*(fstate->vFactIdxs)));
    i = 0;
    while ((charIdx = readGroup(lineBuffer, charIdx, buffer)) != -1) {
      // fact already in the base
      if ((idxFact = getFactIdx(buffer, data)) != -1)
        fstate->vFactIdxs[i++] = idxFact;
      // else add it
      else fstate->vFactIdxs[i++] = addFact(buffer, data);
    }
    addFinalState(fstate, data);
  }
  fclose(fileref);
  return data;
}

/* free the memory allocated for the data : the final states of the game
 */
void freeData(testBase * data) {
  int i;
  for (i = 0; i < data->iNumFacts; i++)
    free(data->vAllFacts[i]);
  free(data->vAllFacts);
  for (i = 0; i < data->iNumFinalStates; i++) {
    free(data->vAllFinalStates[i]->vScore);
    free(data->vAllFinalStates[i]->vFactIdxs);
    free(data->vAllFinalStates[i]);
  }
  free(data->vAllFinalStates);
  free(data);
}

/* read a fact in the form of (...) 
 */
int readGroup(char * text, int textIdx, char * buffer) {
  int bufferIdx = 0;
  // buffer is NULL
  if (buffer == 0) {
    fprintf(stderr, "error in readGroup() : buffer is NULL");
    exit(EXIT_FAILURE);
  }
  // text is NULL
  if (text == NULL || text[textIdx] == '\0') {
    buffer[0] = '\0';
    return -1;
  }
  // find the opening parenthesis
  while (text[textIdx] != '\0' && text[textIdx] != '(')
    textIdx++;
  while (text[textIdx] != '\0' && text[textIdx] == '(')
    textIdx++;
  // copy the text until the closing parenthesis
  while (text[textIdx] != '\0' && text[textIdx] != ')')
    buffer[bufferIdx++] = text[textIdx++];
  // nothing in the buffer
  if (bufferIdx == 0) return -1;
  // else add end string char
  buffer[bufferIdx] = '\0';
  return textIdx;
}

/* get the scores into the final state
 */
int * readScore(char * scoreBuffer, int * iNumPlayers) {
  int i;
  int intBuffer[MAX_PLAYERS_NUMBER];
  int numScores = 0;
  char * number;

  while((number = strsep(&scoreBuffer, " \t")) != NULL) {
    if ((i = sscanf(number, "%d", &intBuffer[numScores++])) < 1) {
      numScores--;
      break;
    }
  }
  if (*iNumPlayers != 0 && numScores != *iNumPlayers) {
    fprintf(stderr, "error in readScore(): different number of scores");
    exit(EXIT_FAILURE);
  }
  *iNumPlayers = numScores;
  if (*iNumPlayers == 0)
    fprintf(stderr, "error in readScore(): no score found");
  int * scores = malloc(numScores * sizeof(*scores));
  memcpy(scores, intBuffer, numScores * sizeof(*scores));
  return scores;
}

/* count the number of facts into the final state 
 */
int getNumberFacts(char * buffer) {
  int i = 0;
  int nb = 0;
  while (buffer[i] != 0)
    if (buffer[i++] == '(') nb++;
  return nb;
}

/* check if the fact already exist and return its index in the vector of all
 * facts or -1
 */
int getFactIdx(char * buffer, testBase * data) {
  int i;
  for (i = 0; i < data->iNumFacts; i++)
    if (strcmp(buffer, data->vAllFacts[i]) == 0)
      return i;
  return -1;
}

/* add a fact into the list of the all differents facts
 * check if the vector is big enougth and reallocate if necessary
 */
int  addFact(char * buffer, testBase * data) {
  if (data->iNumFacts == data->iMaxFacts) {
    data->iMaxFacts += INIT_DATA_VECT_SIZE;
    data->vAllFacts = (char* *) realloc(data->vAllFacts,
                                  data->iMaxFacts * sizeof(*(data->vAllFacts)));
  }
  data->vAllFacts[data->iNumFacts] = strdup(buffer);
  data->iNumFacts++;
  return data->iNumFacts - 1;
}

/* add a new final state into the vector
 * check if the vector is big enougth and reallocate if necessary
 */
void addFinalState(finalState * fstate, testBase * data) {
  if (data->iNumFinalStates == data->iMaxFinalStates) {
    data->iMaxFinalStates += INIT_DATA_VECT_SIZE;
    data->vAllFinalStates = (finalState* *) realloc(data->vAllFinalStates,
                      data->iMaxFinalStates * sizeof(*(data->vAllFinalStates)));
  }
  data->vAllFinalStates[data->iNumFinalStates++] = fstate;
}

/* convenient for debug
 * write the list of all differents facts found into the final states
 */
void dumpFacts(testBase * data) {
  int i;
  printf("%d differents facts :\n", data->iNumFacts);
  for (i = 0; i < data->iNumFacts; i++)
    printf("%d - (%s)\n", i, data->vAllFacts[i]);
  printf("\n");
}

/* convenient for debug
 * write the list of final states with the score
 */
void dumpFinalStates(testBase * data) {
  int i, j;
  puts("All final states :");
  for (i = 0; i < data->iNumFinalStates; i++) {
    for (j = 0; j < data->iNumPlayers; j++)
      printf("%d ", data->vAllFinalStates[i]->vScore[j]);
    printf("=> ");
    for (j = 0; j < data->vAllFinalStates[i]->iNumFacts; j++)
      printf("(%s) ", data->vAllFacts[data->vAllFinalStates[i]->vFactIdxs[j]]);
    printf("\n");
  }
  printf("\n");
}

/*******************************************************************************
 * a simple data base in the form of 2D tables of inputs and outputs 
 * to feed NEAT
 ******************************************************************************/

NeatDataBase * testToNeatDataBase(testBase * data) {
  NeatDataBase * simpleData = malloc(sizeof(*simpleData));
  // number of entries, sumber of input values and output values
  simpleData->iNumEntries = data->iNumFinalStates;
  simpleData->iInputSize = data->iNumFacts;
  simpleData->iOutputSize = data->iNumPlayers;
  // allocate memory for each input/output entry
  simpleData->mInputs =
                 malloc(simpleData->iNumEntries * sizeof(*simpleData->mInputs));
  simpleData->mOutputs =
                malloc(simpleData->iNumEntries * sizeof(*simpleData->mOutputs));
  int i, j;
  for (i = 0; i < simpleData->iNumEntries; i++) {
    // allocate memory for the input values and the output values
    simpleData->mInputs[i] =
                calloc(simpleData->iInputSize, sizeof(*simpleData->mInputs[i]));
    simpleData->mOutputs[i] =
              calloc(simpleData->iOutputSize, sizeof(*simpleData->mOutputs[i]));
    // for each fact present in the final state, input is 1, else 0
    for (j = 0; j < data->vAllFinalStates[i]->iNumFacts; j++)
      simpleData->mInputs[i][data->vAllFinalStates[i]->vFactIdxs[j]] = 1;
    // convert score [0, 100] to output [0, 1]
    for (j = 0; j < simpleData->iOutputSize; j++)
      simpleData->mOutputs[i][j] =
                              (double) data->vAllFinalStates[i]->vScore[j] / 100;
  }
  return simpleData;
}

void freeSimpleData(NeatDataBase * simpleData) {
  int i;
  for (i = 0; i < simpleData->iNumEntries; i++) {
    free(simpleData->mInputs[i]);
    free(simpleData->mOutputs[i]);
  }
  free(simpleData->mInputs);
  free(simpleData->mOutputs);
  free(simpleData);
}

void dumpSimpleData(NeatDataBase * simpleData) {
  puts("All final states in simple form :");
  int i, j;
  for (i = 0; i < simpleData->iNumEntries; i++) {
    for (j = 0; j < simpleData->iOutputSize; j++)
      printf("%4.2f ", simpleData->mOutputs[i][j]);
    printf("=> ");
    for (j = 0; j < simpleData->iInputSize; j++)
      printf("%d ", (int) simpleData->mInputs[i][j]);
    printf("\n");
  }
}

