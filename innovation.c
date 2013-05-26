//
//  innovation.c
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "global.h"

// checks to see if this innovation has already occurred. If it has it
// returns the innovation ID. If not it returns a negative value.
int checkInnovation(sInnovTable * innovTable, int in, int out, innov_type type){
  int i;
  for (i = 0; i < innovTable->iNumInnovs; i++) {
    if ( innovTable->vInnovs[i]->iNeuronIn == in &&
         innovTable->vInnovs[i]->iNeuronOut == out &&
         innovTable->vInnovs[i]->eInnovType == type ) {
			// found a match so assign this innovation number to id
			return i;
		}
  }
  return -1;
}

int createNewInnov(sInnovTable * innovTable, int in, int out, innov_type type,
                   neuron_type ntype) {
  sInnovation * innov = malloc(sizeof(*innov));
  innov->eInnovType = type;
  innov->iInnovId = innovTable->iNumInnovs;
  innov->iNeuronIn = in;
  innov->iNeuronOut = out;
  innov->iNeuronId = 0;
  innov->eNeuronType = ntype;

  if (type == NEW_NEURON)
    innov->iNeuronId = innovTable->iNextNeuronId++;

  if (innovTable->iNumInnovs == innovTable->iTotalInnovs) {
    innovTable->iTotalInnovs += INIT_VECT_SIZE;
    innovTable->vInnovs = realloc(innovTable->vInnovs,
                     innovTable->iTotalInnovs * sizeof(*(innovTable->vInnovs)));
  }
  innovTable->vInnovs[innovTable->iNumInnovs++] = innov;
  return innov->iNeuronId;
}

sInnovTable * createNewInnovTable( sNeuronGene ** neurons, int numNeurons,
                                   sLinkGene ** links, int numLinks ) {
  sInnovTable * innovTable = malloc(sizeof(*innovTable));
  innovTable->iNumInnovs = 0;
  innovTable->iNextNeuronId = 0;
  innovTable->iTotalInnovs = INIT_VECT_SIZE;
  // allocate memory for vectors of innovations
  innovTable->vInnovs = calloc( innovTable->iTotalInnovs,
                               sizeof(*(innovTable->vInnovs)) );
  int i;
  // add the neurons
  for (i = 0; i < numNeurons; i++)
    createNewInnov(innovTable, -1, -1, NEW_NEURON, neurons[i]->eNeuronType);
  // add the links
  for (i = 0; i < numLinks; i++)
    createNewInnov(innovTable, links[i]->iFromNeuron, links[i]->iToNeuron,
                   NEW_LINK, NONE);
  return innovTable;
}

void freeInnovTable(sInnovTable * innovTable) {
  int i;
  for (i = 0; i < innovTable->iNumInnovs; i++)
    free(innovTable->vInnovs[i]);
  free(innovTable->vInnovs);
  free(innovTable);
}

/*******************************************************************************
 * convenient for debug - visualization
 ******************************************************************************/

void dumpInnovTable(FILE * out, sInnovTable * innovTable) {
  fprintf(out, "innovTable (%d innovations, next neuron ID = %d) :\n",
         innovTable->iNumInnovs, innovTable->iNextNeuronId);
  int i;
  for (i = 0; i < innovTable->iNumInnovs; i++) {
    fprintf(out, "innov %-2d ", innovTable->vInnovs[i]->iInnovId);
    switch(innovTable->vInnovs[i]->eInnovType) {
      case NEW_NEURON:
        fprintf(out, "NEW_NEURON ");
        break;
      case NEW_LINK:
        fprintf(out, "NEW_LINK   ");
        break;
      default:
        fprintf(out, "UNKNOWN_TYPE ");
        break;
    }
    fprintf(out, "in=%-2d out=%-2d neuronID=%-2d neuronType=",
           innovTable->vInnovs[i]->iNeuronIn,
           innovTable->vInnovs[i]->iNeuronOut,
           innovTable->vInnovs[i]->iNeuronId);
    switch(innovTable->vInnovs[i]->eNeuronType) {
      case INPUT:
        fprintf(out, "INPUT\n");
        break;
      case HIDDEN:
        fprintf(out, "HIDDEN\n");
        break;
      case OUTPUT:
        fprintf(out, "OUTPUT\n");
        break;
      case BIAS:
        fprintf(out, "BIAS\n");
        break;
      case NONE:
        fprintf(out, "NONE\n");
        break;
      default:
        fprintf(out, "UNKNOWN_TYPE\n");
        break;
    }
  }
}
