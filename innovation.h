//
//  innovation.h
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_innovation_h
#define NeatTest_0_1_innovation_h

#include "genes.h"

/*******************************************************************************
 * innovation
 ******************************************************************************/

typedef enum {NEW_NEURON, NEW_LINK} innov_type;

typedef struct {
  innov_type eInnovType;
  int iInnovId;
  int iNeuronIn;
  int iNeuronOut;
  int iNeuronId;
  neuron_type eNeuronType;
} sInnovation;

/*******************************************************************************
 * table of innovations
 ******************************************************************************/

typedef struct {
  sInnovation * vInnovs;
  int iNumInnovs;
  int iTotalInnovs;
  int iNextNeuronId;
} sInnovTable;

/*******************************************************************************
 * prototypes
 ******************************************************************************/

int checkInnovation(sInnovTable * innovTable, int in, int out, innov_type type);
int createNewInnov(sInnovTable * innovTable, int in, int out, innov_type type,
                   neuron_type ntype);
sInnovTable createNewInnovTable(sNeuronGene * neurons, int numGenes,
                                sLinkGene * links, int numLinks);
void dumpInnovTable(sInnovTable innovTable);

#endif
