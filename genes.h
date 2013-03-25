//
//  genes.h
//  NeatTest_0.1
//
//  Created by dexter on 24/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_genes_h
#define NeatTest_0_1_genes_h

typedef enum {FALSE=0, TRUE} bool;

/*******************************************************************************
 * neuron gene
 ******************************************************************************/

// NONE type is used into Innovation table when a link is specified, neuron type
// is then not usefull
typedef enum {INPUT, HIDDEN, OUTPUT, BIAS, NONE} neuron_type;

// dSliptX, dSplitY : usefull for graphic representation of the genotype and
// also to determines if a link is reccurent
typedef struct neuronGene {
  int iId;
  neuron_type eNeuronType;
  bool bRecurrent;
  double dSigmoidCurvature; // activation rate AR  y = 1 / (1 + e^(-x/AR))
  double dSplitX, dSplitY;
} sNeuronGene;

/*******************************************************************************
 * link gene
 ******************************************************************************/

typedef struct linkGene {
  int iInnovId;
  int iFromNeuron;
  int iToNeuron;
  double dWeight;
  bool bEnabled;
  bool bRecurrent;
} sLinkGene;

#endif
