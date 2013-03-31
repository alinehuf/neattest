//
//  phenotype.h
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_phenotype_h
#define NeatTest_0_1_phenotype_h

#include "global.h"

typedef struct neuron sNeuron;

/*******************************************************************************
 * link of the final ANN (phenotype)
 ******************************************************************************/

typedef struct {
  //pointers to the neurons this link connects
  sNeuron * pNeuronIn;
  sNeuron * pNeuronOut;
  //the connection weight
  double dWeight;
  //is this link a recurrent link?
  bool bRecurrent;
} sLink;

/*******************************************************************************
 * neuron of the final ANN (phenotype)
 ******************************************************************************/

struct neuron {
  int iID;                    // its identification number
  neuron_type eNeuronType;    // what type of neuron is this?
  double dActivationResponse; // sets the curvature of the sigmoid function
  sLink * vLinksIn;           // all the links coming into this neuron
  int iNumLinksIn;
  int iTotalLinksIn;
  sLink * vLinksOut;          // and out
  int iNumLinksOut;
  int iTotalLinksOut;
  double dSumActivation;      // sum of weights x inputs
  double dOutput;             // the output from this neuron
  //used in visualization of the phenotype
  int iPosX, iPosY;
  double dSplitX, dSplitY;
};

/*******************************************************************************
 * the final ANN (phenotype)
 ******************************************************************************/

typedef struct {
  sNeuron * vNeurons;
  int iNumNeurons;
  int iDepth; //the depth of the network
} sPhenotype;

/*******************************************************************************
 * prototypes
 ******************************************************************************/

sNeuron createNeuron(int id, neuron_type type, double x, double y, double act);
sLink createLink(double weight, sNeuron * in, sNeuron * out, bool rec);
void phenotypeAddLinkOut(sNeuron * neuron, sLink * tmpLink);
void phenotypeAddLinkIn(sNeuron * neuron, sLink * tmpLink);

#endif
