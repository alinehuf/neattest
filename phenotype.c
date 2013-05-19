//
//  phenotype.c
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "phenotype.h"

/*******************************************************************************
 * neuron of the final ANN (phenotype)
 ******************************************************************************/

sNeuron * createNeuron(int id, neuron_type type, double x, double y,
                       double act) {
  sNeuron * neuron = malloc(sizeof(*neuron));
  neuron->iID = id;
  neuron->eNeuronType = type;
  neuron->dActivationResponse = act;
  neuron->iNumLinksIn = 0;
  neuron->iTotalLinksIn = INIT_VECT_SIZE;
  neuron->vLinksIn = calloc(INIT_VECT_SIZE, sizeof(*neuron->vLinksIn));
  neuron->iNumLinksOut = 0;
  neuron->iTotalLinksOut = INIT_VECT_SIZE;
  neuron->vLinksOut = calloc(INIT_VECT_SIZE, sizeof(*neuron->vLinksOut));
  neuron->dSumActivation = 0;
  neuron->dOutput = 0;
  neuron->iPosX = 0;
  neuron->iPosY = 0;
  neuron->dSplitX = x;
  neuron->dSplitY = y;
  return neuron;
}

/*******************************************************************************
 * link of the final ANN (phenotype)
 ******************************************************************************/

sLink * createLink(double weight, sNeuron * in, sNeuron * out, bool rec) {
  sLink * link = malloc(sizeof(*link));
  link->pNeuronIn = in;
  link->pNeuronOut = out;
  link->dWeight = weight;
  link->bRecurrent = rec;
  return link;
}

/*******************************************************************************
 * the final ANN (phenotype)
 ******************************************************************************/

/* add a link to the list of outcomming links, reallocate if necessary */
void phenotypeAddLinkOut(sNeuron * neuron, sLink * tmpLink) {
  if (neuron->iNumLinksOut == neuron->iTotalLinksOut) {
    neuron->iTotalLinksOut += INIT_VECT_SIZE;
    neuron->vLinksOut = realloc(neuron->vLinksOut,
                         neuron->iTotalLinksOut * sizeof(*(neuron->vLinksOut)));
  }
  neuron->vLinksOut[neuron->iNumLinksOut++] = tmpLink;
}

/* add a link to the list of incomming links, reallocate if necessary */
void phenotypeAddLinkIn(sNeuron * neuron, sLink * tmpLink) {
  if (neuron->iNumLinksIn == neuron->iTotalLinksIn) {
    neuron->iTotalLinksIn += INIT_VECT_SIZE;
    neuron->vLinksIn = realloc(neuron->vLinksIn,
                           neuron->iTotalLinksIn * sizeof(*(neuron->vLinksIn)));
  }
  neuron->vLinksIn[neuron->iNumLinksIn++] = tmpLink;
}

/* creation and deletion of the phenotype corresponding to a particular 
 * genotype => into genome.c
 */

/*******************************************************************************
 * update ANN response
 ******************************************************************************/

double * updateANNresponse(sPhenotype * phen, double * inputs, int nbInputs,
                             int nbOutputs, run_type type) {
  int i, lnk;
  // create a vector to put the outputs into
  double * outputs = (double *) calloc(nbOutputs, sizeof(*outputs));
  int ouputIdx = 0;
  // if the mode is snapshot then we require all the neurons to be
  // iterated through as many times as the network is deep. If the
  // mode is set to active the method can return an output after
  // just one iteration
  int flushCount = 0;

  if (type == SNAPSHOT)
    flushCount = phen->iDepth;
  else
    flushCount = 1;

  // iterate through the network FlushCount times
  for (i = 0; i < flushCount; i++) {
    // clear the output vector
    // outputs.clear();
    ouputIdx = 0;

    // this is an index into the current neuron
    int cNeuron = 0;

    // first set the outputs of the 'input' neurons to be equal
    // to the values passed into the function in inputs
    while (phen->vNeurons[cNeuron]->eNeuronType == INPUT) {
      phen->vNeurons[cNeuron]->dOutput = inputs[cNeuron];
      cNeuron++;
    }
    
    // set the output of the bias to 1
    phen->vNeurons[cNeuron++]->dOutput = 1;

    // then we step through the network a neuron at a time
    while (cNeuron < phen->iNumNeurons) {
      // this will hold the sum of all the inputs x weights
      double sum = 0;
      // sum this neuron's inputs by iterating through all the links into
      // the neuron
      for (lnk = 0; lnk < phen->vNeurons[cNeuron]->iNumLinksIn; lnk++) {
        // get this link's weight
        double weight = phen->vNeurons[cNeuron]->vLinksIn[lnk]->dWeight;
        // get the output from the neuron this link is coming from
        double neuronOutput =
                     phen->vNeurons[cNeuron]->vLinksIn[lnk]->pNeuronIn->dOutput;
        // add to sum
        sum += weight * neuronOutput;
      }
      // now put the sum through the activation function and assign the
      // value to this neuron's output
      phen->vNeurons[cNeuron]->dOutput =
                     sigmoid(sum, phen->vNeurons[cNeuron]->dActivationResponse);

      if (phen->vNeurons[cNeuron]->eNeuronType == OUTPUT) {
        //add to our outputs
        outputs[ouputIdx++] = phen->vNeurons[cNeuron]->dOutput;
      }

      // next neuron
      cNeuron++;
    }
  } // next iteration through the network

  // the network needs to be flushed if this type of update is performed
  // otherwise it is possible for dependencies to be built on the order
  // the training data is presented
  if (type == SNAPSHOT) {
    for (i = 0; i < phen->iNumNeurons; i++)
      phen->vNeurons[i]->dOutput = 0;
  }
  
  // return the outputs
  return outputs;
}

/*******************************************************************************
 * convenient for debug - visualization
 ******************************************************************************/

/* dump phenotype */
void dumpPhenotype(sPhenotype * phen, int idx) {
  int i, j;
  printf("phenotype of genome %d (%d nodes, depth : %d)\n",
         idx, phen->iNumNeurons, phen->iDepth);
  for (i = 0; i < phen->iNumNeurons; i++) {
    printf("%2d - id=%d - pos(%d,%d) - split(%f,%f) - ", i,
           phen->vNeurons[i]->iID,
           phen->vNeurons[i]->iPosX, phen->vNeurons[i]->iPosY,
           phen->vNeurons[i]->dSplitX, phen->vNeurons[i]->dSplitY);
    switch(phen->vNeurons[i]->eNeuronType) {
      case BIAS:
        printf("type=bias   ");
        break;
      case INPUT:
        printf("type=input  ");
        break;
      case HIDDEN:
        printf("type=hidden ");
        break;
      case OUTPUT:
        printf("type=output ");
        break;
      case NONE:
        printf("type=none   ");
        break;
      default:
        printf("TYPE ERROR !");
    }
    printf(" - activationResponse=%f - sumActivation=%f - output=%f",
           phen->vNeurons[i]->dActivationResponse,
           phen->vNeurons[i]->dSumActivation,
           phen->vNeurons[i]->dOutput);
    printf("\n");
    puts("\tLINKS IN:");
    for (j = 0; j < phen->vNeurons[i]->iNumLinksIn; j++) {
      printf("\tlink %d -> %d - weight : %f - recurent : %s\n",
             phen->vNeurons[i]->vLinksIn[j]->pNeuronIn->iID,
             phen->vNeurons[i]->vLinksIn[j]->pNeuronOut->iID,
             phen->vNeurons[i]->vLinksIn[j]->dWeight,
             (phen->vNeurons[i]->vLinksIn[j]->bRecurrent == TRUE)? "yes":"no");
    }
    puts("\tLINKS OUT:");
    for (j = 0; j < phen->vNeurons[i]->iNumLinksOut; j++) {
      printf("\tlink %d -> %d - weight : %f - recurent : %s\n",
             phen->vNeurons[i]->vLinksOut[j]->pNeuronIn->iID,
             phen->vNeurons[i]->vLinksOut[j]->pNeuronOut->iID,
             phen->vNeurons[i]->vLinksOut[j]->dWeight,
             (phen->vNeurons[i]->vLinksOut[j]->bRecurrent == TRUE)? "yes":"no");
    }
  }
}

/*******************************************************************************
 * activation function / sigmoid
 ******************************************************************************/

double sigmoid(double netinput, double activationResponse) {
	return ( 1 / ( 1 + exp(- netinput / activationResponse)));
}
