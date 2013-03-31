//
//  phenotype.c
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "phenotype.h"

/*******************************************************************************
 * neuron of the final ANN (phenotype)
 ******************************************************************************/

sNeuron createNeuron(int id, neuron_type type, double x, double y, double act) {
  sNeuron neuron;
  neuron.iID = id;
  neuron.dActivationResponse = act;
  neuron.iNumLinksIn = 0;
  neuron.iTotalLinksIn = INIT_VECT_SIZE;
  neuron.vLinksIn = (sLink *) calloc(INIT_VECT_SIZE, sizeof(*neuron.vLinksIn));
  neuron.iNumLinksOut = 0;
  neuron.iTotalLinksOut = INIT_VECT_SIZE;
  neuron.vLinksOut = (sLink *) calloc(INIT_VECT_SIZE, sizeof(*neuron.vLinksIn));
  neuron.dSumActivation = 0;
  neuron.dOutput = 0;
  neuron.iPosX = 0;
  neuron.iPosY = 0;
  neuron.dSplitX = x;
  neuron.dSplitY = y;
  return neuron;
}

/*******************************************************************************
 * link of the final ANN (phenotype)
 ******************************************************************************/

sLink createLink(double weight, sNeuron * in, sNeuron * out, bool rec) {
  sLink link;
  link.pNeuronIn = in;
  link.pNeuronOut = out;
  link.dWeight = weight;
  link.bRecurrent = rec;
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
  neuron->vLinksOut[neuron->iNumLinksOut++] = *tmpLink;
}

/* add a link to the list of incomming links, reallocate if necessary */
void phenotypeAddLinkIn(sNeuron * neuron, sLink * tmpLink) {
  if (neuron->iNumLinksIn == neuron->iTotalLinksIn) {
    neuron->iTotalLinksIn += INIT_VECT_SIZE;
    neuron->vLinksIn = realloc(neuron->vLinksIn,
                           neuron->iTotalLinksIn * sizeof(*(neuron->vLinksIn)));
  }
  neuron->vLinksIn[neuron->iNumLinksIn++] = *tmpLink;
}

/* creation and deletion of the phenotype corresponding to a particular 
 * genotype => into genome.c
 */




