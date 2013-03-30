//
//  genome.c
//  NeatTest_0.1
//
//  Created by dexter on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "genome.h"

/*******************************************************************************
 * neurons, links and genome creation
 ******************************************************************************/ 

/* Creates a neuron gene.
 */
sNeuronGene createNeuronGene(int id, neuron_type type, bool r,
                             double x, double y) {
  sNeuronGene ng;
  ng.iId = id;
  ng.eNeuronType = type;
  ng.bRecurrent = r;
  ng.dSigmoidCurvature = 1;
  ng.dSplitX = x;
  ng.dSplitY = y;
  return ng;
}

/* Creates a link gene.
 */
sLinkGene createLinkGene(int innov, int from, int to, double w, bool e, bool r){
  sLinkGene lg;
  lg.iInnovId = innov;
  lg.iFromNeuron = from;
  lg.iToNeuron = to;
  lg.dWeight = w;
  lg.bEnabled = e;
  lg.bRecurrent = r;

  return lg;
}

void genomeAddNeuron(sGenome * gen, sNeuronGene * ng) {
  if (gen->iNumNeurons == gen->iTotalNeurons) {
    gen->iTotalNeurons += INIT_VECT_SIZE;
    gen->vNeurons = realloc(gen->vNeurons,
                            gen->iTotalNeurons * sizeof(*(gen->vNeurons)));
  }
  gen->vNeurons[gen->iNumNeurons++] = *ng;
}

void genomeAddLink(sGenome * gen, sLinkGene * lg) {
  if (gen->iNumLinks == gen->iTotalLinks) {
    gen->iTotalLinks += INIT_VECT_SIZE;
    gen->vLinks = realloc(gen->vLinks,
                          gen->iTotalLinks * sizeof(*(gen->vLinks)));
  }
  gen->vLinks[gen->iNumLinks++] = *lg;
}

/* Creates a minimal genome where there are output + input neurons + bias
 * and each input neuron is connected to each output neuron.
 */
sGenome createInitialGenome(int id, int nbInputs, int nbOutputs) {
  sGenome gen;
  gen.id = id;
  // initialize fitness
  gen.dFitness = 0;
  gen.dAjustedFitness = 0;
  // remember number of inputs and outputs
  gen.iNumInputs = nbInputs;
  gen.iNumOuputs = nbOutputs;

  // init total and current number of neurons and links
  gen.iTotalNeurons = gen.iTotalLinks = INIT_VECT_SIZE;
  gen.iNumNeurons = gen.iNumLinks = 0;

  // allocate memory for vectors of neurons and links
  gen.vNeurons = calloc(gen.iTotalNeurons, sizeof(*(gen.vNeurons)));
  gen.vLinks = calloc(gen.iTotalLinks, sizeof(*(gen.vLinks)));

  int i, j;
  sNeuronGene ng;
  sLinkGene lg;
  //create the input neurons
  double inSlice = 1 / (double) nbInputs;
	for (i = 0; i < nbInputs; i++) {
    ng = createNeuronGene(i, INPUT, FALSE, i * inSlice, 0);
    genomeAddNeuron(&gen, &ng);
  }
  //create the bias
  ng = createNeuronGene(nbInputs, BIAS, FALSE, 1, 0);
  genomeAddNeuron(&gen, &ng);

  //create the output neurons
  double outSlice = 1 / (double) (nbOutputs - 1);
  for (i = 0; i < nbOutputs; i++) {
    ng = createNeuronGene(i + nbInputs + 1, OUTPUT, FALSE, i * outSlice, 1);
		genomeAddNeuron(&gen, &ng);
  }
	//create the link genes, connect each input neuron to each output neuron and
	//assign a random weight -1 < w < 1
	for (i = 0; i < nbInputs + 1; i++) {
		for (j = 0; j < nbOutputs; j++) {
      lg = createLinkGene(nbInputs + nbOutputs + 1 + gen.iNumLinks,
                          gen.vNeurons[i].iId,
                          gen.vNeurons[nbInputs + j + 1].iId,
                          randClamped(), TRUE, FALSE);
			genomeAddLink(&gen, &lg);
    }
  }
  return gen;
}

sGenome createEmptyGenome(int id, int nbInputs, int nbOutputs) {
  sGenome gen;
  gen.id = id;
  // initialize fitness
  gen.dFitness = 0;
  gen.dAjustedFitness = 0;
  // remember number of inputs and outputs
  gen.iNumInputs = nbInputs;
  gen.iNumOuputs = nbOutputs;
  // init total and current number of neurons and links
  gen.iTotalNeurons = gen.iTotalLinks = INIT_VECT_SIZE;
  gen.iNumNeurons = gen.iNumLinks = 0;
  // allocate memory for vectors of neurons and links
  gen.vNeurons = calloc(gen.iTotalNeurons, sizeof(*(gen.vNeurons)));
  gen.vLinks = calloc(gen.iTotalLinks, sizeof(*(gen.vLinks)));
  return gen;
}


/* Frees memory allocated for the genome
 */
void freeGenome(sGenome * gen) {
  free(gen->vNeurons);
  free(gen->vLinks);
}

/*******************************************************************************
 * convenient for debug - visualization
 ******************************************************************************/

/* dump genome */
void dumpGenome(sGenome gen) {
  printf("genome %d (%d inputs, %d outputs, %d nodes, %d links) :\n",
        gen.id, gen.iNumInputs, gen.iNumOuputs, gen.iNumNeurons, gen.iNumLinks);
  int i;
  for (i = 0; i < gen.iNumNeurons; i++) {
    printf("Neurons %-2d - id=%-2d : ", i,
           gen.vNeurons[i].iId);
    switch(gen.vNeurons[i].eNeuronType) {
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
        printf("TYPE ERROR !!!");
    }
    printf("recur=%d activation=%f splitX=%f splitY=%f\n",
           gen.vNeurons[i].bRecurrent, gen.vNeurons[i].dSigmoidCurvature,
           gen.vNeurons[i].dSplitX, gen.vNeurons[i].dSplitY);
  }
  for (i = 0; i < gen.iNumLinks; i++) {
    printf("link %-2d - inov=%d : from %-2d to %-2d enabled=%i recur=%i weight=%+f\n",
           i, gen.vLinks[i].iInnovId, gen.vLinks[i].iFromNeuron,
           gen.vLinks[i].iToNeuron, gen.vLinks[i].bEnabled,
           gen.vLinks[i].bRecurrent, gen.vLinks[i].dWeight);
  }
}

/*******************************************************************************
 * mutations : weight and sigmoidCurvature mutation, add node or add link
 ******************************************************************************/

void mutateWeigth(sGenome * gen, double weightMutationRate,
                  double probabilityWeightReplaced,
                  double maxWeightPertubation) {
  int i;
  for (i = 0; i < gen->iNumLinks; i++) {
    //do we mutate this gene?
		if (randFloat() < weightMutationRate) {
      //do we change the weight to a completely new weight?
			if (randFloat() < probabilityWeightReplaced) {
        //change the weight using the random distribtion defined by 'type'
				gen->vLinks[i].dWeight = randClamped();
      } else {
				//perturb the weight
				gen->vLinks[i].dWeight += randClamped() * maxWeightPertubation;
      }
    }
  }
}

void mutateActivationResponse(sGenome * gen, double activationMutationRate,
                              double maxActivPertubation) {
  int i;
  for (i = 0; i < gen->iNumNeurons; i++) {
    //do we mutate this neuron gene ?
		if (randFloat() < activationMutationRate)
      //perturb the sigmoid curvature
      gen->vNeurons[i].dSigmoidCurvature += randClamped() * maxActivPertubation;
  }
}

void addLink(sGenome * gen, double chanceAddLink, double chanceAddRecurrentLink,
             sInnovTable * innovTable, int numTrysToFindLoop,
             int numTrysToAddLink) {
  //do we really add a link ?
  if (randFloat() > chanceAddLink) return;
  
  //define holders for the two neurons to be linked. If we have find two
  //valid neurons to link these values will become >= 0.
  int neuron_id1 = -1;
  int neuron_id2 = -1;
  //flag set if a recurrent link is selected (looped or normal)
  bool bRecurrent = FALSE;

  //first test to see if an attempt should be made to create a
  //link that loops back into the same neuron
  if (randFloat() < chanceAddRecurrentLink) {
    //YES: try numTrysToFindLoop times to find a neuron that is not an
    //input or bias neuron and that does not already have a loopback
    //connection
    while(numTrysToFindLoop--) {
      //grab a random neuron that is not an input or bias neuron
      int neuronPos = randInt(gen->iNumInputs + 1, gen->iNumNeurons - 1);
      //check to make sure the neuron does not already have a loopback link
      if (!gen->vNeurons[neuronPos].bRecurrent) {
        // DANS LE CODE DE BUKLAND... INUTILE ??? -------------------------------------???
        if ( gen->vNeurons[neuronPos].eNeuronType == INPUT ||
             gen->vNeurons[neuronPos].eNeuronType == BIAS ) {
          fprintf(stderr, "addLink() : error, try to make recurrent " \
                 "an input or bias neuron!");
        }
        neuron_id1 = neuron_id2 = gen->vNeurons[neuronPos].iId;
        gen->vNeurons[neuronPos].bRecurrent = TRUE;
        bRecurrent = TRUE;
        numTrysToFindLoop = 0;
      }
    }
  } else {
    //No: try to find two unlinked neurons. Make numTrysToAddLink attempts
    while(numTrysToAddLink--) {
      //choose two neurons, the second must not be an input or a bias
      neuron_id1 = gen->vNeurons[randInt(0, gen->iNumNeurons - 1)].iId;
      neuron_id2 =
        gen->vNeurons[randInt(gen->iNumInputs + 1, gen->iNumNeurons -1)].iId;

      // DANS LE CODE DE BUKLAND... INUTILE ??? -------------------------------------???
      if (neuron_id2 == 2)
        fprintf(stderr, "addLink() : error, neuron_id2 == 2");

      //make sure these two are not already linked and that they are
      //not the same neuron
      if ( !( duplicateLink(*gen, neuron_id1, neuron_id2) ||
              neuron_id1 == neuron_id2 ) ) {
        numTrysToAddLink = 0;
      }  else {
        neuron_id1 = neuron_id2 = -1;
      }
    }
  }
  //return if unsuccessful in finding a link
  if ( (neuron_id1 < 0) || (neuron_id2 < 0) ) return;

  //check to see if we have already created this innovation
  int id = checkInnovation(innovTable, neuron_id1, neuron_id2, NEW_LINK);
  //is this link recurrent?
  if ( gen->vNeurons[getNeuronPos(*gen, neuron_id1)].dSplitY >
       gen->vNeurons[getNeuronPos(*gen, neuron_id2)].dSplitY ) {
    bRecurrent = TRUE;
  }
  if (id < 0) { //we need to create a new innovation
    createNewInnov(innovTable, neuron_id1, neuron_id2, NEW_LINK, NONE);
    //then create the new gene
    id = innovTable->iNumInnovs - 1;
    sLinkGene lg = createLinkGene(id, neuron_id1, neuron_id2, randClamped(),
                                  TRUE, bRecurrent);
    genomeAddLink(gen, &lg);
  } else {
    //the innovation has already been created so all we need to
    //do is create the new gene using the existing innovation ID
    sLinkGene lg = createLinkGene(id, neuron_id1, neuron_id2, randClamped(),
                                  TRUE, bRecurrent);
    genomeAddLink(gen, &lg);
  }
}

// this function adds a neuron to the genotype by examining the network,
// splitting one of the links and inserting the new neuron.
void addNeuron(sGenome * gen, double dChanceAddNode, sInnovTable * innovTable,
               int numTrysToFindOldLink) {
  //just return dependent on mutation rate
  if (randFloat() > dChanceAddNode) return;

  //if a valid link is found into which to insert the new neuron
  //this value is set to true.
  bool bDone = FALSE;

  //this will hold the index into gen->vLinks of the chosen link gene
  int chosenLink = 0;

  //first a link is chosen to split. If the genome is small the code makes
  //sure one of the older links is split to ensure a chaining effect does
  //not occur. Here, if the genome contains less than 5 hidden neurons it
  //is considered to be too small to select a link at random
  if (gen->iNumNeurons < gen->iNumInputs + gen->iNumOuputs + 5) {
    while(numTrysToFindOldLink--) {
      //choose a link with a bias towards the older links in the genome
      chosenLink = randInt(0, gen->iNumLinks - 1 - (int) sqrt(gen->iNumLinks));
      //make sure the link is enabled and that it is not a recurrent link
      //or has a bias input
      int fromNeuron = gen->vLinks[chosenLink].iFromNeuron;
      if ( gen->vLinks[chosenLink].bEnabled &&
           !gen->vLinks[chosenLink].bRecurrent &&
           gen->vNeurons[getNeuronPos(*gen, fromNeuron)].eNeuronType != BIAS) {
        bDone = TRUE;
        numTrysToFindOldLink = 0;
      }
    }
    if (!bDone) //failed to find a decent link
      return;
  } else {
    //the genome is of sufficient size for any link to be acceptable
    while (!bDone) {
      chosenLink = randInt(0, gen->iNumLinks - 1);
      //make sure the link is enabled and that it is not a recurrent link
      //or has a BIAS input
      int fromNeuron = gen->vLinks[chosenLink].iFromNeuron;

      if ( gen->vLinks[chosenLink].bEnabled &&
           !gen->vLinks[chosenLink].bRecurrent &&
           gen->vNeurons[getNeuronPos(*gen, fromNeuron)].eNeuronType != BIAS ) {
        bDone = TRUE;
      }
    }
  }
  // disable this gene
  gen->vLinks[chosenLink].bEnabled = FALSE;
  // grab the weight from the gene (we want to use this for the weight of
  // one of the new links so that the split does not disturb anything the
  // NN may have already learned...
  double originalWeight = gen->vLinks[chosenLink].dWeight;
  // identify the neurons this link connects
  int from = gen->vLinks[chosenLink].iFromNeuron;
  int to = gen->vLinks[chosenLink].iToNeuron;
  // calculate the depth and width of the new neuron. We can use the depth
  // to see if the link feeds backwards or forwards
  double new_split_x = ( gen->vNeurons[getNeuronPos(*gen, from)].dSplitX +
                      gen->vNeurons[getNeuronPos(*gen, to)].dSplitX ) / 2;
  double new_split_y = ( gen->vNeurons[getNeuronPos(*gen, from)].dSplitY +
                      gen->vNeurons[getNeuronPos(*gen, to)].dSplitY ) / 2;
  // Now to see if this innovation has been created previously by
  // another member of the population
  int id = checkInnovation(innovTable, from, to, NEW_NEURON);
  /* it is possible for NEAT to repeatedly do the following:
   1. Find a link. Lets say we choose link 1 to 5
   2. Disable the link,
   3. Add a new neuron and two new links
   4. The link disabled in Step 2 maybe re-enabled when this genome
   is recombined with a genome that has that link enabled.
   5  etc etc
   Therefore, this function must check to see if a neuron ID is already
   being used. If it is then the function creates a new innovation
   for the neuron. */
  if (id >= 0) {
    int neuron_id = innovTable->vInnovs[id].iNeuronId;
    if (alreadyHaveThisNeuronId(*gen, neuron_id))
      id = -1;
  }
  if (id < 0) {
    // add the innovation for the new neuron
    int neuron_id = createNewInnov(innovTable, from, to, NEW_NEURON, HIDDEN);
    // create the new neuron gene and add it.
    sNeuronGene ng = createNeuronGene( neuron_id, HIDDEN, FALSE, new_split_x,
                                       new_split_y );
		genomeAddNeuron(gen, &ng);
    // Two new link innovations are required, one for each of the
    // new links created when this gene is split.
    //-----------------------------------first link
    //get the next innovation ID
    int idLink1 = innovTable->iNumInnovs;
    //create the new innovation
    createNewInnov(innovTable, from, neuron_id, NEW_LINK, NONE);
    //create the new link gene
    sLinkGene link1 = createLinkGene(idLink1, from, neuron_id, 1.0, TRUE,FALSE);
    genomeAddLink(gen, &link1);
    //-----------------------------------second link
    //get the next innovation ID
    int idLink2 = innovTable->iNumInnovs;
    //create the new innovation
    createNewInnov(innovTable, neuron_id, to, NEW_LINK, NONE);
    //create the new link gene
    sLinkGene link2 = createLinkGene(idLink2, neuron_id, to, originalWeight,
                                     TRUE, FALSE);
    genomeAddLink(gen, &link2);
  } else {
    // this innovation has already been created so grab the relevant neuron
    // and link from the innovation table
    int neuron_id = innovTable->vInnovs[id].iNeuronId;
    //get the innovation IDs for the two new link genes.
    int idLink1 = checkInnovation(innovTable, from, neuron_id, NEW_LINK);
    int idLink2 = checkInnovation(innovTable, neuron_id, to, NEW_LINK);
    // this should never happen because the innovations *should* have already
    // occurred
    if (idLink1 < 0 || idLink2 < 0) {
      fprintf(stderr, "addNeuron() : error links around a neurons unknown in " \
                      "innovTable\n");
      exit(EXIT_FAILURE);
    }
    //now we need to create 2 new genes to represent the new links
    sLinkGene link1 = createLinkGene(idLink1, from, neuron_id, 1.0, TRUE,FALSE);
    sLinkGene link2 = createLinkGene(idLink2, neuron_id, to, originalWeight,
                                     TRUE, FALSE);
    genomeAddLink(gen, &link1);
    genomeAddLink(gen, &link2);
    //create the new neuron
    sNeuronGene ng = createNeuronGene( neuron_id, HIDDEN, FALSE, new_split_x,
                                       new_split_y );
    //and add it
    genomeAddNeuron(gen, &ng);
  }
}

/*******************************************************************************
 * tools to manipulate the genome
 ******************************************************************************/


// returns true if the link is already part of the genome
bool duplicateLink(sGenome gen, const int neuron_id1, const int neuron_id2) {
  int i;
  for (i = 0; i < gen.iNumLinks; i ++) {
    if ( gen.vLinks[i].iFromNeuron == neuron_id1 &&
        gen.vLinks[i].iToNeuron == neuron_id2 ) {
      return TRUE;
    }
  }
  return FALSE;
}

//	given a neuron ID just finds its position in gen.vNeurons
int getNeuronPos(sGenome gen, const int id) {
  int i;
  for (i = 0; i < gen.iNumNeurons; i++)
    if (gen.vNeurons[i].iId == id)
      return i;
  fprintf(stderr, "getNeuronPos : neuron id %d doesn't exist !", id);
  return -1;
}

// tests to see if the parameter is equal to any existing neuron ID's.
// Returns true if this is the case.
bool alreadyHaveThisNeuronId(sGenome gen, const int id) {
  int i;
  for (i = 0; i < gen.iNumNeurons; i++)
    if (gen.vNeurons[i].iId == id)
      return TRUE;
  return FALSE;
}

