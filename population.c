//
//  population.c
//  NeatTest_0.1
//
//  Created by dexter on 25/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "genome.h"
#include "population.h"

sGenome crossover(sGenome mum, sGenome dad, int babyId) {
  typedef enum {MUM, DAD} parent_type;
  // first, calculate the genome we will using the disjoint/excess genes from.
  // This is the fittest genome.
  parent_type best;

  // if they are of equal fitness use the shorter (because we want to keep
  // the networks as small as possible)
  if (mum.dFitness == dad.dFitness) {
    // if they are of equal fitness and length just choose one at random
    if (mum.iNumLinks == dad.iNumLinks)
      best = (parent_type) randInt(0, 1);
    else if (mum.iNumLinks < dad.iNumLinks)
      best = MUM;
    else
      best = DAD;
  } else if (mum.dFitness > dad.dFitness) {
    best = MUM;
  } else {
    best = DAD;
  }

  // create en empty genome for the baby
  sGenome baby = createEmptyGenome(babyId, mum.iNumInputs, mum.iNumOuputs);
  // this will hold a copy of the gene we wish to add at each step
  sLinkGene selectedGene;
  parent_type selectedFrom;

  int curMum = 0;
  int curDad = 0;
  // step through each parents genes until we reach the end of both
  while (curMum < mum.iNumLinks || curDad < dad.iNumLinks) {
    
    // the end of mum's genes have been reached
    if (curMum >= mum.iNumLinks && curDad < dad.iNumLinks) {
      // if dad is fittest
      if (best == DAD) { //add dads genes
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto dad's next gene
      curDad++;
    }
    // the end of dads's genes have been reached
    else if ( curDad >= dad.iNumLinks && curMum < mum.iNumLinks) {
      // if mum is fittest
      if (best == MUM) { //add mums genes
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
       }
      // move onto mum's next gene
      curMum++;
    }
    // if mums innovation number is less than dads
    else if (mum.vLinks[curMum].iInnovId < dad.vLinks[curDad].iInnovId) {
      // if mum is fittest add gene
      if (best == MUM) {
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
      }
      // move onto mum's next gene
      curMum++;
    }
    // if dads innovation number is less than mums
    else if (dad.vLinks[curDad].iInnovId < mum.vLinks[curMum].iInnovId) {
      // if dad is fittest add gene
      if (best == DAD) {
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto dad's next gene
      curDad++;
    }
    // if innovation numbers are the same
    else if (mum.vLinks[curMum].iInnovId == dad.vLinks[curDad].iInnovId)
    {
      // grab a gene from either parent
      if (randFloat() < 0.5f) {
        selectedGene = mum.vLinks[curMum];
        selectedFrom = MUM;
      } else {
        selectedGene = dad.vLinks[curDad];
        selectedFrom = DAD;
      }
      // move onto next gene of each parent
      curMum++;
      curDad++;
    }

    // add the selected gene if not already added
    if ( baby.iNumLinks == 0 ||
         baby.vLinks[baby.iNumLinks-1].iInnovId != selectedGene.iInnovId )
      genomeAddLink(&baby, &selectedGene);

    // Check if we already have the nodes referred to in SelectedGene.
    // If not, they need to be added.
    if (!alreadyHaveThisNeuronId(baby, selectedGene.iFromNeuron)) {
      if (selectedFrom == MUM)
        genomeAddNeuron(&baby,
                   &mum.vNeurons[getNeuronPos(mum, selectedGene.iFromNeuron)] );
      else
        genomeAddNeuron(&baby,
                   &dad.vNeurons[getNeuronPos(dad, selectedGene.iFromNeuron)] );
    }
    if (!alreadyHaveThisNeuronId(baby, selectedGene.iToNeuron)) {
      if (selectedFrom == MUM)
        genomeAddNeuron(&baby,
                     &mum.vNeurons[getNeuronPos(mum, selectedGene.iToNeuron)] );
      else
        genomeAddNeuron(&baby,
                     &dad.vNeurons[getNeuronPos(dad, selectedGene.iToNeuron)] );
    }
  } // end while

  //now sort the neurons into order
  qsort(baby.vNeurons, baby.iNumNeurons, sizeof(sNeuronGene),
        (int (*) (const void *, const void *)) cmpNeuronsByIds);

  return baby;
}

/*******************************************************************************
 * convenient
 ******************************************************************************/

bool idNotIntoVect(int id, int * vect, int size) {
  int i;
  for (i = 0; i < size; i++) {
    if (id == vect[i]) return FALSE;
  }
  return TRUE;
}

int cmpNeuronsByIds(const sNeuronGene * a, const sNeuronGene * b) {
  if (a->iId == b->iId) return 0;
  else if (a->iId < b->iId) return -1;
  else return 1;
}
