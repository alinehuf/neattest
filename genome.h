//
//  genome.h
//  NeatTest_0.1
//
//  Created by dexter on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_genome_h
#define NeatTest_0_1_genome_h

#include "utils.h"
#include "innovation.h"
#include "params.h"

/*******************************************************************************
 * genome
 ******************************************************************************/

typedef struct genome {
  int id;
  sNeuronGene * vNeurons;
  int iNumNeurons;
  int iTotalNeurons;
  sLinkGene * vLinks;
  int iNumLinks;
  int iTotalLinks;
  //NeuralNet * sPhenotype;
  double dFitness;
  double dAjustedFitness;
  //double dAmountToSpawn;
  int iNumInputs, iNumOuputs;
} sGenome;

/*******************************************************************************
 * prototypes
 ******************************************************************************/

// neurons, links and genome creation
sNeuronGene createNeuronGene(int id, neuron_type t, bool r, double x, double y);
sLinkGene createLinkGene(int inov, int from, int to, double w, bool e, bool r);
void genomeAddNeuron(sGenome * gen, sNeuronGene * ng);
void genomeAddLink(sGenome * gen, sLinkGene * lg);
sGenome createInitialGenome(int id, int nbInputs, int nbOutputs);
sGenome createEmptyGenome(int id, int nbInputs, int nbOutputs);
void freeGenome(sGenome * gen);
// convenient for debug - visualization
void dumpGenome(sGenome gen);
// mutations : weight and sigmoidCurvature mutation, add node or add link
void mutateWeigth(sGenome * gen, double weightMutationRate,
                  double probabilityWeightReplaced,
                  double maxWeightPertubation);
void mutateActivationResponse(sGenome * gen, double mut_rate,
                              double maxPertubation);
void addLink(sGenome * gen, double mutationRate, double chanceOfLooped,
             sInnovTable * innovTable, int numTrysToFindLoop,
             int numTrysToAddLink);
void addNeuron(sGenome * gen, double mutationRate, sInnovTable * innovTable,
               int numTrysToFindOldLink);
// genome compatibility distance
double getCompatibilityScore(const sGenome gen1, const sGenome gen2,
                                                               const sParams p);
// tools to manipulate the genome
int getNeuronPos(sGenome gen, const int id);
bool duplicateLink(sGenome gen, const int neuron_id1, const int neuron_id2);
bool alreadyHaveThisNeuronId(sGenome gen, const int id);
#endif
