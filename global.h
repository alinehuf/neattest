//
//  global.h
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_global_h
#define NeatTest_0_1_global_h

#include "params.h"

#define INIT_VECT_SIZE 10

typedef enum {E_FALSE = 0, E_TRUE} bool;

// NONE type is used into Innovation table when a link is specified,
// neuron type is then not usefull
typedef enum {INPUT, HIDDEN, OUTPUT, BIAS, NONE} neuron_type;

// phenotype function mode
typedef enum {SNAPSHOT, ACTIVE} run_type;


/*******************************************************************************
 * random generators
 ******************************************************************************/

// xcode refuse de compiler avec inline ???
//inline
double randFloat();
//inline
double randClamped();
//inline
int randInt(int x,int y);

/*******************************************************************************
 * neuron gene
 ******************************************************************************/

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
  sInnovation ** vInnovs;
  int iNumInnovs;
  int iTotalInnovs;
  int iNextNeuronId;
} sInnovTable;

/*******************************************************************************
 * link of the final ANN (phenotype)
 ******************************************************************************/
typedef struct neuron sNeuron;

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
  sLink ** vLinksIn;          // all the links coming into this neuron
  int iNumLinksIn;
  int iTotalLinksIn;
  sLink ** vLinksOut;         // and out
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
  sNeuron ** vNeurons;
  int iNumNeurons;
  int iDepth; //the depth of the network
} sPhenotype;

/*******************************************************************************
 * genome
 ******************************************************************************/

typedef struct genome {
  int iId;
  sNeuronGene ** vNeurons;
  int iNumNeurons;
  int iTotalNeurons;
  sLinkGene ** vLinks;
  int iNumLinks;
  int iTotalLinks;
  sPhenotype * pPhenotype;
  double dFitness;
  double dAjustedFitness;
  double dAmountToSpawn;
  int iNumInputs, iNumOuputs;
  int iNumRecur; // to iterate few more times to take account of recurrent input
} sGenome;

/*******************************************************************************
 * species
 ******************************************************************************/

typedef struct {
  sGenome * sLeader;   // keep a local copy of the first member of this species
  sGenome ** vMembers; // pointers to all the genomes within this species
  int iNumMembers;
  int iTotalMembers;
  int iSpeciesId;      // the species needs an identification number
  double dBestFitness;    // best fitness found so far by this species
  int iGensNoImprovement; // generations since fitness has improved
  int iAge;
  double dSpawnsRqd; // how many of this species should be spawned for next pop
} sSpecies;

/*******************************************************************************
 * population
 ******************************************************************************/

// this structure is used in the creation of a network depth lookup table.
typedef struct splitDepth {
  double val;
  struct splitDepth * cdr;
} sSplitDepth;

typedef struct {
  sGenome ** vGenomes;      // current population (fixed size = iNumIndividuals)
  int iNumGenomes;
  int iNextGenomeId;
  // all the species
  sSpecies ** vSpecies;
  int iNumSpecies;
  int iNextSpeciesId;
  int iTotalSpecies;
  sInnovTable * sInnovTable; // to keep track of innovations
  // current generation
  int iGeneration;
  // adjusted fitness scores
  double dTotFitAdj;
  double dAvFitAdj;
  // index into the genomes for the fittest genome
  int iFittestGenome;
  double dBestEverFitness;
  sParams * sParams;
  // this holds the precalculated split depths. They are used
  // for rendering and also for calculating the flush depth of the network
  // when a phenotype is working in 'snapshot' mode.
  double * vDepths;
  int iNumDepth;
  int iTotalDepths;
} sPopulation;

/*******************************************************************************
 * prototypes population.c
 ******************************************************************************/

// population
sPopulation * createPopulation(int nbIn, int nbOut, sParams * params);
void freePopulation(sPopulation * pop);
// epoch
void epoch(sPopulation * pop);
void resetAndKill(sPopulation * pop);
void sortGenomes(sPopulation * pop);
sGenome * tournamentSelection(sPopulation * pop, const int numComparisons);
// convenient
//bool idNotIntoVect(int id, int * vect, int size);
int cmpLinksByInnovIds(const sLinkGene ** a, const sLinkGene ** b);
int cmpNeuronsByIds(const sNeuronGene ** a, const sNeuronGene ** b);
int cmpGenomesByFitness(const sGenome ** a, const sGenome ** b);
// depth lookup table
void addDepth(sPopulation * pop, double splitY);
int calculateNetDepth(sGenome * gen);

/*******************************************************************************
 * prototypes genome.c
 ******************************************************************************/

// neurons, links and genome creation
sNeuronGene * createNeuronGene(int id, neuron_type t, bool r,
                               double x, double y);
sNeuronGene * copyNeuronGene(sNeuronGene * neuron);
sLinkGene * createLinkGene(int inov, int from, int to, double w, bool e,bool r);
sLinkGene * copyLinkGene(sLinkGene * link);
void genomeAddNeuron(sGenome * gen, sNeuronGene * ng);
void genomeAddLink(sGenome * gen, sLinkGene * lg);
sGenome * createInitialGenome(int id, int nbInputs, int nbOutputs);
sGenome * createEmptyGenome(int id, int nbInputs, int nbOutputs);
sGenome * copyGenome(sGenome * model);
void freeGenome(sGenome * gen);
// convenient for debug - visualization
void dumpGenome(FILE * out, sGenome * gen);

// mutations : weight and sigmoidCurvature mutation, add node or add link
void mutateWeigth(sGenome * gen, sParams * params);
void mutateActivationResponse(sGenome * gen, sParams * params);
void addLink(sGenome * gen, sPopulation * pop);
void addNeuron(sGenome * gen, sPopulation * pop);
// crossover
sGenome * crossover(sGenome * mum, sGenome * dad);
// genome compatibility distance
double getCompatibilityScore(sGenome * gen1, sGenome * gen2, sParams * p);

// tools to manipulate the genome
int getNeuronPos(sGenome * gen, const int id);
bool duplicateLink(sGenome * gen, const int neuron_id1, const int neuron_id2);
bool alreadyHaveThisNeuronId(sGenome * gen, const int id);

/*******************************************************************************
 * prototypes phenotype.c
 ******************************************************************************/

sNeuron * createNeuron(int id, neuron_type type, double x, double y, double act);
sLink * createLink(double weight, sNeuron * in, sNeuron * out, bool rec);
void phenotypeAddLinkOut(sNeuron * neuron, sLink * tmpLink);
void phenotypeAddLinkIn(sNeuron * neuron, sLink * tmpLink);
// create / freePhenotype
sPhenotype * createPhenotype(sGenome * gen, int depth);
void freePhenotype(sGenome * gen);
// update ANN
double * updateANNresponse(sPhenotype * phen, double * inputs, int nbInputs,
                           int nbOutputs, run_type type);
double sigmoid(double netinput, double activationResponse);
// convenient for debug
void dumpPhenotype(sPhenotype * phen, int idx);

/*******************************************************************************
 * prototypes species.c
 ******************************************************************************/

sSpecies * createSpecies(sGenome * firstOrg, int specId, int iNumIndividuals);
void freeSpecies(sPopulation * pop);
void addMember(sSpecies * spec, sGenome * newMember);
void purgeSpecies(sSpecies * spec);
// convenient for debug
void dumpSpecies(FILE * out, sPopulation * pop);
// speciation
void addOneSpecies(sPopulation * pop, sGenome * firstOrg);
void removeOneSpecies(sPopulation * pop, int id);
void speciateAndCalculateSpawnLevels(sPopulation * pop);
// tools to manipulate the species
void adjustFitnesses(sSpecies * spec, sParams * param);
void speciesSpawnAmount(sSpecies * spec);
sGenome * randomAmongBest(sSpecies * spec, double dSurvivalRate);

/*******************************************************************************
 * prototypes innovation.c
 ******************************************************************************/

int checkInnovation(sInnovTable * innovTable, int in, int out, innov_type type);
int createNewInnov(sInnovTable * innovTable, int in, int out, innov_type type,
                   neuron_type ntype);
sInnovTable * createNewInnovTable( sNeuronGene ** neurons, int numNeurons,
                                  sLinkGene ** links, int numLinks );
void freeInnovTable(sInnovTable * innovTable);
void dumpInnovTable(FILE * out, sInnovTable * innovTable);

#endif
