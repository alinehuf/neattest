//
//  main.c
//  NeatTest_0.1
//
//  Created by AlineHUF on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//
//  largely based on Mat Buckland 2002

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "genome.h"
#include "population.h"

int
main(int argc, char * argv[])
{
  srand((unsigned int) time(NULL)); // inits random

  sParams params = loadConf(NULL);
  sGenome gen1 = createInitialGenome(0, 8, 3);
  sGenome gen2 = createInitialGenome(1, 8, 3);
  sInnovTable innovTable = createNewInnovTable(gen1.vNeurons, gen1.iNumNeurons,
                                               gen1.vLinks, gen1.iNumLinks);
  dumpGenome(gen1);
  dumpGenome(gen2);
  int i;
  for (i = 0; i < 100; i ++) {
    addNeuron(&gen1, params.dChanceAddNode, &innovTable,
              params.iNumTrysToFindOldLink);
    addNeuron(&gen2, params.dChanceAddNode, &innovTable,
              params.iNumTrysToFindOldLink);
  }
  for (i = 0; i < 100; i ++) {
    addLink(&gen1, params.dChanceAddLink, params.dChanceAddRecurrentLink,
            &innovTable, params.iNumTrysToFindLoop, params.iNumTrysToAddLink);
    addLink(&gen2, params.dChanceAddLink, params.dChanceAddRecurrentLink,
            &innovTable, params.iNumTrysToFindLoop, params.iNumTrysToAddLink);
  }
  puts("------------MUTATION--------");
  dumpGenome(gen1);
  dumpGenome(gen2);
  printf("compatibility : %f \n",
         getCompatibilityScore(gen1, gen2, &params));
  sGenome gen3 = crossover(gen1, gen2);
  gen3.iId = 2;
  puts("------------CROSSOVER--------");
  dumpGenome(gen3);
  
  puts("------------INNOVATIONS--------");
  dumpInnovTable(innovTable);
  
  // free memory
  freeGenome(&gen1);
  freeGenome(&gen2);
  return 0;
}