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
#include "params.h"
#include "genome.h"

int
main(int argc, char * argv[])
{
  srand((unsigned int) time(NULL)); // inits random

  sParams params = loadConf(NULL);
  sGenome gen = createInitialGenome(0, 8, 3);
  sInnovTable innovTable = createNewInnovTable(gen.vNeurons, gen.iNumNeurons,
                                               gen.vLinks, gen.iNumLinks);
  dumpGenome(gen);
  int i;
  for (i = 0; i < 100; i ++)
    addNeuron(&gen, params.dChanceAddNode, &innovTable,
              params.iNumTrysToFindOldLink);
  for (i = 0; i < 100; i ++)
    addLink(&gen, params.dChanceAddLink, params.dChanceAddRecurrentLink,
            &innovTable, params.iNumTrysToFindLoop, params.iNumTrysToAddLink);
  puts("------------MUTATION--------");
  dumpGenome(gen);
  // free memory
  freeGenome(&gen);
  puts("------------INNOVATIONS--------");
  dumpInnovTable(innovTable);
  return 0;
}