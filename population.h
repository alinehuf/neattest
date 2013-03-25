//
//  population.h
//  NeatTest_0.1
//
//  Created by dexter on 25/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_population_h
#define NeatTest_0_1_population_h

sGenome crossover(sGenome mum, sGenome dad, int babyId);
bool idNotIntoVect(int id, int * vect, int size);
int cmpNeuronsByIds(const sNeuronGene * a, const sNeuronGene * b);

#endif
