//
//  global.h
//  NeatTest_0.1
//
//  Created by dexter on 31/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_global_h
#define NeatTest_0_1_global_h

#define INIT_VECT_SIZE 10

typedef enum {FALSE=0, TRUE} bool;

// NONE type is used into Innovation table when a link is specified, neuron type
// is then not usefull
typedef enum {INPUT, HIDDEN, OUTPUT, BIAS, NONE} neuron_type;

#endif
