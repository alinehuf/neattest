//
//  utils.h
//  NeatTest_0.1
//
//  Created by dexter on 23/03/13.
//  Copyright (c) 2013 dex. All rights reserved.
//

#ifndef NeatTest_0_1_utils_h
#define NeatTest_0_1_utils_h

// xcode refuse de compiler avec inline ???

//returns a random float in the range 0 <= n < 1
//inline
double randFloat()	{ return rand() / (RAND_MAX + 1.0); }
//returns a random float in the range -1 < n < 1
//inline
double randClamped() { return randFloat() - randFloat(); }
//returns a random integer between x and y
//inline
int randInt(int x,int y) { return rand()%(y-x+1)+x; }


#endif
