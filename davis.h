#ifndef DAVIS_H
#define DAVIS_H

#include "vec3d.h"
#include "helpers.h"

/*
  Multithread programming is hard.
  In the first version I had direct links in Particle->next and
  Cells->cells to Particle pointers. But for the different threads I
  simply copied the whole particle array for each thread not realising, 
  that the pointers do point to the original locations.
  So I now use relative positions stored as int. This is a bit slower.
  
 */

typedef struct {
  vec r;
  vec v;
  vec a;
  long next;
} Particle;

typedef struct {
  double dr;
  long binning;
  long num_cells;
  long *cells;
} Cells;

typedef struct {
  long ww_counter;
  long real_ww_counter;
  double E_pot;
} Stats;

#endif
