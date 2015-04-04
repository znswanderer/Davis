/*
C-part for 
Molecular Dynamics Simulation on a Sphere
To be used with python.

Particles on a unit-sphere.
Radius of sphere is R = 1.

MIT license

Copyright(c) 2015 Tim Scheffler

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "davis.h"


Cells *Cells_new(int binning, int num_particles)
{
  Cells *self = (Cells*)malloc(sizeof(Cells));
  // Box is 2.0, not 1.0 because we need to embed sphere with R=1
  self->dr = 2.0 / binning;
  self->binning = binning;
  self->num_cells = binning * binning * binning;
  self->cells = (int*)malloc(sizeof(int) * self->num_cells);
  return self;
}

void Cells_free(Cells *self)
{
  free(self->cells);
  free(self);
}

void Cells_clear(Cells *self)
{
  for (int i=0; i<self->num_cells; i++) {
    self->cells[i] = -1;
  }
}


void dvs_advance(int nparticles, Particle *ps, double dt)
{
  double dt_hlf = 0.5 * dt;
  for (int i=0; i<nparticles; i++) {
    Particle *p = &(ps[i]);
    vec old_r = p->r;
    vec_inc(p->v, vec_sclr(dt_hlf, p->a));
    vec_inc(p->r, vec_sclr(dt, p->v));
    
    // RATTLE_r (for sphere of R=1)
    double r0_dot_r = vec_dot(old_r, p->r);
    double r0_sqr = vec_dot(p->r, p->r);
    // TODO beware of negative sqrt arguments!
    // (But in this case the timestep dt is way too large anyway)
    double lambda = -r0_dot_r + sqrt(1.0 - r0_sqr + r0_dot_r*r0_dot_r);
    
    vec_inc(p->r, vec_sclr(lambda, old_r));
    vec_inc(p->v, vec_sclr(lambda/dt, old_r));
  }
}


void dvs_correct(int nparticles, Particle *ps, double dt)
{
  double dt_hlf = 0.5 * dt;
  for (int i=0; i<nparticles; i++) {
    Particle *p = &(ps[i]);
    vec_inc(p->v, vec_sclr(dt_hlf, p->a));
    // RATTLE_v
    double lambda = -vec_dot(p->v, p->r);
    vec_inc(p->v, vec_sclr(lambda, p->r));
  }
}

inline void dvs_calc_force(Particle *p, Particle *q, 
			   double cutoff, double gamma, 
			   Stats *stats)
{
  vec dr = vec_sub(p->r, q->r);
  double r2 = vec_magnitude2(dr);
  stats->ww_counter++;
  double cutoff2 = cutoff*cutoff;
  if (r2 < cutoff2) {
    double r = sqrt(r2);
    stats->real_ww_counter++;
    
    // Coulomb force
    double force_mag = 1.0/r2 - 1.0/cutoff2;
    stats->E_pot += 1.0/r + r/cutoff2 - 2.0/cutoff;
    vec force = vec_sclr(force_mag / r, dr);
    // Damping
    // Note: the damping is done with v(t + dt/2).  To be correct, we
    // might should use v(t+dt) (which is not available at this
    // point). But as we just model this dissipative force in a
    // handwaving sort anyway this should not matter. The damping has
    // not real physical meaning and is just introduced to cool the
    // system down.
    vec dv = vec_sub(p->v, q->v);
    vec damping = vec_sclr(-gamma , dv);
    vec_inc(force, damping);

    // TODO Force is still completely 3d
    vec_inc(p->a, force);
    vec_dec(q->a, force);
  }
}

void dvs_populate_cells(int nparticles, Particle *ps, Cells *cells)
{
  int bin_x, bin_y, bin_z;
  int L = cells->binning;
  int L2 = L*L;

  Cells_clear(cells);
  int i = 0;
  for (Particle *p=ps; p<ps+nparticles; p++) {
    p->a = VEC_ZERO;
    // The simulation box is [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]]
    bin_x = MAX(0, MIN((int)floor((p->r.x+1.0)/cells->dr), L-1));
    bin_y = MAX(0, MIN((int)floor((p->r.y+1.0)/cells->dr), L-1));
    bin_z = MAX(0, MIN((int)floor((p->r.z+1.0)/cells->dr), L-1));
    int cellNum = bin_x + bin_y*L + bin_z*L2;
    p->next = cells->cells[cellNum];
    cells->cells[cellNum] = i;
    i++;
  }
}

void dvs_calc_forces(Particle *ps, Cells *cells, int cell0, int cell1,
		     double cutoff, double gamma,
		     Stats *stats)
{
  int L = cells->binning;
  int L2 = L*L;

  // Cells without periodic boundary conditions!
  for (int z=0; z<L; z++) {
    for (int y=0; y<L; y++) {
      for (int x=0; x<L; x++) {
	int this_cell = x + L*y + L2*z;
	if ((this_cell < cell0) || (this_cell >= cell1)) continue;
	for (int i=cells->cells[this_cell]; i!=-1; i=ps[i].next) {
	  // other particles in main cell (q > p)
	  for (int j=ps[i].next; j!=-1; j=ps[j].next) {
	    dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
	  }
	  // next neaerest cells
	  int nnx, nny, nnz;
	  nny = y;
	  nnz = z;
	  // to the right
	  nnx = x + 1;
	  if (nnx < L) {
	    for (int j=cells->cells[nnx + nny*L + nnz*L2]; j!=-1; j=ps[j].next) {
	      dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
	    }
	  }
	  // upper row
	  nny = y + 1;
	  if (nny < L) {
	    // ... upper-left
	    nnx = x - 1;
	    if (nnx >= 0) {
	      for (int j=cells->cells[nnx + nny*L + nnz*L2]; j!=-1; j=ps[j].next) {
		dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
	      }
	    }
	    // ... upper-middle
	    nnx = x;
	    for (int j=cells->cells[nnx + nny*L + nnz*L2]; j!=-1; j=ps[j].next) {
	      dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
	    }
	    // ... upper-right
	    nnx = x + 1;
	    if (nnx < L) {
	      for (int j=cells->cells[nnx + nny*L + nnz*L2]; j!=-1; j=ps[j].next) {
		dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
	      }
	    }
	  }
	  // plus possibly 9 NN boxes in the z+1 direction
	  nnz = z + 1;
	  if (nnz < L ) {
	    for (nny=y-1; nny<=y+1; nny++) {
	      if ((nny >= 0) && (nny < L)) {
		for (nnx=x-1; nnx<=x+1; nnx++) {
		  if ((nnx >= 0) && (nnx < L)) {
		    for (int j=cells->cells[nnx + nny*L + nnz*L2]; j!=-1; j=ps[j].next) {
		      dvs_calc_force(ps + i, ps + j, cutoff, gamma, stats);
		    }
		  } // nasty ...
		}
	      } // ... curly ...
	    }
	  }
	} //... braces.
      }
    }
  }
}

void dvs_visualise_positions(int nparticles, Particle *ps, double *target)
{
  for (int i=0; i < nparticles; i++) {
      target[i*3    ] = ps[i].r.x;
      target[i*3 + 1] = ps[i].r.y;
      target[i*3 + 2] = ps[i].r.z;
  }
}


void dvs_copy_particles(int nparticles, Particle *src, Particle *dst)
{
  memcpy(dst, src, nparticles * sizeof(Particle));
}


void dvs_collect_forces(int nparticles, Particle *accu, Particle *part)
{
  for (int i=0; i<nparticles; i++) {
    vec_inc(accu[i].a, part[i].a);
  }
}

