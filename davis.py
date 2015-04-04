"""
Python part for 
Molecular Dynamics Simulation on a Sphere

Please see README.md for details.

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

"""


import ctypes as c
import numpy as np
import math
import sys
import Queue
import threading
import thread
import itertools as it
import random
import time
from scipy.spatial import ConvexHull
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

clib = c.CDLL("./davis.so")

from Visualiser.gui import Window
import Visualiser.simulation

COLORING_INTERVAL = 4.0      # in seconds
PARALLEL = 2                 # number of parallel workers
#NUM_PARTICLES = 12           # Icosahedron
NUM_PARTICLES = 40000          
POINT_SIZE = 3

class Vec(c.Structure):
    _fields_ = [
        ("x", c.c_double),
        ("y", c.c_double),
        ("z", c.c_double),
    ]
    def __str__(self):
        return "(%f, %f, %f)" % (self.x, self.y, self.z)


class Particle(c.Structure):
    _fields_ = [
        ("r", Vec),
        ("v", Vec),
        ("a", Vec),
        ("next", c.c_void_p),
    ]

class Stats(c.Structure):
    _fields_ = [
        ("ww_counter", c.c_long),
        ("real_ww_counter", c.c_long),
        ("E_pot", c.c_double)
    ]
    def __str__(self):
        names = [x[0] for x in self._fields_]
        return ", ".join( ("%s: %s" % (key, getattr(self, key)))
                          for key in names )
    def reset(self):
        names = [x[0] for x in self._fields_]
        for name in names:
            setattr(self, name, 0)

    def collect(self, other):
        names = [x[0] for x in self._fields_]
        for name in names:
            my_val = getattr(self, name)
            other_val = getattr(other, name)
            setattr(self, name, my_val + other_val)
                            
        
class Cells(object):
    def __init__(self, binning, num_particles):
        _make_cells = clib.Cells_new
        _make_cells.restype = c.c_void_p
        maker = lambda: c.cast(_make_cells(binning, num_particles), c.c_void_p)
        self.pointer = maker()
    
    def __del__(self):
        clib.Cells_free(self.pointer)


class Simulation(Visualiser.simulation.Simulation):

    def __init__(self, num_particles, dt, cutoff, gamma, binning=1):
        super(Simulation, self).__init__()
        self.ArrayType = Particle * num_particles
        self.particles = self.ArrayType()
        self.num_particles = c.c_int(num_particles)
        self.dt = c.c_double(dt)
        self.gamma = c.c_double(gamma)
        self.cutoff = c.c_double(cutoff)
        self.binning = binning
        self.num_cells = binning**3
        self.cells = Cells(binning, self.num_particles)
        self.pos_ArrayType = c.c_double * (3 * num_particles)
        self.positions = self.pos_ArrayType()
        self.stats = Stats()
        self.last_coloring = 0
        self.is_coloring = False
        self.point_size = POINT_SIZE
        self.observable_filename = "davis_observables.csv"
        self.has_new_NN_data = False

    def timestep(self):
        if self.has_new_NN_data:
            self.write_observables()
            self.has_new_NN_data = False
        self.integrate()

    def integrate(self):
        self.stats.reset()
        clib.dvs_advance(self.num_particles, self.particles, self.dt)
        clib.dvs_populate_cells(self.num_particles, self.particles,
                             self.cells.pointer)
        clib.dvs_calc_forces(self.particles, 
                             self.cells.pointer, 0, self.num_cells,
                             self.cutoff, self.gamma, 
                             c.pointer(self.stats))
        clib.dvs_correct(self.num_particles, self.particles, self.dt)

    def do_print_stats(self):
        s = self.stats
        print str(s),
        print "N**2/2:", self.num_particles.value**2/2,
        print "WW/(N**2/2):", s.ww_counter / (self.num_particles.value**2/2.0),
        print "WW/Real_WW:", s.ww_counter / (1.0*s.real_ww_counter + 0.00001),
        print "Real_WW/N", s.real_ww_counter / (self.num_particles.value*1.0)

    def retrieve_visual_data(self):
        # retrieve data from simulation kernel for color_data and position_data
        self.get_3dPositions()
        a = np.frombuffer(self.positions)    
        self.position_data = a.reshape(-1, 3)        
        self.has_data = True
        if (not self.is_coloring 
            and time.time() - self.last_coloring > COLORING_INTERVAL):
            self.is_coloring = True
            thread.start_new_thread(self.do_NN_coloring, ())

    def get_3dPositions(self):
        clib.dvs_visualise_positions(self.num_particles, self.particles,
                                     self.positions)

    def render(self):
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)     
        glutSolidSphere(.99,128,128)
        glDisable(GL_LIGHTING)
        glColor3f(0., 1.0, 0.2)
        super(Simulation, self).render()

    def do_NN_coloring(self):
        if not self.has_data:
            return

        hull = ConvexHull(self.position_data) 
        # TODO the following could be made more efficient?
        county = np.zeros(len(hull.points), int)
        for v in hull.simplices:
            for i in v:
                county[i] += 1
        data = np.zeros(shape=(len(county), 3), dtype='float32')
        self.count6 = self.countless = self.countmore = 0
        for i, c in enumerate(county):
            if c == 6:
                data[i,0] = 0.0
                data[i,1] = 1.0
                data[i,2] = 0.0
                self.count6 += 1
            elif c < 6:
                data[i,0] = 0.0
                data[i,1] = 0.0
                data[i,2] = 1.0
                self.countless += 1
            else:
                data[i,0] = 1.0
                data[i,1] = 0.0
                data[i,2] = 0.0
                self.countmore += 1
        self.color_data = data
        self.has_color_data = True
        self.is_coloring = False
        self.last_coloring = time.time()
        self.has_new_NN_data = True
    
    def write_observables(self):
        # TODO with csv module
        with open(self.observable_filename, "a") as f:
            f.write("%d, %f, %d, %d, %d\n" % 
                    (self.steps, self.stats.E_pot, 
                     self.count6, self.countless, self.countmore))



class ParallelSimulation(Simulation):
    def __init__(self, *args, **kwargs):
        self.num_workers = kwargs.setdefault('num_workers', 1)
        print "Parallel with", self.num_workers, "workers."
        del kwargs['num_workers']
        assert(self.num_workers > 1)
        super(ParallelSimulation, self).__init__(*args, **kwargs)

        self.queue = Queue.Queue()
        self.particles_pw = [self.particles]
        self.stats_pw = [self.stats]
        for i in range(self.num_workers - 1):
            self.particles_pw.append(self.ArrayType())
            self.stats_pw.append(Stats())

        delta = int(math.ceil(self.num_cells / (1.0*self.num_workers)))
        self.intervals = []
        cell0 = 0
        for i in range(self.num_workers):
            cell1 = min(self.num_cells, cell0 + delta)
            self.intervals.append((cell0, cell1))
            t = threading.Thread(target=self.worker_starter(i))
            t.daemon = True
            t.start()
            cell0 = cell1
        print self.intervals

    def worker_starter(self, i):
        # separate method for closure handling
        def start():
            print "started worker", i
            while True:
                start_cell, end_cell = self.queue.get()
                clib.dvs_calc_forces(self.particles_pw[i], self.cells.pointer, 
                                     start_cell, end_cell,
                                     self.cutoff, self.gamma, 
                                     c.pointer(self.stats_pw[i]))
                self.queue.task_done()
        return start

    def integrate(self):
        clib.dvs_advance(self.num_particles, self.particles, self.dt)
        clib.dvs_populate_cells(self.num_particles, self.particles,
                             self.cells.pointer)
        for i in range(1, self.num_workers):
            clib.dvs_copy_particles(self.num_particles, 
                                    self.particles, self.particles_pw[i])
        for i in range(self.num_workers):
            self.stats_pw[i].reset()
        for i in self.intervals:
            self.queue.put(i)
        self.queue.join()
        for i in range(1, self.num_workers):
            clib.dvs_collect_forces(self.num_particles, 
                                    self.particles, self.particles_pw[i])
            self.stats.collect(self.stats_pw[i])

        clib.dvs_correct(self.num_particles, self.particles, self.dt)


def parallel_simulation(num_workers):
    def maker(*args, **kwargs):
        kwargs['num_workers'] = num_workers
        return ParallelSimulation(*args, **kwargs)
    return maker


# from https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere/26127012#26127012
def fibonacci_sphere(samples=1,randomize=True):
    rnd = 1.
    if randomize:
        rnd = random.random() * samples

    points = []
    offset = 2./samples
    increment = math.pi * (3. - math.sqrt(5.));

    for i in range(samples):
        y = ((i * offset) - 1) + (offset / 2);
        r = math.sqrt(1 - pow(y,2))

        phi = ((i + rnd) % samples) * increment

        x = math.cos(phi) * r
        z = math.sin(phi) * r

        points.append([x,y,z])

    return points
    
def demo(N=100, dt=0.0001, cutoff=0.25, gamma=0.01, 
         binning=1, vel_scale=0.0, simu_type=Simulation):
    s = simu_type(N, dt=dt, cutoff=cutoff, gamma=gamma, binning=binning)
    points = fibonacci_sphere(int(N*1.0))
    for p, point in it.izip(s.particles, points):
        p.r.x, p.r.y, p.r.z = point
        # # velocity
        # if abs(p.r.z) > abs(p.r.x):
        #     p.v.x = (random.random() - 0.5) * vel_scale
        #     p.v.y = (random.random() - 0.5) * vel_scale
        #     p.v.z = - (p.r.x*p.v.x + p.r.y*p.v.y) / p.r.z
        # else:
        #     p.v.z = (random.random() - 0.5) * vel_scale
        #     p.v.y = (random.random() - 0.5) * vel_scale
        #     p.v.x = - (p.r.z*p.v.z + p.r.y*p.v.y) / p.r.x
    # TODO: Reset all velocities to zero angular momentum?
    return s

        
    
if __name__ == '__main__':
    N = NUM_PARTICLES
    A_perPartcile = 4*math.pi / N
    cutoff = math.sqrt(A_perPartcile) * 2.5
    print "cutoff", cutoff
    # binning with regard to R=1 (so Lx=Ly=Lz=2.0)
    binning = max(1, int(2.0/cutoff))
    vel_scale = 15.0*cutoff    # purely handwaving
    if PARALLEL > 1: 
        simu_type = parallel_simulation(PARALLEL)
        print "Parallel simulation"
    else:
        simu_type = Simulation
        print "Single threaded simulation"
    simu = demo(N=N, dt=0.0001, cutoff=cutoff, gamma=0.01, 
                vel_scale=cutoff, binning=binning,
                simu_type=simu_type)
    Window(simu, "Davis Sphere Simulation, N=%d" % N)
        

