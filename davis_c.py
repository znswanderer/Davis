"""
Interface module for the C-version of davis.
Molecular Dynamics Simulation on a Sphere

Please see README.md for details.

MIT license

Copyright(c) 2015 - 2019 Tim Scheffler

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
import queue
import threading
import fibonacci_sphere
import glob
import os

possible_clibs = glob.glob("davis_clib*")
if len(possible_clibs) == 1:
    lib_name = possible_clibs[0]
    print("Loading %s" % lib_name)
    clib = c.CDLL(os.path.join(".", lib_name))
elif len(possible_clibs) == 0:
    print("Cannot find the davis clib part. Have you run the 'python setup.py build_ext -i' step?")
    exit(-1)
else:
    print("Found %d possible davis_clib modules (%s). Must be unique!" % (
            len(possible_clibs), ", ".join(possible_clibs)))
    exit(-1)

class Vec(c.Structure):
    _fields_ = [("x", c.c_double), ("y", c.c_double), ("z", c.c_double)]

    def __str__(self):
        return "(%f, %f, %f)" % (self.x, self.y, self.z)

    def magnitude2(self):
        return self.x**2 + self.y**2 + self.z**2


class Particle(c.Structure):
    _fields_ = [("r", Vec), ("v", Vec), ("a", Vec), ("next", c.c_long)]


class Stats(c.Structure):
    _fields_ = [("ww_counter", c.c_long), ("real_ww_counter", c.c_long),
                ("E_pot", c.c_double)]

    def __str__(self):
        names = [x[0] for x in self._fields_]
        return ", ".join( ("%s: %s" % (key, getattr(self, key))) for key in names)

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


class World:
    """Single-threaded simulation with force calculations based on cells."""

    def __init__(self, num_particles, dt=0.1, gamma=0.0, cutoff=0.0, binning=None, *args, **kwargs):
        self.ArrayType = Particle * num_particles
        self.particles = self.ArrayType()
        self.num_particles = c.c_long(num_particles)
        self.dt = c.c_double(dt)
        self.gamma = c.c_double(gamma)
        self.cutoff = c.c_double(cutoff)
        self.binning = binning
        self.num_cells = binning**3
        self.cells = Cells(binning, self.num_particles)
        self.pos_ArrayType = c.c_double * (3 * num_particles)
        self.positions = self.pos_ArrayType()
        self.stats = Stats()

        points = fibonacci_sphere.create_points(num_particles)
        for p, point in zip(self.particles, points):
            p.r.x, p.r.y, p.r.z = point

    def timestep(self):
        clib.dvs_advance(self.num_particles, self.particles, self.dt)
        clib.dvs_populate_cells(self.num_particles, self.particles, self.cells.pointer)
        clib.dvs_calc_forces(self.particles, self.cells.pointer, 0, self.num_cells,
                             self.cutoff, self.gamma, c.pointer(self.stats))
        clib.dvs_correct(self.num_particles, self.particles, self.dt)

    def get_3dPositions(self):
        clib.dvs_visualise_positions(self.num_particles, self.particles, self.positions)
        a = np.frombuffer(self.positions)
        return a.reshape(-1, 3)
    
    def stop(self):
        pass


class BruteForceWorld(World):
    """Single-threaded simulation with O(N^2) force calculations"""

    i_0 = c.c_long(0)

    def timestep(self):
        clib.dvs_advance(self.num_particles, self.particles, self.dt)
        clib.dvs_calc_brute_forces(self.num_particles, self.particles,
                self.i_0, self.num_particles, self.cutoff, self.gamma, c.pointer(self.stats))
        clib.dvs_correct(self.num_particles, self.particles, self.dt)


class ParallelWorld(World):
    """Multi-threaded simulation with force calculations based on cells."""

    uses_cells = True

    def __init__(self, *args, **kwargs):
        self.num_workers = kwargs.setdefault('num_workers', 1)
        del kwargs['num_workers']
        assert self.num_workers > 1
        super().__init__(*args, **kwargs)

        self.queue = queue.Queue()
        self.particles_pw = [self.particles]
        #self.stats_pw = [self.stats]
        self.stats_pw = [Stats()]
        for i in range(self.num_workers - 1):
            self.particles_pw.append(self.ArrayType())
            self.stats_pw.append(Stats())

        delta = int(math.ceil(self.num_cells / (1.0*self.num_workers)))
        self.intervals = []
        cell0 = 0
        for i in range(self.num_workers):
            if i == self.num_workers - 1:
                cell1 = self.num_cells
            else:
                cell1 = cell0 + delta
            self.intervals.append((cell0, cell1))
            t = threading.Thread(target=self.worker_starter(i))
            t.daemon = True
            t.start()
            cell0 = cell1

    def stop(self):
        print("STOP")
        for _ in range(self.num_workers):
            self.queue.put('STOP')


    def worker_starter(self, i):
        # separate method for closure handling
        def start():
            print("started worker", i)
            while True:
                q = self.queue.get()
                if q == 'STOP':
                    break
                else:
                    start_cell, end_cell = q
                    clib.dvs_calc_forces(self.particles_pw[i], self.cells.pointer,
                                        start_cell, end_cell,
                                        self.cutoff, self.gamma,
                                        c.pointer(self.stats_pw[i]))
                self.queue.task_done()
        return start

    def timestep(self):
        clib.dvs_advance(self.num_particles, self.particles, self.dt)
        if self.uses_cells:
            clib.dvs_populate_cells(self.num_particles, self.particles, self.cells.pointer)
        for i in range(1, self.num_workers):
            clib.dvs_copy_particles(self.num_particles, self.particles, self.particles_pw[i])
        for i in self.intervals:
            self.queue.put(i)
        self.queue.join()
        for i in range(1, self.num_workers):
            clib.dvs_collect_forces(self.num_particles, self.particles, self.particles_pw[i])
        self.stats.reset()
        for i in range(0, self.num_workers):
            self.stats.collect(self.stats_pw[i])

        # ATTENTION!
        # Calculating E_kin below will decrease performance massively! 
        # (approx 10x! in ParallelSimulation and 2x in BruteParallelSimulation)
        # self.E_kin += sum(p.v.magnitude2() for p in self.particles) / 2.0

        clib.dvs_correct(self.num_particles, self.particles, self.dt)


class BruteParallelWorld(ParallelWorld):
    """Multi-threaded simulation with O(N^2) force calculations"""
    uses_cells = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        r"""
        For N particles we have to calculate N * (N - 1) / 2 interactions.
        So the "volume" of the integration is two-dimensional in the
        number of particles, O(N^2), as we move along the first particle
        axis (the "i"-axis) because of the i<j shortcut (actio=reactio),
        the number of interactions to calculate becomes smaller.

        If we want to distribute the number of interactions evenly 
        onto all worker threads, the number interval [i_s, i_e]
        for each worker will not be of the same length, with
        d_i = i_e - i_s we get the number of interactions 

            V = d_i * (d_i + 1) / 2 + (N - i_e - 1) * d_i

        (Index is zero-based, therefore \max(i) = N - 1. If we 
        have i_s = 0 and i_e = N - 1, the formula gives 
        V = N * (N - 1) / 2 as expected.)
        
        If we run the particles the outer index i \in [i_s, i_e] 
        and the inner index j \in [i, N].

        We would like to calculate 

            V_W = (N * (N - 1) / 2) / W

        Interactions on each of the W worker threads.

        So if we start the interval for a given thread at i_s, the 
        end index i_e for this thread is given by

            i_e = i_s + b - \sqrt{ b^2 - 2 V_W}
        with 
            b = N - i_s - 1/2
        
        """
        N = int(self.num_particles.value)
        V  = N * (N - 1) / 2
        V_W = V / self.num_workers
        self.intervals = []
        i_s = 0
        for n in range(self.num_workers):
            b = N - i_s - 0.5
            i_e = int(i_s + b - math.sqrt(b**2 - 2*V_W))
            if n == self.num_workers - 1:
                i_e = N
            self.intervals.append( (c.c_long(i_s), c.c_long(i_e)) )
            i_s = i_e

    def worker_starter(self, i):
        # separate method for closure handling
        def start():
            print(("started worker", i))
            while True:
                q = self.queue.get()
                if q == 'STOP':
                    break
                else:                
                    i_s, i_e = q
                    clib.dvs_calc_brute_forces(
                        self.num_particles,
                        self.particles_pw[i],
                        i_s, i_e,
                        self.cutoff, self.gamma,
                        c.pointer(self.stats_pw[i]))
                self.queue.task_done()
        return start

