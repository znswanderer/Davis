"""
Python part for
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


import os
import sys
import time

import numpy as np
import OpenGL.GL as GL
import OpenGL.GLUT as GLUT
import scipy.spatial

import _thread
import simulation

COLORING_INTERVAL = 1.0      # in seconds
POINT_SIZE = 3

GREEN = np.array([0, 1, 0], dtype='float32')
RED = np.array([1, 0, 0], dtype='float32')
BLUE = np.array([0, 0, 1], dtype='float32')
YELLOW = np.array([1, 1, 0], dtype='float32')

class Simulation(simulation.Simulation):

    def __init__(self, world):
        super(Simulation, self).__init__()
        self.world = world
        self.last_coloring = 0
        self.is_coloring = False
        self.point_size = POINT_SIZE
        self.simulation_running = True   
        self.position_data = None


    def timestep(self):
        self.integrate()

    def integrate(self):
        self.world.timestep()

    def retrieve_visual_data(self):
        """Retrieve data from simulation kernel for color_data and position_data"""
        self.position_data = self.world.get_3dPositions()
        if (not self.is_coloring and time.time() - self.last_coloring > COLORING_INTERVAL):
            self.is_coloring = True
            _thread.start_new_thread(self.do_NN_coloring, ())

    def render(self):
        GL.glEnable(GL.GL_LIGHTING)
        GL.glEnable(GL.GL_LIGHT0)
        GLUT.glutSolidSphere(.99,128,128)
        GL.glDisable(GL.GL_LIGHTING)
        GL.glColor3f(0., 1.0, 0.2)
        super(Simulation, self).render()

        if not self.position_data is None:
            i0 = int(len(self.position_data) / 2)
            x, y, z = self.position_data[i0]
            GL.glPushMatrix()
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
            GL.glTranslated(x, y, z)
            GL.glColor4f(1, 1, 0, 0.5)
            GLUT.glutSolidSphere(self.world.cutoff, 16, 16)
            GL.glPopMatrix()

    def do_NN_coloring(self):
        """Find the color for a particle based on the number of nearest neighbours."""
        if self.position_data is None:
            return

        hull = scipy.spatial.ConvexHull(self.position_data)
        a = hull.simplices.reshape(-1)
        counts = np.bincount(a)

        data = ( (counts == 6)[:, np.newaxis] * GREEN
               + (counts < 6)[:, np.newaxis] * BLUE
               + (counts > 6)[:, np.newaxis] * RED )

        i0 = int(len(data) / 2)
        data[i0] = YELLOW

        self.count6 = np.sum(counts == 6)
        self.countless = np.sum(counts < 6)
        self.countmore = np.sum(counts > 6)

        self.color_data = data
        self.has_color_data = True
        self.is_coloring = False
        self.last_coloring = time.time()



