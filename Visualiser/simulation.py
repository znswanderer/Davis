import threading
import time
import OpenGL.GL as GL
import OpenGL.GLU as GLU
import numpy as np

from gui import Window

class Simulation(threading.Thread):
    LOCK = threading.Lock()
    DATA_LOCK = threading.Lock()
    PRINT_INTERVAL = 2

    def __init__(self):
        super(Simulation, self).__init__()
        self.needs_new_data = False
        self.has_data = False
        self.has_color_data = False
        self.last_print = 0
        self.color_data = None
        self.position_data = None
        self.point_size = 3
        self.simulation_running = False
        self.print_stats = True
        self.daemon = True
        self.steps = 0

    def get_color_data(self):
        return self.color_data

    def get_position_data(self):
        return self.position_data

    def retrieve_visual_data(self):
        # retrieve data from simulation kernel
        # for color_data and position_data
        pass

    def save(self):
        # save state of simulation kernel
        print "Saving..."

    def timestep(self):
        # advance simulation kernel
        time.sleep(0.01)    # overwrite in sub-class

    def render(self):
        # called from main thread (OpenGL)
        self.DATA_LOCK.acquire()
        if self.has_data:
            positions = self.get_position_data()
            colors = self.get_color_data()
            indices = np.arange(len(positions), dtype='i')
            GL.glEnable(GL.GL_POINT_SMOOTH)
            GL.glPointSize(self.point_size)
            GL.glEnable(GL.GL_BLEND)
            GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
            GL.glVertexPointerd(positions)
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            if self.has_color_data:
                GL.glColorPointerd(colors)
                GL.glEnableClientState(GL.GL_COLOR_ARRAY)
            GL.glDrawArrays(GL.GL_POINTS, 0, len(positions))
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
            if self.has_color_data:
                GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisable(GL.GL_BLEND)		
        self.DATA_LOCK.release()

        self.LOCK.acquire()
        self.needs_new_data = True
        self.LOCK.release()

    def start_run(self):
        self.simulation_running = not self.simulation_running

    def handle_key(self, key):
        pass

    def do_print_stats(self):
        print

    def run(self):
        last_steps = self.steps
        data_fetches = 0
        while True:
            self.LOCK.acquire()
            get_new_data = self.needs_new_data
            self.LOCK.release()
            if get_new_data:
                self.DATA_LOCK.acquire()
                self.retrieve_visual_data()
                self.DATA_LOCK.release()
                self.LOCK.acquire()
                self.needs_new_data = False
                self.LOCK.release()
                data_fetches += 1

            if self.simulation_running:
                self.timestep()
                self.steps += 1
            else:
                time.sleep(0.01)
            
            if time.time() > self.last_print + self.PRINT_INTERVAL:
                if self.print_stats:
                    print "steps", self.steps-last_steps, 
                    print "data_fetches", data_fetches,
                    self.do_print_stats()
                self.last_print = time.time()
                last_steps = self.steps
                data_fetches = 0


