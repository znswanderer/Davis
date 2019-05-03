import math
import davis

import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1" 

import davis_tf
from gui import Window


if __name__ == '__main__':
    NUM_PARTICLES = 1000
    DT = 0.0001
    GAMMA = 0.01

    A_perPartcile = 4*math.pi / NUM_PARTICLES
    cutoff = math.sqrt(A_perPartcile) * 2
    #cutoff = 0.1
    print(("cutoff", cutoff))
    binning = max(1, int(2.0/cutoff))
    print("binning", binning)

    world = davis_tf.World(NUM_PARTICLES, DT, GAMMA, cutoff, binning)
    world.build_graph()
    simu = davis.Simulation(world)
        
    Window(simu, "Davis Sphere Simulation, N=%d, dt=%.0le" % (NUM_PARTICLES, DT))
