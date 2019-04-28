import math
import davis
import davis_c
from gui import Window


if __name__ == '__main__':
    PARALLEL = 4                 # number of parallel workers
    NUM_PARTICLES = 20000
    DT = 0.0001
    GAMMA = 0.01

    A_perPartcile = 4*math.pi / NUM_PARTICLES
    cutoff = math.sqrt(A_perPartcile) * 2
    #cutoff = 0.1
    print(("cutoff", cutoff))
    binning = max(1, int(2.0/cutoff))
    print("binning", binning)

    world = davis_c.ParallelWorld(NUM_PARTICLES, DT, GAMMA, cutoff, binning, num_workers=PARALLEL)
    simu = davis.Simulation(world)
        
    Window(simu, "Davis Sphere Simulation, N=%d, dt=%.0le" % (NUM_PARTICLES, DT))
