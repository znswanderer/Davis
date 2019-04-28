"""
This is a little script to measure the time complexity behaviour
for varying number of particles as described in the README.md.
"""

import math
import timeit
import davis_c
import time
import sys


if __name__ == '__main__':
    NUM_PARTICLES = (
           list(range(1000, 20000, 500)) 
         + list(range(20000, 50000, 1000))
         + list(range(50000, 60000, 5000))
         + list(range(60000, 110000, 10000))
         + list(range(150000, 350000, 50000))
         + [400000, 500000, 750000, 1000000,1250000, 1500000, 2000000, 3000000]
    )
    max_particles = max(NUM_PARTICLES)
    DT = 0.00001
    GAMMA = 0.01
    STEPS = 4000
    ORIG_STEPS = 7500

    print(NUM_PARTICLES)

    for world_class in (davis_c.World,):
        print(world_class.__name__)
        print("started...")
        results = []
        try:
            for n in NUM_PARTICLES:
                # Here we are using a flexible cutoff.
                cutoff = math.sqrt(4*math.pi / n) * 1.8
                binning = max(1, int(2.0/cutoff))
                world = world_class(n, DT, GAMMA, cutoff, binning, num_workers=PARALLEL)
                # We vary the number of steps with the number particles, such
                # that for few particles we run more steps in order to get a better
                # statistics and for a large number of particles we use smaller steps
                # to finish in reasonable run-times.
                # For the output all times will be scaled to ORIG_STEPS.
                steps = int(STEPS * math.sqrt(max_particles / n))
                t = timeit.timeit("for i in range(steps):\n  world.timestep()", globals=globals(), number=1)
                world.stop()
                print("n:", n, "time:", t, "cutoff", cutoff, "steps", steps, "time", time.ctime(time.time()))
                results.append((n, t, cutoff, steps))
        except KeyboardInterrupt:
            print("Stopping.")

        print("\n" + 80*"-" + "\n")
        print("world class:", world_class.__name__)
        print("cutoff = math.sqrt(4*math.pi / n) * 1.8")
        print("STEPS:", STEPS)
        print("DT:", DT)
        print("GAMMA", GAMMA, "\n")
        for n, t, cutoff, steps in results:
            print("n:", n, "scaled_time:", t * (ORIG_STEPS / steps), "cutoff", cutoff, "steps", steps)
        print("\n" + 80*"-" + "\n")

