import sys
import matplotlib.pyplot as plt

class Result(object):
    pass

def read_file(fname):
    with open(fname) as f:
        lines = f.readlines()
    res = []
    for line in lines:
        r = Result()
        cols = [x.strip() for x in line.split(",")]
        r.steps = int(cols[0])
        r.E_pot = float(cols[1])
        r.num_6NN = int(cols[2])
        if r.steps > 0:
            res.append(r)
    return res

def main(fname):
    # http://mple.m-artwork.eu/tutorial
    res = read_file(fname)
    steps = [x.steps for x in res]
    nums_6NN = [x.num_6NN for x in res]
    #define plot size in inches (width, height) & resolution(DPI)
    fig = plt.figure(figsize=(5, 3), dpi=100)
    #define font size
    plt.rc("font", size=10)

    plt.plot(steps, nums_6NN, 'ro-', fillstyle='none', markersize=3)
    plt.ylabel("# of particles with 6 NN")
    plt.xlabel("simulation steps")

    #adjust plot
    plt.subplots_adjust(left=0.15,bottom=0.15)
    
    fig.savefig('figure_1.png') 
    
    

if __name__ == '__main__':
    main(sys.argv[1])
