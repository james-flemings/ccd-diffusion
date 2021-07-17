#!/usr/bin/env python3
from scipy.stats import expon
from random import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import multiprocessing as mp
import math

# Constants
T = 163.15 # temperature in kelvin (-110 Celsius)
k = 8.617333e-5 # Boltzman constant in eV K^-1

numOfElect = 10
initial = [-3500, -1500, 0, 1500, 3500] # vary start position in z direction
#initial = [0] # no varying start position

# This function produces a graph of the boltzman distribution
# using our temperature value and boltzmann contstant. 
# It also uses the mean free path values from a text file
def boltzmann(mfp_table):
    energy = []
    r_exp = []
    for i in range(1000):
        prob = random()
        e = expon(0, k*T).rvs()
        mfp = 0

        if e < 0.1:
            energy.append(round(e, 2))
        else:
            energy.append(round(e, 1))
        if energy[i] == 0:
            energy[i] = 0.01

        if prob <= mfp_table[energy[i]][1]:
            mfp = mfp_table[energy[i]][0]
        else:
            mfp = mfp_table[energy[i]][2]
        r_exp.append(expon(0, mfp).rvs())

    plt.subplot(1, 2, 1)
    plt.hist(r_exp, bins='auto')
    plt.title("Path length distribution")
    plt.xlabel("Path length (nm)")
    plt.ylabel("Entries")

    plt.subplot(1, 2, 2)
    plt.hist(r_exp, bins='auto', log=True)
    plt.title("Path length distribution (log scale)")
    plt.xlabel("Path length (nm)")
    plt.ylabel("Entries (log)")
    plt.show()

    plt.hist(energy)
    plt.title("Boltzmann Distribution")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Entries")
    plt.show()

def simulate(z0, lowerBoundary, upperBoundary, mfp_table, x = [], y = [], z = []):
    i = 1
    dx = []
    dy = []
    dz = []
    dx.append(0)
    dy.append(0)
    dz.append(z0)
    shift = 0 # This is used to include points when a electron hits the upper boundary
    exit_x = 0
    exit_y = 0
    exit_r = 0
    while True:
        energy = 0
        prob = random()
        e = expon(0, k*T).rvs() # Boltzman Distribution
        if e < 0.1:
            energy = round(e, 2)
        else:
            energy = round(e, 1)
        if energy == 0: # No electrons with 0 energy
            energy = 0.01

        # probability that electron emits or absorps a phonon, which determines mfp
        if prob <= mfp_table[energy][1]:
            mfp = mfp_table[energy][0]
        else:
            mfp = mfp_table[energy][2]

        r_exp = expon(0, mfp) # exponential distribution with mean equal to mean free path
        theta = np.pi * random() # theta between 0 to 2 * pi
        phi = 2 * np.pi * random() # phi between 0 and pi
        r = r_exp.rvs() # choosing a magnitude from exponetial distribution
        z_check = r * np.cos(theta) + dz[i+shift-1]

        if z_check > upperBoundary: # electron must reflect if it reaches the upper boundary
            # The distance from current position to upper boundary
            r2 = (upperBoundary - dz[i+shift-1])/np.cos(theta)
            dx.append(r2 * np.sin(theta) * np.cos(phi) + dx[i+shift-1])
            dy.append(r2 * np.sin(theta) * np.sin(phi) + dy[i+shift-1])
            dz.append(r2 * np.cos(theta) + dz[i+shift-1])
            shift +=1
            r = r - r2 # remaining distance to travel after reflection
            theta = np.pi - theta
            z_check = r * np.cos(theta) + dz[i+shift-1]

        if z_check >= lowerBoundary:
            dx.append(r * np.sin(theta) * np.cos(phi) + dx[i+shift-1])
            dy.append(r * np.sin(theta) * np.sin(phi) + dy[i+shift-1])
            dz.append(z_check)

        else: # electron has reached the exit boundary
            r2 = (lowerBoundary - dz[i+shift-1])/np.cos(theta)
            dx.append(r2 * np.sin(theta) * np.cos(phi) + dx[i+shift-1])
            dy.append(r2 * np.sin(theta) * np.sin(phi) + dy[i+shift-1])
            dz.append(r2 * np.cos(theta) + dz[i+shift-1])
            exit_x = dx[i+shift]/8000
            exit_y = dy[i+shift]/8000
            exit_r = np.sqrt(exit_x**2 + exit_y**2)
            break
        i += 1

    x.append(dx)
    y.append(dy)
    z.append(dz)
    return exit_x, exit_y, exit_r

def main():
    mfp_table = {} # mfp and probabilities indexed by electron energy
    l = [0] * 5
    f = open('mfp.txt', 'r')
    for line in f:
        sp = line.split()
        j = 0
        for i in sp:
            l[j] = float(i)
            j += 1
        mfp_table[l[0]] = l[1:5]

    upperBoundary = 4000
    lowerBoundary = -4000
    x_hist = []
    y_hist = []
    r_hist = []

    for z0 in initial:
        exit_x = []
        exit_y = []
        exit_r = []
        pool = mp.Pool(mp.cpu_count()) # iniatlize multi-processing
        results = [pool.apply_async(simulate, args=(z0, lowerBoundary, upperBoundary,
            mfp_table)) for j in range(numOfElect)]

        pool.close()
        pool.join()

        for r in results:
            exit_x.append(r.get()[0])
            exit_y.append(r.get()[1])
            exit_r.append(r.get()[2])

        x_hist.append(exit_x)
        y_hist.append(exit_y)
        r_hist.append(exit_r)

        mean_x = round(sum(exit_x)/len(exit_x), 4)
        mean_y = round(sum(exit_y)/len(exit_y), 4)
        mean_r = round(sum(exit_r)/len(exit_r), 4)

        square_x = 0
        square_y = 0
        square_r = 0

        for i in range(len(exit_r)):
            square_x += (exit_x[i]**2)
            square_y += (exit_y[i]**2)
            square_r += (exit_r[i]**2)

        sqr_mean_x = (square_x / float(len(exit_x)))
        sqr_mean_y = (square_y / float(len(exit_y)))
        sqr_mean_r = (square_r / float(len(exit_r)))

        root_x = round(math.sqrt(sqr_mean_x), 4)
        root_y = round(math.sqrt(sqr_mean_y), 4)
        root_r = round(math.sqrt(sqr_mean_r), 4)

        print("z0 = {}\t mean_x = {}\t mean_y = {}\t mean_r = {}".format(
            z0, mean_x, mean_y, mean_r))
        print("z0 = {}\t RMS_x = {}\t RMS_y = {}\t RMS_r = {}\n".format(
            z0, root_x, root_y, root_r))

    #boltzmann(mfp_table)

    # r exit coordinates histogram
    colors = ['b', 'g', 'r', 'c', 'm']
    for i in range(len(initial)):
        plt.hist(r_hist[i],bins='auto', color=colors[i], label=str(initial[i])+" nm", histtype='step')
    plt.title("Exit coordinates for r")
    plt.xlabel("r/dff")
    plt.ylabel("Entries")
    plt.legend(loc='upper right')
    plt.show()

    # x exit coordinates histogram
    for i in range(len(initial)):
        plt.hist(x_hist[i],bins='auto', color=colors[i], label=str(initial[i])+" nm", histtype='step')
    plt.title("Exit coordinates for x")
    plt.xlabel("x/dff")
    plt.ylabel("Entries")
    plt.legend(loc="upper right")
    plt.show()

    # y exit coordinates histogram
    for i in range(len(initial)):
        plt.hist(y_hist[i],bins='auto', color=colors[i], label=str(initial[i])+" nm", histtype='step')
    plt.title("Exit coordinates for y")
    plt.xlabel("y/dff")
    plt.ylabel("Entries")
    plt.legend(loc="upper right")
    plt.show()

    x, y , z = [], [], []
    simulate(0, lowerBoundary, upperBoundary, mfp_table, x, y, z)

    # Animation plot y vs z
    fig = plt.figure()
    ax = plt.axes(xlim=(min(y[0])-10, max(y[0])+10), ylim=(lowerBoundary-100, upperBoundary+100))
    line, = ax.plot([], [], lw = 2)

    def init():
        line.set_data([], [])
        return line, 

    def iterr(i, y, z):
        line.set_data(y[0:i], z[0:i])
        return line, 

    anim = animation.FuncAnimation(fig, iterr, frames=(len(y[0])+2),
                                    init_func=init, fargs=(y[0],z[0]), interval=600,
                                    repeat=True, blit=True)
    plt.show()
    
    # x vs y and y vs z displacement graphs and scatter plots
    plt.subplot(2, 2, 1)
    plt.plot(x[0], y[0], color = 'red')
    plt.title("x vs y displacement")
    plt.xlabel("x (nm)")
    plt.ylabel("y (nm)")

    plt.subplot(2, 2, 2)
    plt.plot(y[0], z[0])
    plt.title("y vs z displacement")
    plt.xlabel("y (nm)")
    plt.ylabel("z (nm)")

    plt.subplot(2, 2, 3)
    plt.scatter(x[0], y[0], s = 2, color = 'red')
    plt.title("x vs y displacement")
    plt.xlabel("x (nm)")
    plt.ylabel("y (nm)")

    plt.subplot(2, 2, 4)
    plt.scatter(y[0], z[0], s = 2)
    plt.title("y vs z displacement")
    plt.xlabel("y (nm)")
    plt.ylabel("z (nm)")

    plt.show()

if __name__=="__main__":
    main()
