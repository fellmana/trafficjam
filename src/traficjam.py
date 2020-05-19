#!/usr/bin/python3
"""
=========================
       MONTE CARLO,
      Aslak Fellman,
      Final project,
      Traffic jam   
          2020
=========================
"""
import numpy as np
import argparse
import scipy.stats as stats
import itertools

#------------------------------------------

def argumentparser():
    """
    --------------------------
    Basic argument parser for
    command line arguments
    --------------------------
    :input:   None
    :output:  dictionary containing
              arguments and values
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-l',type=int,default=200,
                        help="Length of road, default=300")
    parser.add_argument('-n',type=int,default=100,
                        help="Number of cars, default=150") 
    parser.add_argument('-p',action='store_true',default=False,
                        help="Plot situation,default=False")
    parser.add_argument('--rolling',action='store_true',default=False,
                        help="Dissalble rolling in visualmode,default=False")
    parser.add_argument('--random',action='store_true',default=False,
                        help="Random positions,default=False")
    parser.add_argument('-pr',type=int,default=100,
                        help="Print/plot rate, default=100")
    parser.add_argument('-d',type=int,default=None,
                        help="Create initial array with density between [1]/10-[9]/10, default=None")
    parser.add_argument('-m','--max',type=int,default=10**5,
                        help="Maximum monte-carlo steps,default=10**6")
    parser.add_argument('-sl','--speedlimit',type=int,default=None,
                        help="Impose speed limit on cars, Default=None")
    parser.add_argument('--debug',action='store_true',default=False,
                        help="print debug options")
    return parser.parse_args()
#-------------------------------------------
# SIMULATION MAIN METHODS
#-------------------------------------------

def run_simulation(road,args):
    pass

#-------------------------------------------

def run_visual_simulation(road,args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    fig = plt.figure()
    ax  = fig.add_subplot(121, polar=True)
    ax2 = fig.add_subplot(122)
    plt.setp(ax.yaxis.get_ticklabels(), visible=False) 
    ax.set_yticklabels([])
    ax2.set_xlabel("MC steps",fontsize=20)
    ax2.set_ylabel("Mean queue length",fontsize=20)
    ax2.set_xlim(0, args.max)
    ax2.set_ylim(0.0, 2*(args.l/args.n))
   
    circle, = ax.plot([],[],marker='s',color="k",linestyle="",ms=7)
    info,   = ax2.plot([],[],color='r')
    index = np.array([i for i in range(road.size)])
    title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
    steps = 0
    steps_taken = np.array([])
    mqls        = np.array([])

    def animate(i):
        nonlocal road
        nonlocal steps
        nonlocal steps_taken, mqls
        for _ in range(args.pr):
            road = monte_carlo(road,index,args)
            steps += 1
        if args.rolling:
            road = np.roll(road,1)
        mask  = road == 1
        steps_taken = np.append(steps_taken,steps)
        mqls        = np.append(mqls,mean_queue_length(road))
        circle.set_data(index[mask]*(2*np.pi/len(index)),road[mask])
        info.set_data(steps_taken,mqls)
        
        
        title.set_text("FRAME " + str(i))
        return circle,info,title,

    anim = FuncAnimation(fig, animate,frames=100,interval=200, blit=True)
    plt.show()

#-------------------------------------------
# MONTE CARLO SIMULATION METHOD
#-------------------------------------------

def monte_carlo(road,index,args):
    """
    :Description: 
    function that handles a monte carlo step
    :input: road = np.array containing information
    """
    mask = road == 1
    nonzeros = index[mask]

    rand = np.random.randint(len(nonzeros))
    speed = abs(nonzeros[(rand+ 1) % len(nonzeros)]
                - nonzeros[rand]
                - len(road)*round((nonzeros[(rand+ 1) % len(nonzeros)]-nonzeros[rand] )/len(road)) - 1)
    
    if args.speedlimit != None and speed > args.speedlimit:
        speed = args.speedlimit

    if speed == 0:
        return road

    if np.random.rand() < np.exp(-(1.0/speed)):
        road[nonzeros[rand]] = 0
        road[(nonzeros[rand] + 1) % len(road)] = 1
        return road

    return road

#-------------------------------------------
# CONVENIENCE FUNCTIONS
#-------------------------------------------

def mean_queue_length(road):
    return np.mean([sum(g) for b, g in itertools.groupby(road) if b])


#-------------------------------------------
# LATTICE CREATION OPTIONS BELOW
#-------------------------------------------

def lattice_even(L,N):
    road = np.zeros(L)
    road[::2] = 1
    return road

#-------------------------------------------

def lattice_random(L,N):
    road = np.zeros(L)
    mask = np.zeros(L);mask[0:N] = 1
    np.random.shuffle(mask)
    road[mask.astype(np.bool)] = 1
    return road

def lattice_density_10(L,denominator):
    stencil = np.zeros(10)
    stencil[:denominator] = 1
    road = np.tile(stencil,round(L/10))
    remainder = L % 10
    road = np.append(road,stencil[0:remainder])
    return road


#-------------------------------------------


#-------------------------------------------
# IF PROGRAM IS RUN AS __main__
#-------------------------------------------
if __name__ == "__main__":
    # Parse commandline argumet
    args = argumentparser()
    #Create initial condition
    if args.random:
        road = lattice_random(args.l,args.n)
    elif args.d != None:
        road = lattice_density_10(args.l,args.d)
        print(road)
    else:
        road = lattice_even(args.l,args.n)
    
    if args.p:
        run_visual_simulation(road,args)
    else:
        run_simulation(road,args)
    