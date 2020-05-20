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
class car:

    def __init__(self,pos,index, lane, track=False):
        self.position = pos
        self.index    = index
        self.lane     = lane
        self.speed    = 0
        if track:
            self.trace    = []
            self.time     = []

    def set_speed(self,val):
        self.speed = val

    def __str__(self):
        return "Car: " + str(self.index) + ",position: " + str(self.position) + ",speed: " + str(self.speed)


#------------------------------------------

class simulation:

    def __init__(self, array, track):
        self.cars = []
        self.R          = 0
        self.length = len(array)
        n = 0
        for i,elem in enumerate(array):
            if elem == 1:
                self.cars.append(car(i,n,0,track))
                n += 1
        self.num   = n
        self.cumulative = np.zeros(self.num)
        self.global_time = 0

    def get_array(self):
        arr = np.zeros(self.length)
        for elem in self.cars:
            arr[elem.position] = 1
        return arr

    def distance(self,a,b):
        return abs(self.cars[b].position - self.cars[a].position - self.length*round((self.cars[b].position - self.cars[a].position)/self.length))


    def get_random(self):
        rand = np.random.rand()*self.cumulative[-1]
        #print(rand)
        for i in range(self.num):
            if i == 0:
                if rand <= self.cumulative[0]:
                    return i
            else:
                if (self.cumulative[i-1] < rand <= self.cumulative[i]):
                    return i

    def calculate_speeds(self):
        for i in range(self.num):
            if i == self.num - 1:
                self.cars[i].set_speed(self.distance(i,0))
            else:
                self.cars[i].set_speed(self.distance(i,i+1))

    def calculate_cumulative(self):
        summa = 0
        for i in range(self.num):
            summa += self.cars[i].speed
            self.cumulative[i] = summa

    def move_car(self,a):
        if self.cars[a].speed > 1:
            if self.cars[a].position == self.length - 1:
                self.cars[a].position = 0
            else:
                self.cars[a].position += 1

    def update_history(self,t):
        for c in self.cars:
            c.trace.append(c.position)
            c.time.append(t)
    
    def update_time(self):
        self.global_time += -(np.log(np.random.rand()))/(self.cumulative[-1])

    def print_cars(self):
        print("Road: lengt: " + str(self.length) + ",cars: " + str(self.num))
        for c in self.cars:
            print(c)



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
    parser.add_argument('-l',type=int,default=300,
                        help="Length of road, default=300")
    parser.add_argument('-n',type=int,default=150,
                        help="Number of cars, default=150") 
    parser.add_argument('-p',action='store_true',default=False,
                        help="Plot situation,default=False")
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
    parser.add_argument('--follow',action="store_true",default=False,
                        help="Follow a single element during simulation")
    parser.add_argument('--debug',action='store_true',default=False,
                        help="print debug options")
    return parser.parse_args()
#-------------------------------------------
# SIMULATION MAIN METHODS
#-------------------------------------------

def run_simulation(r,args):

    lane = simulation(r)
    lane.print_cars()
    lane.calculate_speeds()
    lane.print_cars()
    lane.calculate_cumulative()
    print(lane.cumulative)
    c = lane.get_random()
    print(c)
    print(lane.get_array())

    lane.move_car(c)
    lane.calculate_speeds()
    lane.print_cars()
    print(lane.get_array())
    
    lane.calculate_cumulative()
    print(lane.cumulative)
#-------------------------------------------

def run_visual_simulation(r,args):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from collections import defaultdict
    
    index = np.array([i for i in range(r.size)])
    road = simulation(r,track=True)

    anim_running = True

    fig = plt.figure()
    ax  = fig.add_subplot(121, polar=True)
    ax2 = fig.add_subplot(122)
    plt.setp(ax.yaxis.get_ticklabels(), visible=False) 
    ax.set_yticklabels([])
    ax2.set_xlabel("Position on road [m]",fontsize=20)
    ax2.set_ylabel("time [s]",fontsize=20)
    line, = ax2.plot([], [], lw=2)
    ax2.set_xlim(0, args.l)
    ax2.set_ylim(0.0, 100)
   
    circle, = ax.plot([],[],marker='s',color="k",linestyle="",ms=7)
    title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")

    traces = []
    for j in range(args.n):
        trace = ax2.plot([],[],color="k",marker="s",ms=2,linestyle="")[0]
        traces.append(trace)

    steps = 0
    steps_taken = np.array([])
    mqls        = np.array([])

    def onClick(event):
        nonlocal anim_running
        if anim_running:
            anim.event_source.stop()
            anim_running = False
        else:
            anim.event_source.start()
            anim_running = True   

    def animate(i):
        nonlocal road
        nonlocal steps
        nonlocal steps_taken, mqls
        for _ in range(args.pr):
            road.calculate_speeds()
            road.calculate_cumulative()
            rand_car_index = road.get_random() 
            road.move_car(rand_car_index)
            road.update_time()

        road.update_history(road.global_time)
        r = road.get_array() 
        mask  = r == 1
        steps_taken = np.append(steps_taken,steps)
        mqls        = np.append(mqls,mean_queue_length(r))
        circle.set_data(index[mask]*(2*np.pi/len(index)),r[mask])
        title.set_text("FRAME " + str(i)) 

        for tn,trace in enumerate(traces):
            #print(tn,cars[tn][0],cars[tn][1])
            trace.set_data((road.cars[tn].trace),road.cars[tn].time)

        if i == 500:
            anim.event_source.stop()
        return circle,traces,title,
        #return traces

    fig.canvas.mpl_connect('button_press_event', onClick)
    anim = FuncAnimation(fig, animate,frames=10000,interval=200, blit=False)
    plt.show()

#-------------------------------------------
# CONVENIENCE FUNCTIONS
#-------------------------------------------

def mean_queue_length(road):
    return np.mean([sum(g) for b, g in itertools.groupby(road) if b])


def queue_locations(road,length_filter=0):
    queues = []
    in_queue =  False
    start    = 0
    length   = 0
    for i,elem in enumerate(road):
        if elem == 1 and not in_queue:
            in_queue = True
            start    = i 
            length  += 1
        elif elem == 1 and in_queue:
            length += 1
        elif elem == 0 and in_queue:
            if length >= length_filter:
                queues.append((start,i-1,length))
            in_queue = False
            length = 0
    return queues



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

#-------------------------------------------

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
    #np.random.seed(0)
    if args.random:
        r = lattice_random(args.l,args.n)
    elif args.d != None:
        r = lattice_density_10(args.l,args.d)
        print(r)
    else:
        r = lattice_even(args.l,args.n)
    
    if args.p:
        run_visual_simulation(r,args)
    else:
        run_simulation(r,args)
    