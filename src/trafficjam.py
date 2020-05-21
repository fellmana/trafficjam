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

    def calculate_speeds(self,max_speed):
        for i in range(self.num):
            if i == self.num - 1:
                self.cars[i].set_speed(min(self.distance(i,0),max_speed))
            else:
                self.cars[i].set_speed(min(self.distance(i,i+1),max_speed))

    def NS_speed(self,max_speed, p):
        for i,c in enumerate(self.cars):
            c.set_speed(min(c.speed +1,max_speed))
            if i == self.num -1:
                c.set_speed(min(c.speed,self.distance(i,0) - 1))
            else:
                c.set_speed(min(c.speed,self.distance(i,i+1) - 1))
            if np.random.rand() <= p:
                c.set_speed(max(0,c.speed - 1))

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

    def kmc(self,sl):
        self.calculate_speeds(sl)
        self.calculate_cumulative() 
        self.move_car(self.get_random())
        self.update_time()


    def Nagel_Schreckenberg(self,sl,prob):
        self.NS_speed(sl,prob)
        for a in range(self.num):
                self.cars[a].position = (self.cars[a].position + self.cars[a].speed) % self.length

        self.global_time += 1


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
                        help="length of road, default=300")
    parser.add_argument('-n',type=int,default=None,
                        help="number of cars, default=150") 
    parser.add_argument('--noplot',action='store_true',default=False,
                        help="plot situation,default=False")
    parser.add_argument('--random',action='store_true',default=False,
                        help="random positions,default=False")
    parser.add_argument('-pr',type=int,default=100,
                        help="print/plot rate, default=100")
    
    parser.add_argument('-d',type=int,default=None,
                        help="create initial array with density between [1]/10-[9]/10, default=None")
    parser.add_argument('-m','--max',type=int,default=10**5,
                        help="maximum monte-carlo steps,default=10**6")
    parser.add_argument('-sl','--speedlimit',type=int,default=10,
                        help="impose speed limit on cars, Default=None")
    parser.add_argument('--type',choices=["kmc","ns"],default="kmc",
                        help="select simulation type"
                             "kmc = Basic kinetic monte carlo [default]"
                             "ns  = Nagel Schreckenberg method")
    parser.add_argument('-p',type=float,default=0.5,
                        help="braking probability in NS method, default= 0.5")
    parser.add_argument('--seed',type=int,default=None,
                        help="give seed to random number generator")
    return parser.parse_args()
    
#-------------------------------------------
# SIMULATION MAIN METHODS
#-------------------------------------------

def run_simulation(r,args):

    lane = simulation(r, track=False)
    lane.print_cars()
    road = lane.get_array()
    print(road)
    print(mean_queue_length(road))



    # lane.calculate_speeds()
    # lane.print_cars()
    # lane.calculate_cumulative()
    # print(lane.cumulative)
    # c = lane.get_random()
    # print(c)
    # print(lane.get_array())

    # lane.move_car(c)
    # lane.calculate_speeds()
    # lane.print_cars()
    # print(lane.get_array())
    
    # lane.calculate_cumulative()
    # print(lane.cumulative)
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
    for j in range(road.num):
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
        if args.type == "kmc": 
            for _ in range(args.pr):
                road.kmc(args.speedlimit)
        elif args.type == "ns":
            road.Nagel_Schreckenberg(args.speedlimit,args.p)

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
    """
    ----------------------------------------
    Calculates the mean queue length
    from given array
    ----------------------------------------
    :input:  road = array where 0 = empty
                   and 1 = occupied
    :output: mean queue length

    NOTES
    ----------------------------------------
    - Rolling is used to guarantee that a queue
      is not splitted over periodic border
    - With large dense arrays rolling might be 
      very inefficient
    - If you can accept reasonable accuracy 
      while loop could be turned of for better
      performance, SUGGESTION: if check for this?
    - Is is probably a smartet way to handel 
      periodic borders

    >>> mean_queue_length([1. 0. 0. 0. 1. 1. 0. 1. 0. 1.])
    1.6666666666666667
    """
    n = 0
    while True:
        if road[0] == 0 or n == len(road):
            break
        else:
            road = np.roll(road,1)
            n += 1
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

def lattice_even(L):
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

    if args.seed != None:
        np.random.seed(args.seed)
    if args.random or args.n != None:
        if args.n == None:
            args.n = int(args.l/2)
        r = lattice_random(args.l,args.n)
    elif args.d != None:
        r = lattice_density_10(args.l,args.d)
        print(r)
    else:
        r = lattice_even(args.l)
    
    if args.noplot:
        run_simulation(r,args)
    else:
        run_visual_simulation(r,args)
    