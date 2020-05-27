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
    """
    ---------------------
    Base class for car
    properties
    ---------------------
    :Notes:
    - tracking is quite computationally
      expensive in this implementation
    - set_speed function unnecessary:
      refactor out
    """

    def __init__(self,pos,index, lane, track=False):
        """ Initialization of class
        :input: pos   = position
                index = index number specific for car
                lane  = what lane car is at
                track = allow tracking options
        """
        self.position = pos
        self.index    = index
        self.lane     = lane
        self.speed    = 0
        self.moved    = 0
        self.prefer   = 0
        if track:
            self.trace    = []
            self.time     = []
            self.x        = []
            self.lane_history = []

    def set_speed(self,val):
        #unnecessary 
        self.speed = val

    def __str__(self):
        return "Car: " + str(self.index) + ",position: " + str(self.position) + ",speed: " + str(self.speed) + ",lane: " + str(self.lane)

#------------------------------------------

class simulation:
    """
    ---------------------------------------
    Class for handling simulation options
    ---------------------------------------
    :Notes:
    - initialization of lanes 
    - constructing arrays from car data
    - Has 3 different simulation modes
        - Nagel_Schreckenberg method:
          K. Nagel and M. Schreckenberg, 
          A cellular automaton model for freeway traffic, J. Phys.
          I, 2 (1992), pp. 2221–2229
        - Sigle lane kinetic monte carlo
        - multilane kinetic monte carlo 
    """

    def __init__(self, array, track,num_lanes = 0):
        """ Initialization of simulation eviroment
        :input: array = single lane or array of lanes
                track = allow tracking
                num_lanes = number of lanes
        :Notes: 
        - num_lanes != 0 allows multilane initialization
        - initializes arrays for used during simulation
        - cumulative array not necessary for NS method 
          so creationg of this is unnecessary with this
          option
        """
        self.cars = []
        self.R          = 0
        n = 0
        if num_lanes == 0:
            self.length = len(array)
            for i,elem in enumerate(array):
                if elem == 1:
                    self.cars.append(car(i,n,0,track))
                    n += 1
        else: 
            self.length = len(array[0])
            for l,lane in enumerate(array):
                for i,elem in enumerate(lane):
                    if elem == 1:
                        self.cars.append(car(i,n,l,track))
                        n += 1

        self.num   = n
        self.cumulative = np.zeros(self.num)
        self.global_time = 0
        self.num_lanes = num_lanes

    def get_array(self):
        """
        Construct road array from cars in the
        case of single lane simulation
        """
        arr = np.zeros(self.length)
        for elem in self.cars:
            arr[elem.position] = 1
        return arr
    
    def get_lanes(self):
        """
        Construnc array of lanes from cars in the 
        case of multiple lanes
        """
        lanes = []
        for _ in range(self.num_lanes):
            lanes.append(np.zeros(self.length))
        for c in self.cars:
            lanes[c.lane][c.position] = 1
        return lanes

    def distance(self,a,b):
        """ 
        Calculates distance to next car while
        also handling periodic borders
        """
        return abs(self.cars[b].position - self.cars[a].position - self.length*round((self.cars[b].position - self.cars[a].position)/self.length))

    def arrhenius(self,speed,rel):
        """
        Defines rates, uses a modified arrhenius 
        equation, this can have different forms
        :input: speed = speed of car
                rel   = relaxation or characteristic
                        time
        """
        return (rel)*np.exp(-1/speed)

    def get_random(self):
        """
        Get random event from possible events.
        In this case get index of car that is to be 
        moved
        """
        rand = np.random.rand()*self.cumulative[-1]
        for i in range(self.num):
            if i == 0:
                if rand <= self.cumulative[0]:
                    return i
            else:
                if (self.cumulative[i-1] < rand <= self.cumulative[i]):
                    return i

    def calculate_speeds(self,max_speed):
        """
        Calculate speed of cars in single lane case
        Also implements speed limiting
        :input: max_speed = speed limit
        """
        for i in range(self.num):
            if i == self.num - 1:
                self.cars[i].set_speed(min(self.distance(i,0),max_speed))
            else:
                self.cars[i].set_speed(min(self.distance(i,i+1),max_speed))


    def multilane_speed(self,lanes,num_lanes,max_speed):
        """
        Calculate speed of cars in multilane case. 
        This is a quite innefficient approach 
        :input:
        """
        for i,c in enumerate(self.cars):
            if lanes[c.lane][(c.position + 1) % self.length] == 1:
                e_right = 0
                e_left  = 0
                init_pos = c.position
                if c.lane == 0:    
                    while e_right < self.length:
                        if lanes[c.lane + 1][(init_pos + 1) % self.length] == 1:
                            break
                        else:
                            e_right  += 1
                            init_pos += 1
                elif c.lane == num_lanes - 1:
                    while e_left < self.length:
                        if lanes[c.lane - 1][(init_pos + 1) % self.length] == 1:
                            break
                        else:
                            e_left    += 1
                            init_pos  += 1
                else: 
                    while e_left < self.length:
                        if lanes[c.lane - 1][(init_pos + 1) % self.length] == 1:
                            break
                        else:
                            e_left        += 1
                            init_pos      += 1
                    while e_right < self.length:
                        if lanes[c.lane + 1][(init_pos + 1) % self.length] == 1:
                            break
                        else:
                            e_right       += 1
                            init_pos      += 1
                if e_left >= e_right:
                    c.prefer = -1
                else:
                    c.prefer = 1
                c.set_speed(min(max(e_left,e_right) + 1,max_speed))
            else:
                empty_lattice = 0
                init_pos = c.position
                while True:
                    if lanes[c.lane][(init_pos +1) % self.length] == 1:
                        break
                    else:
                        empty_lattice += 1
                        init_pos      += 1
                c.prefer = 0
                c.set_speed(min(empty_lattice + 1,max_speed))



    def NS_speed(self,max_speed, p):
        """
        NS method speed definitions 
        4 step process
        """
        for i,c in enumerate(self.cars):
            c.set_speed(min(c.speed +1,max_speed))
            if i == self.num -1:
                c.set_speed(min(c.speed,self.distance(i,0) - 1))
            else:
                c.set_speed(min(c.speed,self.distance(i,i+1) - 1))
            if np.random.rand() <= p:
                c.set_speed(max(0,c.speed - 1))

    def calculate_cumulative(self,rel):
        """
        Create cumulative rates array
        """
        summa = 0
        for i in range(self.num):
            summa += self.arrhenius(self.cars[i].speed,rel)
            self.cumulative[i] = summa

    def move_car(self,a,fixed):
        """
        Move car
        :input: a     = index of car
                fixed = {
                    True = one unit max jumo
                    }
        """
        if self.cars[a].speed > 1:
            if fixed:
                self.cars[a].position = (self.cars[a].position + 1) % self.length
                self.cars[a].moved   += 4
            else:
                self.cars[a].position = (self.cars[a].position + (self.cars[a].speed-1)) % self.length
                self.cars[a].moved   += 4

    def move_car_multilane(self,a,fixed):
        """
        Move car in multilane case, handles lane changes as well
        :input: a     = index of car
                fixed = {
                    True = one unit max jumo
                    }
        """
        if self.cars[a].speed > 1:
            if fixed:
                self.cars[a].position = (self.cars[a].position + 1) % self.length
                self.cars[a].moved   += 4
            else:
                self.cars[a].position = (self.cars[a].position + (self.cars[a].speed-1)) % self.length
                self.cars[a].moved   += 4
            self.cars[a].lane += self.cars[a].prefer
            self.cars[a].prefer = 0

    def update_history(self,t):
        """
        If cars allow tracking this updated 
        tracking data
        """
        for c in self.cars:
            c.trace.append(c.position)
            c.time.append(t)
            c.x.append(c.moved)
            c.lane_history.append(c.lane)
    
    def update_time(self):
        """
        Updates global time with following function
        u random number in [0,1]
        \Delta t  = - \frac{\log u}{R_{tot}}
        """
        self.global_time += -(np.log(np.random.rand()))/(self.cumulative[-1])

    def print_cars(self):
        """
        Print car information to stdout
        """
        print("Road: lengt: " + str(self.length) + ",cars: " + str(self.num))
        for c in self.cars:
            print(c)

    def kmc(self,sl,rel,fixed=True):
        """
        Single iteration of kinetic monte carlo
        simulation (single lane)
        """
        self.calculate_speeds(sl)
        self.calculate_cumulative(rel) 
        self.move_car(self.get_random(),fixed)
        self.update_time()

    def kmc_multilane(self,sl,rel,fixed=True):
        """
        Single iteration of kinetic monte carlo
        simulation (multilane)
        """
        self.multilane_speed(self.get_lanes(),self.num_lanes,sl)
        self.calculate_cumulative(rel) 
        self.move_car_multilane(self.get_random(),fixed)
        self.update_time()


    def Nagel_Schreckenberg(self,sl,prob):
        """
        Single iteration of NS method
        :Note:
          - global time update is not verified
        """
        self.NS_speed(sl,prob)
        for a in range(self.num):
                self.cars[a].position = (self.cars[a].position + self.cars[a].speed) % self.length

        self.global_time += 1

    def line(self,x,a,b):
        """
        Basic line equation
        """
        return a*x + b

    def flow(self):
        """
        Used to quantify free flow rates numerically
        """
        from scipy.optimize import curve_fit
        vels = []
        for c in self.cars:
            popt,pcov = curve_fit(self.line,c.x,c.time)
            vels.append(1.0/popt[0])
        print('Mean velocity:',np.mean(vels),np.std(vels),'kmh:',np.mean(vels)*3.6)


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
    parser.add_argument('--release',action='store_false',default=True,
                        help="Limit jump rate to 1 unit per MCstep")
    parser.add_argument('--random',action='store_true',default=False,
                        help="random positions,default=False")
    parser.add_argument('-pr',type=int,default=100,
                        help="print/plot rate, default=100")
    parser.add_argument('-d',type=int,default=None,
                        help="create initial array with density between [1]/10-[9]/10, default=None")
    parser.add_argument('-m','--max',type=int,default=1000,
                        help="maximum time steps,default=100")
    parser.add_argument('-sl','--speedlimit',type=int,default=10,
                        help="impose speed limit on cars, Default=None")
    parser.add_argument('--type',choices=["kmc","ns"],default="kmc",
                        help="select simulation type"
                             "kmc = Basic kinetic monte carlo [default]"
                             "ns  = Nagel Schreckenberg method")
    parser.add_argument('--units',choices=["meters","cells"],default="meters",
                        help="Units in plot"
                             "meters  [default]"
                             "cells")
    parser.add_argument('--case',choices=["1","2","3"],default=None)
    parser.add_argument('-p',type=float,default=0.5,
                        help="braking probability in NS method, default= 0.5")
    parser.add_argument('-w',type=float,default=6.0,
                        help="Characteristic time w_0 in arrhenius equation,"
                             "Should correspond to ≈ 80 km/h in sparse case default=6.0")
    parser.add_argument('--seed',type=int,default=None,
                        help="give seed to random number generator")
    return parser.parse_args()
    
#-------------------------------------------
# SIMULATION MAIN METHODS
#-------------------------------------------

def run_simulation(r,args):
    """ 
    Run single non visual simulation no plotting
    """
    index = np.array([i for i in range(r.size)])
    road = simulation(r,track=True)
    mqls        = []
    for ii in range(args.max):
        for _ in range(args.pr):
            road.kmc(args.speedlimit,args.w)
        road.update_history(road.global_time)
        mqls.append(mean_queue_length(road.get_array()))
    print("Mean queue length:", args.d,np.mean(mqls[len(mqls)//2:]),np.std(mqls[len(mqls)//2:]),len(road.get_array()),np.sum(road.get_array()))


def run_multilane(r,args):
    """
    Run simulation in multilane mode that draws plot at the end
    if chosen
    """
    roads = simulation(r,track=True,num_lanes=len(r))
    mqls  = []
    #roads.print_cars()
    #print(roads.get_lanes())
    for ii in range(args.max):
        for _ in range(args.pr):
            roads.kmc_multilane(args.speedlimit,args.w)
        roads.update_history(roads.global_time)
        temp = []
        for i in range(roads.num_lanes):
            temp.append(mean_queue_length(roads.get_lanes()[i]))
        mqls.append(np.mean(temp))

    print("Mean queue length:", args.d,np.mean(mqls[len(mqls)//2:]),np.std(mqls[len(mqls)//2:]),roads.num_lanes,roads.length,roads.num)
    import matplotlib.pyplot as plt
    if args.case == "2":
        fig = plt.figure()
        ax  = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)


        #roads.print_cars()
        for j in range(roads.num):
            m = np.array(roads.cars[j].lane_history,dtype=int) == 0
            if args.units == "meters":
                ax.plot([4*a for a in np.array(roads.cars[j].trace)[m]],np.array(roads.cars[j].time)[m],color="k",marker="s",ms=2,linestyle="")
            else:
                ax.plot(np.array(roads.cars[j].trace)[m],np.array(roads.cars[j].time)[m],color="k",marker="s",ms=2,linestyle="")
            m = np.array(roads.cars[j].lane_history,dtype=int) == 1
            if args.units == "meters":
                ax2.plot([4*a for a in np.array(roads.cars[j].trace)[m]],np.array(roads.cars[j].time)[m],color="k",marker="s",ms=2,linestyle="")
            else:
                ax2.plot(np.array(roads.cars[j].trace)[m],np.array(roads.cars[j].time)[m],color="k",marker="s",ms=2,linestyle="")
        if args.units == "meters":
            ax2.set_xlim(0, 4*args.l)
            ax.set_xlim(0, 4*args.l)
            ax2.set_xlabel("Position on road [m]",fontsize=20)
            ax.set_xlabel("Position on road [m]",fontsize=20)
        elif args.units == "cells":
            ax2.set_xlim(0, args.l)
            ax.set_xlim(0, args.l)
            ax2.set_xlabel("Position on road [m]",fontsize=20)
            ax.set_xlabel("Position on road [m]",fontsize=20)

        ax.set_ylabel("time [s]",fontsize=20)
        ax2.set_ylabel("time [s]",fontsize=20)   
        ax.set_title("Lane 0",fontsize=20)
        ax2.set_title("Lane 1",fontsize=20) 
        ax2.set_ylim(0,100)
        ax.set_ylim(0,100)
        plt.show()

        


#-------------------------------------------

def run_visual_simulation(r,args):
    """
    Run simulation in visual mode that updates the plot 
    in real time
    """
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
    if args.units == "meters":
        ax2.set_xlim(0, 4*args.l)
    elif args.units == "cells":
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
        road.flow()
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
                road.kmc(args.speedlimit,args.w,args.release)
        elif args.type == "ns":
            road.Nagel_Schreckenberg(args.speedlimit,args.p)

        road.update_history(road.global_time)
        r = road.get_array() 
        mask  = r == 1
        circle.set_data(index[mask]*(2*np.pi/len(index)),r[mask])
        title.set_text("FRAME " + str(i)) 

        for tn,trace in enumerate(traces):
            #print(tn,cars[tn][0],cars[tn][1])
            if args.units=="cells":
                trace.set_data((road.cars[tn].trace),road.cars[tn].time)
            elif args.units=="meters":
                trace.set_data([4*a for a in (road.cars[tn].trace)],road.cars[tn].time)
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
    """
    Helper function 
    Used in benchmarking
    """
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

def lattice_density(L,denominator,numerator = 10):
    stencil = np.zeros(numerator)
    stencil[:denominator] = 1
    road = np.tile(stencil,round(L/numerator))
    remainder = L % numerator
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
        r = lattice_density(args.l,args.d)
    else:
        r = lattice_even(args.l)
    
    if args.noplot:
        if args.case == "1":
            for i in range(1,10):
                args.d = i
                r = lattice_density(args.l,args.d)
                for _ in range(1):
                    run_simulation(r,args)
        elif args.case == "2":
            r = [lattice_even(args.l),np.zeros(args.l)]
            args.n = int(args.l/2)
            run_multilane(r,args)
        elif args.case == "3":
            for i in range(2,6):
                for d in range(1,10):
                    args.d = d
                    r = lattice_density(args.l,args.d)
                    roads = []
                    for _ in range(i):
                        roads.append(r)
                    run_multilane(roads,args)
        else:
            run_simulation(r,args)
    else:
        run_visual_simulation(r,args)
    