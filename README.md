# Monte Carlo trafficjam simulator

Simple monte carlo traffic jam simulator


## Instruction and using the code

There are program supports 3 different modes:
1. Single lane 
2. Multilane 
3. Visual

supports in single lane in both kinetic monte carlo


## Options 
usage: trafficjam.py [-h] [-l L] [-n N] [--noplot] [--release] [--random]
                     [-pr PR] [-d D] [-m MAX] [-sl SPEEDLIMIT]
                     [--type {kmc,ns}] [--units {meters,cells}]
                     [--case {1,2,3}] [-p P] [-w W] [--seed SEED]

optional arguments:
  -h, --help            show this help message and exit
  -l L                  length of road, default=300
  -n N                  number of cars, default=150
  --noplot              plot situation,default=False
  --release             Limit jump rate to 1 unit per MCstep
  --random              random positions,default=False
  -pr PR                print/plot rate, default=100
  -d D                  create initial array with density between
                        [1]/10-[9]/10, default=None
  -m MAX, --max MAX     maximum time steps,default=100
  -sl SPEEDLIMIT, --speedlimit SPEEDLIMIT
                        impose speed limit on cars, Default=None
  --type {kmc,ns}       select simulation typekmc = Basic kinetic monte carlo
                        [default]ns = Nagel Schreckenberg method
  --units {meters,cells}
                        Units in plotmeters [default]cells
  --case {1,2,3}
  -p P                  braking probability in NS method, default= 0.5
  -w W                  Characteristic time w_0 in arrhenius equation,Should
                        correspond to ≈ 80 km/h in sparse case default=6.0
  --seed SEED           give seed to random number generator
