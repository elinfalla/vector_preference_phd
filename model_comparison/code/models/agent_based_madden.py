#### AGENT-BASED NPT VIRUS MODEL - based on Madden ####

# This an agent-based version of the 4 vector state Madden model, where tau = phi.
# It tracks a population of plants and aphids. Aphid feeding is dealt with implicitly 
# through setting phi as 1/w*eta.

# %%

from random import randint, uniform
from time import sleep
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %%
MAX_X = 3
MAX_Y = 3
NUM_PLANTS = MAX_Y * MAX_X
NUM_INIT_INFECTED = 1 # plant(s)
NUM_APHIDS = NUM_PLANTS * 3
TIMEFRAME = 20

eta = 1
omega = 0.2 # NOTE: 1/(eta*omega) needs to be an integer
a = 1
b = 1
epsilon = 1
v = 1
c_d = 4 # this is the plant recovery rate, i.e. c+d

class Aphid:
    def __init__(self, max_x, max_y, plants):
        self.x_pos = randint(0, max_x - 1) # -1 because coords are indexed from 0
        self.y_pos = randint(0, max_y - 1)
        self.infective = False
        self.on_I_plant = False
        self.flown_since_infective = False
        self.plants = plants

    def fly(self, num_plants, num_I_plants):

        if uniform(0,1) < (num_I_plants / num_plants): # would need changing if v!=1. e.g random.choices(["I", "S"], cum_weights = [vI/S+vI, 1])
            
            selected_plants = [plant for plant in self.plants if plant.infected] # infected plants
            self.on_I_plant = True
        else:
            selected_plants = [plant for plant in self.plants if not plant.infected] # healthy plants
            self.on_I_plant = False

        new_plant = selected_plants[randint(0, len(selected_plants)-1)] # randomly select new plant and move to it
        self.x_pos = new_plant.x_pos
        self.y_pos = new_plant.y_pos

    def probe(self, a, epsilon, omega, b):
        if self.infective == True:
            self.flown_since_infective = True # already made a flight while infective

        if self.on_I_plant == True: # landed on I plant
            if uniform(0, 1) <= a and uniform(0, 1) <= (1 - epsilon * omega): # (re-)aquire virus
                self.infective = True
                self.flown_since_infective = False # just been (re-)infected with virus

        else: # landed on S plant
            if self.infective == True: # if infective
                if uniform(0, 1) <= b: # inoculate plant
                    self.on_I_plant = True
                    current_plant = [plant for plant in self.plants if plant.x_pos == self.x_pos and plant.y_pos == self.y_pos]
                    #current_plant = get_current_plant(self, self.plants)
                    if len(current_plant) > 1:
                        print("ERROR!!")
                    current_plant[0].infected = True

    def lose_infectivity(self):
        if self.flown_since_infective == True:
            self.infective = False # lose infectivity
            self.flown_since_infective = False # re-set

    

class Plant:
    def __init__(self, x_pos, y_pos):
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.infected = False
        self.time_infected = 0

    def increment_time_infected(self, phi):
        if self.infected == True:
            self.time_infected += 1/phi 

    def recovery(self, c_d):
        if self.time_infected >= 1/c_d: # infectious period = 1/(c+d)
            self.infected = False
            self.time_infected = 0
# %%

""" def get_current_plant(Aphid, plants):
    aphid = Aphid
    current_plant = [plant for plant in plants if plant.x_pos == Aphid.x_pos and plant.y_pos == Aphid.y_pos]
    
    return(current_plant[0])   """  
# %%

def run_simulation(aphids, plants, init_states, parms):
        
    states = init_states # will append to states over timeframe
        
    for time in range(1, TIMEFRAME):

        num_I_plants = len([plant for plant in plants if plant.infected])
        #num_S_plants = NUM_PLANTS - num_I_plants

        for dispersal in range(0, int(parms["phi"])): # per timestep, do phi number of flights + probes
            
            # aphid actions
            for aphid in aphids:
                aphid.fly(NUM_PLANTS, num_I_plants)
                aphid.probe(parms["a"], parms["epsilon"], parms["omega"], parms["b"])
                aphid.lose_infectivity() # as tau = phi, lose_infectivity() happens at same rate as dispersal
                num_I_plants = len([plant for plant in plants if plant.infected]) # update num I plants
            
            # plant actions
            for plant in plants:
                plant.increment_time_infected(parms["phi"])
                plant.recovery(parms["c_d"]) # recover if been infected for 1/c_d amount of time

            num_I_plants = len([plant for plant in plants if plant.infected]) # update number infected plants

        states["I"].append(num_I_plants)
        states["Xs"].append(len([aphid for aphid in aphids if not aphid.infective and not aphid.on_I_plant]))
        states["Xi"].append(len([aphid for aphid in aphids if not aphid.infective and aphid.on_I_plant]))
        states["Zs"].append(len([aphid for aphid in aphids if aphid.infective and not aphid.on_I_plant]))
        states["Zi"].append(len([aphid for aphid in aphids if aphid.infective and aphid.on_I_plant]))
    
    return(states)

if __name__ == '__main__':
    
    phi = 1/(eta*omega)

    parms = {
        "phi": phi,
        "a": a, 
        "b": b,
        "omega": omega,
        "epsilon": epsilon, 
        "eta": eta, 
        "c_d": c_d
    }


    # %%
    times = [x for x in range(0, TIMEFRAME)]
    num_runs = 200
    all_runs_df = pd.DataFrame()

    for run in range(0, num_runs):
        print(f'Run: {run}')

        # initialise plants
        plants = []
        count = 0
        for x in range(0, MAX_X):
            for y in range(0, MAX_Y):
                count += 1
                new_plant = Plant(x, y)
                if count <= NUM_INIT_INFECTED:
                    new_plant.infected = True # by default is false
                plants.append(new_plant)
    
        # initialise aphids     # todo: make a certain number of aphids start on each plant
        aphids = [Aphid(MAX_X, MAX_Y, plants) for _ in range(0, NUM_APHIDS)]

        for aphid in aphids: # correctly set on_I_plant attribute for all aphids
            current_plant = [plant for plant in plants if plant.x_pos == aphid.x_pos and plant.y_pos == aphid.y_pos]
            if current_plant[0].infected == True:
                aphid.on_I_plant = True

        init_states = {
            "I": [NUM_INIT_INFECTED],
            "Xs": [len([aphid for aphid in aphids if aphid.on_I_plant == False])],
            "Xi": [len([aphid for aphid in aphids if aphid.on_I_plant == True])],
            "Zs": [0], 
            "Zi": [0]
        }   
        
        res = run_simulation(aphids, plants, init_states, parms)
        res["time"] = times

        res_df = pd.DataFrame(data=res)
        res_df_long = pd.melt(res_df, id_vars='time', value_vars=['I','Xs','Xi','Zs','Zi'], var_name='state', value_name='number')
        res_df_long = res_df_long.assign(run = [run + 1] * res_df_long.shape[0]) # new column for run number
        
        all_runs_df = all_runs_df.append(res_df_long) # append result to df with all runs

    all_runs_means = all_runs_df.groupby(['time', 'state']).agg(mean=('number', np.mean)).reset_index() #reset_index means groupedby cols stay as cols not indexes
    '''all_runs_means = all_runs_df.pivot_table(index=['time', 'state'],
                        values='number',
                        aggfunc=np.mean).reset_index()''' # this does the same thing, but mean col is called 'number'
    
    # plot result (change df from long to wide first)
    all_runs_means.pivot_table(index='time', columns='state', values='mean').plot()
    plt.show()
    

""" time_plt = times + [TIMEFRAME]
    plt.plot(time_plt, num_I)
    plt.axis([0,TIMEFRAME,0,NUM_PLANTS])
    plt.ylabel("Number infected plants")
    plt.xlabel("Time")
    plt.show() """

"""     print(num_Xi)
    Xs_plot = plt.plot(time_plt, num_Xs, label="Xs")
    Xi_plot = plt.plot(time_plt, num_Xi, "r", label="Xi")
    Zs_plot = plt.plot(time_plt, num_Zs, "g", label="Zs")
    Zi_plot = plt.plot(time_plt, num_Zi, "b", label="Zi")
    plt.axis([0, TIMEFRAME, 0, NUM_APHIDS])
    plt.legend(loc="upper right")
    plt.ylabel("Number of aphids")
    plt.xlabel("Time")
    plt.show()   """  


# %%
