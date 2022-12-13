
"""
import numpy as np
import matplotlib.pyplot as p
import bisect
import pandas as pd

eta = 1
numAphids = 10
tMax=100

 allCounts = []
for zzzz in range(50):
    numDisp = 0
    nextJump = []
    for i in range(numAphids):
        nextJump.append(np.random.exponential(eta))
        
    t = 0
    tMax = 100
    while t < tMax:
        minEle = np.argmin(nextJump)
        t = nextJump[minEle]
        if nextJump[minEle] < tMax:
            numDisp = numDisp + 1
        nextJump[minEle] = nextJump[minEle] + np.random.exponential(eta)
    #print(f"{zzzz} after {tMax} have {numDisp} dispersals from {numAphids} aphids")
    
    allCounts.append(numDisp)
print(f"Average = {np.mean(allCounts)}")
 
vals = np.random.exponential(eta, size=1000)
p.hist(vals) # PLOT DISTRIBUTION OF TIME TO NEXT DISPERSAL WITH ETA = 1

 
nextEvent = list(np.random.exponential(eta, size = numAphids))
nextEvent.sort()

t = 0
numDispersal = 0
res1 = pd.DataFrame()
while t < tMax:
    min = nextEvent.pop(0)
    t += min
    nextEvent = [x - min for x in nextEvent]
    numDispersal += 1
    new = np.random.exponential(eta)
    bisect.insort(nextEvent, new)
    
    res1 = res1.append({'time': t, 'num_disp': numDispersal},
                       ignore_index=True)
    
t = 0
numDispersal = 0
res2 = pd.DataFrame()

nextEvent = list(np.random.exponential(eta, size = numAphids))
nextEvent.sort()

while t < tMax:
    min = nextEvent.pop(0)
    t += min
    numDispersal += 1
    bisect.insort(nextEvent, new)
    
    res2 = res2.append({'time': t, 'num_disp': numDispersal},
                       ignore_index=True)
    
p.plot(res1['time'], res1['num_disp'], label='new')
p.plot(res2['time'], res2['num_disp'], label='old')
x=np.linspace(0,tMax, num=50)
p.plot(x, x*numAphids*eta, label='expected')
p.legend()
    
    
 """
 
#### REASLISTIC AGENT-BASED NPT VIRUS MODEL ####

# This agent-based model tracks aphids and plants in a Donnelly-model style
# of plant infection, i.e. aphids lose infectivity after one probe, and keep
# flying and probing until they settle down to feed.

# %%

from random import uniform, choices, randint, seed
import bisect
import sys
import time
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.integrate import odeint
import matplotlib.pylab as p



class Agent:
    def __init__(self):
        self.time_of_next_event = float()

    # below functions mean that can compare 2 agents by their time to next event (e.g. using bisect)
    def __lt__(self, other): # less than
        return self.time_of_next_event < other.time_of_next_event
    
    def __le__(self, other): # less than or equal to
        return self.time_of_next_event <= other.time_of_next_event
    
    def __eq__(self, other): # equal to
        return self.time_of_next_event == other.time_of_next_event
        

class Aphid(Agent):
    def __init__(self, max_x, max_y, plants, time_feeding):
        self.x_pos = randint(0, max_x - 1) # -1 because coords are indexed from 0
        self.y_pos = randint(0, max_y - 1)
        self.infective = False
        self.on_I_plant = False
        self.flown_since_infective = False
        self.feeding = True # assume all aphids start feeding on a plant
        self.plants = plants
        self.time_of_next_event = time_feeding # next event = stop feeding + initiate dispersal

    def fly(self, num_plants, num_I_plants):

        if choices(["I", "S"], cum_weights = [v * num_I_plants / (num_plants - num_I_plants + v * num_I_plants), 1]) == ["I"]:
            
            selected_plants = [plant for plant in self.plants if plant.infected] # infected plants
            self.on_I_plant = True
        else:
            selected_plants = [plant for plant in self.plants if not plant.infected] # healthy plants
            self.on_I_plant = False

        new_plant = selected_plants[randint(0, len(selected_plants)-1)] # randomly select new plant and move to it
        self.x_pos = new_plant.x_pos
        self.y_pos = new_plant.y_pos


    def probe(self, a, epsilon, omega, b, agents, time, num_I_plants):
        if self.infective:
            self.flown_since_infective = True # already made a flight while infective

        if self.on_I_plant: # landed on I plant
            feed = uniform(0, 1)
            if uniform(0, 1) <= a and feed <= (1 - epsilon * omega): # (re-)aquire virus
                self.infective = True
                self.flown_since_infective = False # just been (re-)infected with virus
            elif feed > (1 - epsilon * omega): # i.e. <= epsilon*omega
                self.feeding = True

        else: # landed on S plant
            if self.infective: # if infective
                if uniform(0, 1) <= b: # inoculate plant
                    num_I_plants += 1
                    self.on_I_plant = True
                    current_plant = [plant for plant in self.plants if plant.x_pos == self.x_pos and plant.y_pos == self.y_pos][0]
                    current_plant.infected = True
                    current_plant.time_of_next_event = time + np.random.exponential(1/c_d) # set time to recovery
                    agents.pop(agents.index(current_plant)) # remove current_plant from the list of agents
                    bisect.insort_left(agents, current_plant) # re-insert current_plant into new correct place in list

                    
            if uniform(0, 1) <= omega:
                self.feeding = True
                
        return num_I_plants


    def lose_infectivity(self):
        if self.flown_since_infective:
            self.infective = False # lose infectivity
            self.flown_since_infective = False # re-set


class Plant(Agent):
    def __init__(self, x_pos, y_pos):
        self.x_pos = x_pos
        self.y_pos = y_pos
        self.infected = False
        self.time_of_next_event = TIMEFRAME * 2 # next event = recovery - plant by default is S, so large number for time so plant is not selected in Gillespie alg

    def recovery(self):
        self.infected = False
        

def binary_search(agents, agent, start, end):
    
    """ Function that could be used as an alternative to bisect.insort, as it finds the index that the agent should be put into in agents
    to maintain the sorted list, according to time_of_next_event. Would be used in conjunction with agents.insert(index returned from binary_search, agent)"""
    
    # in the case that the search gets smaller until start == end, if agent.time_of_next_event <= start, that's it's index,
    # otherwise it's start+1
    if start == end:
        if agent > agents[start]:
            return start + 1
        # else agent.time_of_next_event > start
        return start
        
    # in the case that the agent.time_of_next_event is smaller or larger than any in agents, the start point may become larger 
    # than the end (due to the mid-1 and mid+1 in the search respectively). in this case, return the start (index of 0 or len respectively)
    if start > end:
        return start
    
    # algorithm
    mid = (start + end) // 2
    if agent < agents[mid]: # if time_of_next event of agent is less than the midpoint
        binary_search(agents, agent, start, mid - 1) # do another binary search left of the midpoint
    elif agent > agents[mid]: # else if time_of_next_event of agent is greater than midpoint
        binary_search(agents, agent, mid + 1, end) # do another binary search right of midpoint
    else: # else agent.time_of_next_event == midpoint, so return it
        return mid

def initialise_simulation():
    # initialise plants
    plants = []
    count = 0
    for x in range(0, MAX_X):
        for y in range(0, MAX_Y):
            count += 1
            new_plant = Plant(x, y)
            if count <= NUM_INIT_INFECTED:
                new_plant.infected = True # by default is false
                new_plant.time_of_next_event = np.random.exponential(1/c_d) # this exponential func takes a scale parm = 1/rate parm
            plants.append(new_plant)

    # initialise aphids   
    aphids = [Aphid(MAX_X, MAX_Y, plants, np.random.exponential(eta)) for _ in range(NUM_APHIDS)] # this exponential func takes a scale parm = 1/rate parm

    for aphid in aphids: # correctly set on_I_plant attribute for all aphids
        current_plant = [plant for plant in plants if plant.x_pos == aphid.x_pos and plant.y_pos == aphid.y_pos][0]
        if current_plant.infected:
            aphid.on_I_plant = True

    # create initial states dict
    init_states = {
        "I": [NUM_INIT_INFECTED],
        "Xs": [len([aphid for aphid in aphids if aphid.on_I_plant == False])],
        "Xi": [len([aphid for aphid in aphids if aphid.on_I_plant == True])]#,
        #"Zs": [0],
        #"Zi": [0]
    }
    
    return plants, aphids, init_states

# %%

def run_simulation():
    time = 0

    plants, aphids, init_states = initialise_simulation()
    
    agents = plants + aphids
    agents.sort(key=lambda agent: agent.time_of_next_event) # sort all plants + aphids by time of next event

    states = pd.DataFrame(data=init_states) # will add to states over the timeframe
    states = states.assign(time=time) # add in time column
        
    transmissions_df = pd.DataFrame(columns = ["num_I_plants_before", "num_transmissions"])
    
    num_dispersals = 0
    dispersals_df = pd.DataFrame(data = [[0,0]], columns = ["time", "num_dispersals"])
    aphid_status_df = pd.DataFrame(columns = ["time", "num_infective_flights", "num_healthy_flights"])
        
    while time < TIMEFRAME:
        
        num_I_plants = len([plant for plant in plants if plant.infected])
        
        agent = agents.pop(0) # pop selects element 0 from agents and also removes it
        
        # set time as smallest time of next event
        time = agent.time_of_next_event
        
        # if is an aphid, initiate dispersal
        if type(agent) == Aphid:
            agent.feeding = False
            
            num_I_before = num_I_plants
            num_dispersals += 1 # increment dispersals count
            
            num_infective_flights = 0
            num_healthy_flights = 0
            
            while agent.feeding == False:
                if agent.infective == True:
                    num_infective_flights += 1
                else:
                    num_healthy_flights += 1
                
                agent.fly(NUM_PLANTS, num_I_plants)
                num_I_plants = agent.probe(a, epsilon, omega, b, agents, time, num_I_plants)
                agent.lose_infectivity()
                num_I_plants = len([plant for plant in plants if plant.infected])
                
            # add aphid infection status data to df
            aphid_status_df = aphid_status_df.append({"time": time,
                                                      "num_infective_flights": num_infective_flights,
                                                      "num_healthy_flights": num_healthy_flights},
                                                     ignore_index=True)
            
            num_transmissions = num_I_plants - num_I_before
            transmissions_df = transmissions_df.append({"num_I_plants_before": num_I_before, 
                                                        "num_transmissions": num_transmissions},
                                                       ignore_index=True)
            
            # reset time to next dispersal event
            agent.time_of_next_event += np.random.exponential(eta)
        
        else: # type == Plant
            agent.recovery() # plant recovers
            agent.time_of_next_event = TIMEFRAME * 2 # plant is S again so set very high (as S plants don't recover)

        # add num_dispersals to df
        dispersals_df = dispersals_df.append({"time": time,
                                              "num_dispersals": num_dispersals},
                                              ignore_index=True)

        # put agent back into sorted list based on new time_of_next_event
        bisect.insort_left(agents, agent)
        
        # add states at this timepoint to states df
        current_states = {
            "I": len([plant for plant in plants if plant.infected]),
            "Xs": len([aphid for aphid in aphids if not aphid.infective and not aphid.on_I_plant]),
            "Xi": len([aphid for aphid in aphids if not aphid.infective and aphid.on_I_plant]),
            #"Zs": len([aphid for aphid in aphids if aphid.infective and not aphid.on_I_plant]),
            #"Zi": len([aphid for aphid in aphids if aphid.infective and aphid.on_I_plant]),
            "time": time
        }

        states = states.append(current_states, ignore_index=True)
                
    transmissions_df = transmissions_df.astype(int) # change data type of cols from object to integer
    dispersals_df = dispersals_df.astype(float) # change data type to float
    
    return(states, transmissions_df, dispersals_df, aphid_status_df)

# %%

if __name__ == "__main__":

    np.random.seed(10)
    seed(10)

    MAX_X = 10
    MAX_Y = 7
    NUM_PLANTS = MAX_Y * MAX_X
    NUM_INIT_INFECTED = NUM_PLANTS # plant(s) (NOTE: CANNOT BE > NUM_PLANTS)
    NUM_APHIDS = NUM_PLANTS * 2
    TIMEFRAME = 20
    #EXPECTED_FINAL_I = 18.75
    
    if NUM_INIT_INFECTED > NUM_PLANTS:
        print("ERROR: num_init_infected > num_plants")
        sys.exit()

    eta = 1.3
    omega = 0.2
    a = 1
    b = 1
    epsilon = 1
    v = 1
    c_d = 3 # this is the plant recovery rate, i.e. c+d

    NUM_RUNS = 1
    result_all_runs = pd.DataFrame()
    transmissions_all_runs = pd.DataFrame()
    dispersals_all_runs = pd.DataFrame()
    aphid_status_all_runs = pd.DataFrame()

    start_t = time.time()
    #for run in range(NUM_RUNS):
        #print(f"run: {run + 1}")

    result, transmissions_df, dispersals_df, aphid_status_df = run_simulation()
    run_time = time.time() - start_t
    print(f"run time = {run_time}")
    sys.exit()
"""
        # set timesteps for interpolated result (states over time) and dispersals over time
        time_steps = np.linspace(start=0, stop=TIMEFRAME, num=TIMEFRAME*100)
        
        # interpolate main result (number in each state over time)
        interpolated_result = pd.DataFrame(data=time_steps, columns=["time"])
        for state in result.drop("time", axis=1).columns: # for each state, i.e. Xs, Xi, ((Zs, Zi))

            # get linear interpolation of current state with regular timesteps
            interp = interpolate.interp1d(result["time"], result[state], kind="linear")
            state_linear = interp(time_steps)

            interpolated_result.insert(loc=interpolated_result.shape[1], column=state, value=state_linear) # add to df

        # interpolate aphid dispersals over time
        interp_disp = interpolate.interp1d(dispersals_df["time"], dispersals_df["num_dispersals"], kind="linear")
        interp_disp_df = pd.DataFrame(data = {"time": time_steps, "num_dispersals": interp_disp(time_steps)})
        
        # add all data to respective sdfs of all runs
        result_all_runs = result_all_runs.append(interpolated_result)
        transmissions_all_runs = transmissions_all_runs.append(transmissions_df)
        dispersals_all_runs = dispersals_all_runs.append(interp_disp_df)
        aphid_status_all_runs = aphid_status_all_runs.append(aphid_status_df)


    
    # find mean number of transmissions per aphid dispersal at each stage of the epidemic (each number of I plants)
    mean_transmissions = transmissions_all_runs.pivot_table(index="num_I_plants_before", values="num_transmissions", aggfunc="mean").reset_index()
    mean_transmissions["expected_transmissions"] = mean_transmissions["num_I_plants_before"].transform(lambda I: ((1-omega)/omega) * I*(NUM_PLANTS - I)/NUM_PLANTS**2)
    mean_transmissions["num_samples"] = transmissions_all_runs["num_I_plants_before"].value_counts().sort_index().values    

    # FIGURE 1: MEAN TRANSMISSIONS PER DISPERSAL PER NUMBER OF INFECTED PLANTS
    fig1, ax = p.subplots()
    
    p.plot(mean_transmissions["num_I_plants_before"], mean_transmissions["expected_transmissions"], label = "Expected from Donnelly model")
    p.plot(mean_transmissions["num_I_plants_before"], mean_transmissions["num_transmissions"], label = "Agent-based model")
    p.xlabel("Number of infected plants")
    p.ylabel("Number of transmissions per aphid dispersal")
    p.legend()
    
    # find mean number of dispersals at each interpolated time point
    mean_dispersals = dispersals_all_runs.pivot_table(index="time", aggfunc="mean").reset_index()
    
    mean_dispersals["expected_dispersals"] = mean_dispersals["time"].transform(lambda t: (1/eta)*NUM_APHIDS*t) # 1/eta=theta
    
    # FIGURE 2: NUMBER OF APHID DISPERSALS OVER TIME COMPARED TO EXPECTED FROM DONNELLY
    fig2, ax = p.subplots()
    p.plot(mean_dispersals["time"], mean_dispersals["expected_dispersals"], label = "Expected from Donnelly model")
    p.plot(mean_dispersals["time"], mean_dispersals["num_dispersals"], label = "Agent-based model")
    p.xlabel("Time (arbitrary units")
    p.ylabel("Total number of aphid dispersals")
    p.legend()
        
    # FIGURE 3: NUMBER IN EACH STATE OVER TIME
    fig3, ax = p.subplots()
    
    # find median of each interpolated time point
    median_interpolated_res = result_all_runs.pivot_table(index="time", aggfunc="median")
    p.plot(median_interpolated_res)#.plot()
    p.ylabel("Number")
    p.xlabel("Time")
    first_legend = p.legend(median_interpolated_res,
                            loc = "center right", framealpha = 0.5)  
    p.gca().add_artist(first_legend) # allows 2 legends on same plot

    # horizontal line showing total num plants
    p.axhline(y = NUM_PLANTS, color = "black", linestyle = "--", label = "Total number of plants") 
    #p.text(x = 0.2, y = NUM_PLANTS, s = "number of plants", fontdict = {"size": 9})
    
    # horizonal line showing expected final I according to Donnelly ode
    #p.axhline(y = EXPECTED_FINAL_I, color = "red", linestyle = "--", label = "Expected final I")
    #p.text(x = 0.2, y = EXPECTED_FINAL_I, s = "EXPECTED FINAL I", fontdict = {"size": 9})
    #p.legend(loc = "upper left")

    plot_text = "\n".join((
        r"num runs=%d" % NUM_RUNS,
        r"num plants=%d" % NUM_PLANTS,
        r"init infected=%d" % NUM_INIT_INFECTED,
        r"num aphids=%d" % NUM_APHIDS,
        r"eta=%.1f" % eta,       
        r"omega=%.1f" % omega,
        r"c_d=%.1f" % c_d
    ))
    p.text(0.99, 0.99, plot_text,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
        bbox=dict(facecolor="salmon", alpha=0.5, boxstyle="round")
        )
   
    # FIGURE 4: SHOW I TRAJECTORY COMPARED TO ODEs
    fig4, ax = p.subplots()
    
    # first run donnelly ode
    def donnelly_ode(states, times, c_d, eta, H, v, epsilon, omega, a, b, A):
        I = states
        
        gamma = c_d
        theta = 1/eta
        w = omega
        e = epsilon
        Pacq = a
        Pinoc = b
        
        i_hat = v*I / (H - I + v*I)
        xI = (Pacq*Pinoc*i_hat*(1 - e*w)*(1 - i_hat))/(w*(1 - i_hat*(1 - e)))
        dI = (theta*A)*xI - gamma*I
        
        return dI
    
    don_I = odeint(func=donnelly_ode, y0=NUM_INIT_INFECTED, t=time_steps, args=(c_d, eta, NUM_PLANTS, v, epsilon, omega, a, b, NUM_APHIDS))

    # now run new Donnelly-logic compartmental ode
    def new_compartment_ode(states, times, c_d, eta, H, v, epsilon, omega, a, b, A):
        I = states[0]
        Z = states[1]
        
        S = H - I
        X = A - Z
  
        phi = (S + v*I) / (omega*eta*(S + v*epsilon*I))

        dI = phi*Z*b*S/(S + v*I) - (c_d)*I
        dZ = phi*X*a*(1 - epsilon*omega)*v*I/(S + v*I) - phi*Z*(1 - a*(1 - epsilon*omega))*v*I/(S + v*I) - phi*Z*S/(S + v*I)

        return dI, dZ

    INIT_Z = NUM_APHIDS*a*(1 - epsilon*omega)*v*NUM_INIT_INFECTED / (v*NUM_INIT_INFECTED + NUM_PLANTS - NUM_INIT_INFECTED)
    new_ode_I_Z = pd.DataFrame(odeint(func=new_compartment_ode, y0=[NUM_INIT_INFECTED, INIT_Z], t=time_steps, args=(c_d, eta, NUM_PLANTS, v, epsilon, omega, a, b, NUM_APHIDS)))
    new_ode_I_Z.columns = ["I", "Z"]
    
    p.plot(time_steps, median_interpolated_res["I"], label = "Agent-based simulation")
    p.plot(time_steps, don_I, label = "Donnelly model")
    p.plot(time_steps, new_ode_I_Z["I"], label = "New compartmental model")
    p.xlabel("Time (arbitrary units)")
    p.ylabel("Number of infected plants")
    p.legend()
    
    # %%
    # FIGURE 5: SHOW Z TRAJECTORY COMPARED TO ODEs
    fig5, ax = p.subplots()
    
    # create new column: proportion of flights that were infective
    aphid_status_all_runs["prop_infective_flights"] = aphid_status_all_runs.apply(lambda row: row.num_infective_flights / (row.num_infective_flights + row.num_healthy_flights), axis=1)
    
    # create bins of times to average over
    aphid_status_all_runs["time_bins"] = pd.cut(aphid_status_all_runs["time"], bins = TIMEFRAME*20)
    
    # find mean proportion of infective flights for each bin and * by num aphids to get Z (num infective aphids) for each bin
    Z_aphids_df = aphid_status_all_runs.pivot_table(index="time_bins", values="prop_infective_flights", aggfunc=lambda x: np.mean(x) * NUM_APHIDS).reset_index().rename(columns={"prop_infective_flights": "Z"})
    
    # find mid point of each bin to use as time x axis on plot
    Z_aphids_df["time"] = Z_aphids_df["time_bins"].transform(lambda t: t.mid)
    p.plot(Z_aphids_df["time"], Z_aphids_df["Z"], label="Agent-based simulation")
    
    # add Z over time from Donnelly model to plot
    don_Z = (don_I/NUM_PLANTS)*a*(1 - omega)*NUM_APHIDS
    p.plot(time_steps, don_Z, label = "Donnelly model")

    # add Z over time from new compartmental model to plot
    p.plot(time_steps, new_ode_I_Z["Z"], label="New compartmental model")
    
    p.xlabel("Time (arbitrary units)")
    p.ylabel("Number of infective aphids (Z)")
    p.legend()
    
    if NUM_RUNS > 2:
        fig1.savefig(f"./simulation_mean_transmissions_{NUM_RUNS}_runs.pdf")
        fig2.savefig(f"./simulation_mean_dispersals_{NUM_RUNS}_runs.pdf")
        fig3.savefig(f"./simulation_res_{NUM_RUNS}_runs.pdf")
        fig4.savefig(f"./simulation_ode_comparison_{NUM_RUNS}_runs.pdf")
        fig5.savefig(f"./simulation_Z_comparison_{NUM_RUNS}_runs.pdf")
    
    p.show()

            
# %%
"""