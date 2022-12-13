#### REASLISTIC AGENT-BASED NPT VIRUS MODEL ####

# This agent-based model tracks aphids and plants in a Donnelly-model style
# of plant infection, i.e. aphids lose infectivity after one probe, and keep
# flying and probing until they settle down to feed.

# %%

from random import uniform, choices, randint, seed
import bisect
from re import T
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
        self.time_of_next_event = TIMEFRAME * 2 # plant is S again so set very high (as S plants don't recover)

        

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

def run_simulation(diagnostics=True, track_Z=True):
    time = 0

    plants, aphids, init_states = initialise_simulation()
    
    agents = plants + aphids
    agents.sort(key=lambda agent: agent.time_of_next_event) # sort all plants + aphids by time of next event

    states = pd.DataFrame(data=init_states) # will add to states over the timeframe
    states = states.assign(time=time) # add in time column
    
    num_I_plants = init_states["I"][0]
    
    if diagnostics:
        transmissions_df = pd.DataFrame(columns = ["num_I_plants_before", "num_transmissions"])
    
        num_dispersals = 0
        dispersals_df = pd.DataFrame(data = [[0,0]], columns = ["time", "num_dispersals"])
    
    if track_Z:
        aphid_status_df = pd.DataFrame(columns = ["time", "num_infective_flights", "num_healthy_flights"])
        
    while time < TIMEFRAME:
        
        agent = agents.pop(0) # pop selects element 0 from agents and also removes it
        
        # set time as smallest time of next event
        time = agent.time_of_next_event
        
        # if is an aphid, initiate dispersal
        if type(agent) == Aphid:
            agent.feeding = False
            
            if diagnostics:
                num_I_before = num_I_plants
                num_dispersals += 1 # increment dispersals count
            
            if track_Z:
                num_infective_flights = 0
                num_healthy_flights = 0
            
            while agent.feeding == False:
                
                if track_Z:
                    if agent.infective:
                        num_infective_flights += 1
                    else:
                        num_healthy_flights += 1
                
    
                agent.fly(NUM_PLANTS, num_I_plants)
                num_I_plants = agent.probe(a, epsilon, omega, b, agents, time, num_I_plants)
                agent.lose_infectivity()
                
            # reset time to next dispersal event
            agent.time_of_next_event += np.random.exponential(eta)
                        
            if diagnostics:
                # add num_transmissions to df
                num_transmissions = num_I_plants - num_I_before
                transmissions_df = transmissions_df.append({"num_I_plants_before": num_I_before, 
                                                            "num_transmissions": num_transmissions},
                                                        ignore_index=True)
                # add num_dispersals to df
                dispersals_df = dispersals_df.append({"time": time,
                                                      "num_dispersals": num_dispersals},
                                                    ignore_index=True)

            if track_Z:
                # add aphid infectivity status data to df
                aphid_status_df = aphid_status_df.append({"time": time,
                                                        "num_infective_flights": num_infective_flights,
                                                        "num_healthy_flights": num_healthy_flights},
                                                        ignore_index=True)

        else: # type == Plant
            agent.recovery() # plant recovers
            num_I_plants -= 1
            


        # put agent back into sorted list based on new time_of_next_event
        bisect.insort_left(agents, agent)
        
        # add states at this timepoint to states df
        states = states.append({"I": num_I_plants,
                                "Xs": len([aphid for aphid in aphids if not aphid.infective and not aphid.on_I_plant]),
                                "Xi": len([aphid for aphid in aphids if not aphid.infective and aphid.on_I_plant]),
                                #"Zs": len([aphid for aphid in aphids if aphid.infective and not aphid.on_I_plant]),
                                #"Zi": len([aphid for aphid in aphids if aphid.infective and aphid.on_I_plant]),
                                "time": time
                                },
                               ignore_index=True)
    
    if diagnostics:
        transmissions_df = transmissions_df.astype(int) # change data type of cols from object to integer
        dispersals_df = dispersals_df.astype(float) # change data type to float

        if track_Z:
            return(states, transmissions_df, dispersals_df, aphid_status_df)
        return(states, transmissions_df, dispersals_df)
    
    elif track_Z:
        return(states, aphid_status_df)
    
    return states

# %%

if __name__ == "__main__":

    np.random.seed(10)
    seed(10)
    
    MAX_X = 10
    MAX_Y = 7
    NUM_PLANTS = MAX_Y * MAX_X
    NUM_INIT_INFECTED = NUM_PLANTS # plant(s) (NOTE: CANNOT BE > NUM_PLANTS)
    NUM_APHIDS = NUM_PLANTS * 2
    TIMEFRAME = 7
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

    NUM_RUNS = 200
    TRACK_Z = False # track number of infective aphids during simulation?
    DIAGNOSTICS = False # track number of dispersals and transmissions during simulation?
    TRACK_EXTINCTIONS = False # track whether or not runs go extinct and plot runs without extinctions in additional panel for figs 3,4,5
    
    result_all_runs = pd.DataFrame()
    
    if DIAGNOSTICS:
        transmissions_all_runs = pd.DataFrame()
        dispersals_all_runs = pd.DataFrame()
    if TRACK_Z:
        aphid_status_all_runs = pd.DataFrame()
        aphid_status_no_extinct = pd.DataFrame()
    if TRACK_EXTINCTIONS:
        result_no_extinct_runs = pd.DataFrame()


    # start_t = time.time()
    # res, aphid = run_simulation(diagnostics=False, track_Z=True)
    # run_time = time.time() - start_t
    # print(f"run time = {run_time}")

    num_extinct = 0
    extinct_runs = {}
    
    for run in range(NUM_RUNS):
        
        print(f"run: {run + 1}")
        
        #start_t = time.time()
        if DIAGNOSTICS:
            if TRACK_Z:
                result, transmissions_df, dispersals_df, aphid_status_df = run_simulation(diagnostics=DIAGNOSTICS,
                                                                                          track_Z=TRACK_Z)
            else:
                result, transmissions_df, dispersals_df = run_simulation(diagnostics=DIAGNOSTICS,
                                                                         track_Z=TRACK_Z)
        elif TRACK_Z:
            result, aphid_status_df = run_simulation(diagnostics=DIAGNOSTICS,
                                                     track_Z=TRACK_Z)
        else:
            result = run_simulation(diagnostics=DIAGNOSTICS,
                                    track_Z=TRACK_Z)
        #run_time = time.time() - start_t
        #print(f"run time = {run_time}")
        
   
        
        # set timesteps for interpolated result (states over time) and dispersals over time
        time_steps = np.linspace(start=0, stop=TIMEFRAME, num=TIMEFRAME*100)
        
        # interpolate main result (number in each state over time)
        interpolated_result = pd.DataFrame(data=time_steps, columns=["time"])
        for state in result.drop("time", axis=1).columns: # for each state, i.e. Xs, Xi, ((Zs, Zi))

            # get linear interpolation of current state with regular timesteps
            interp = interpolate.interp1d(result["time"], result[state], kind="linear")
            state_linear = interp(time_steps)

            interpolated_result.insert(loc=interpolated_result.shape[1], column=state, value=state_linear) # add to df

        # add result to df of all runs, and to df of non-extinct runs if final I != 0
        result_all_runs = result_all_runs.append(interpolated_result)
        
        if TRACK_EXTINCTIONS:
            if interpolated_result["I"].loc[interpolated_result.shape[0]-1]: # if final I != 0
                result_no_extinct_runs = result_no_extinct_runs.append(interpolated_result)
                
                if TRACK_Z:
                    aphid_status_no_extinct = aphid_status_no_extinct.append(aphid_status_df) 
            else:
                num_extinct += 1
                print("EXTINCT")
                extinct_runs["run" + str(run)] = interpolated_result
                
        if TRACK_Z:
            aphid_status_all_runs = aphid_status_all_runs.append(aphid_status_df)

        if DIAGNOSTICS:
            # interpolate aphid dispersals over time
            interp_disp = interpolate.interp1d(dispersals_df["time"], dispersals_df["num_dispersals"], kind="linear")
            interp_disp_df = pd.DataFrame(data = {"time": time_steps, "num_dispersals": interp_disp(time_steps)})
        
            # add all data to respective dfs of all runs
            transmissions_all_runs = transmissions_all_runs.append(transmissions_df)
            dispersals_all_runs = dispersals_all_runs.append(interp_disp_df)
        
    print(f"interp={np.mean(interpolated_result['I'].loc[interpolated_result['time']>2])}")
    print(f"all_runs={np.mean(result_all_runs['I'].loc[result_all_runs['time']>2])}")
    print(f"raw={np.mean(result['I'].loc[result['time']>2])}")
    # p.plot(result['time'], result['I'], label="raw")
    # p.plot(interpolated_result['time'], interpolated_result['I'], label="interpolated")
    # p.legend()
    # p.show()

    #sys.exit()

    # FIGURE 1 + 2: DIAGNOSTICS (NUM DISPERSALS AND TRANSMISSIONS PER DISPERSAL)
    if DIAGNOSTICS:
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
        p.xlabel("Time (arbitrary units)")
        p.ylabel("Total number of aphid dispersals")
        p.legend()
# %%
    # FIGURE 3: NUMBER IN EACH STATE OVER TIME
    if TRACK_EXTINCTIONS:
        
        fig3, (ax1, ax2) = p.subplots(nrows=1, ncols=2,
                              sharex=True, sharey=True,
                              tight_layout=True)
        ax = ax1 # so default (with extinctions) plot goes on first axes
    else:
        fig3, ax = p.subplots(nrows=1, ncols=1,
                              sharex=True, sharey=True,
                              tight_layout=True)
    
    # SUBPLOT 1: ALL RUNS
    # find median of each interpolated time point
    median_interpolated_res = result_all_runs.pivot_table(index="time", aggfunc="mean")
    ax.plot(median_interpolated_res)#.plot()
    ax.set_ylabel("Number")
    ax.set_xlabel("Time (arbitrary units)")
    ax.legend(median_interpolated_res,
                             loc = "center right", framealpha = 0.5)
    # first_legend = ax.legend(median_interpolated_res,
    #                         loc = "center right", framealpha = 0.5)
    # p.gca().add_artist(first_legend) # allows 2 legends on same plot

    # horizontal line showing total num plants
    ax.axhline(y = NUM_PLANTS, color = "black", linestyle = "--", label = "Total number of plants") 
    ax.text(x = 0.2, y = NUM_PLANTS, s = "number of plants", fontdict = {"size": 9})
    
    # horizonal line showing expected final I according to Donnelly ode
    #ax1.axhline(y = EXPECTED_FINAL_I, color = "red", linestyle = "--", label = "Expected final I")
    #ax1.text(x = 0.2, y = EXPECTED_FINAL_I, s = "EXPECTED FINAL I", fontdict = {"size": 9})
    #ax1.legend(loc = "upper left")

    plot_text = "\n".join((
        r"num runs=%d" % NUM_RUNS,
        r"num plants=%d" % NUM_PLANTS,
        r"init infected=%d" % NUM_INIT_INFECTED,
        r"num aphids=%d" % NUM_APHIDS,
        r"eta=%.1f" % eta,       
        r"omega=%.1f" % omega,
        r"c_d=%.1f" % c_d
    ))
    ax.text(0.99, 0.99, plot_text,
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
        bbox=dict(facecolor="salmon", alpha=0.5, boxstyle="round"))
    
    if TRACK_EXTINCTIONS:
        # SUBPLOT 2: ALL RUNS WHERE EPIDEMIC PERSISTS
        median_no_extinct = result_no_extinct_runs.pivot_table(index="time", aggfunc="median")
        ax2.plot(median_no_extinct)
        ax2.set_xlabel("Time (arbitrary units)")
        #ax2.legend(median_no_extinct)
   
   
    # FIGURE 4: SHOW I TRAJECTORY COMPARED TO ODEs
    
    if TRACK_EXTINCTIONS:
        fig4, (ax1, ax2) = p.subplots(nrows=1, ncols=2,
                            sharex=True, sharey=True,
                            tight_layout=True)
        ax = ax1
    
    else:
        fig4, ax = p.subplots(nrows=1, ncols=1,
                            sharex=True, sharey=True,
                            tight_layout=True)
    
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
    
    ax.plot(time_steps, median_interpolated_res["I"], label = "Agent-based simulation")
    ax.plot(time_steps, don_I, label = "Donnelly model")
    ax.plot(time_steps, new_ode_I_Z["I"], label = "New compartmental model")
    ax.set_xlabel("Time (arbitrary units)")
    ax.set_ylabel("Number of infected plants")
    ax.legend()
    
    if TRACK_EXTINCTIONS:
        ax2.plot(time_steps, median_no_extinct["I"])
        ax2.plot(time_steps, don_I)
        ax2.plot(time_steps, new_ode_I_Z["I"])
        ax2.set_xlabel("Time (arbitrary units)")
    
    # %%
    
    # FIGURE 5: SHOW Z TRAJECTORY COMPARED TO ODEs
    if TRACK_Z:
       
        Z_aphids_dfs = []
        
        if TRACK_EXTINCTIONS:
            aphid_status_dfs = [aphid_status_all_runs, aphid_status_no_extinct]
            
            fig5, (ax1, ax2) = p.subplots(nrows=1, ncols=2,
                                      sharex=True, sharey=True,
                                      tight_layout=True)
        else:
            aphid_status_dfs = [aphid_status_all_runs]
            
            fig5, ax = p.subplots(nrows=1, ncols=1,
                                  sharex=True, sharey=True,
                                  tight_layout=True)
        
        for df in aphid_status_dfs:
            
            # create new column: proportion of flights that were infective
            df["prop_infective_flights"] = df.apply(lambda row: row.num_infective_flights / (row.num_infective_flights + row.num_healthy_flights), axis=1)
           # aphid_status_no_extinct["prop_infective_flights"] = aphid_status_no_extinct.apply(lambda row: row.num_infective_flights / (row.num_infective_flights + row.num_healthy_flights), axis=1)
            
            # create bins of times to average over
            df["time_bins"] = pd.cut(df["time"], bins = TIMEFRAME*30)
            #aphid_status_no_extinct["time_bins"] = pd.cut(aphid_status_no_extinct["time"], bins = TIMEFRAME*30)
            
            # find mean proportion of infective flights for each bin and * by num aphids to get Z (num infective aphids) for each bin
            Z_aphids_df = df.pivot_table(index="time_bins", values="prop_infective_flights", aggfunc=lambda x: np.mean(x) * NUM_APHIDS).reset_index().rename(columns={"prop_infective_flights": "Z"})
            #Z_aphids_no_extinct = aphid_status_no_extinct.pivot_table(index="time_bins", values="prop_infective_flights", aggfunc=lambda x: np.mean(x) * NUM_APHIDS).reset_index().rename(columns={"prop_infective_flights": "Z"})
            
            # find mid point of each bin to use as time x axis on plot
            Z_aphids_df["time"] = Z_aphids_df["time_bins"].transform(lambda t: t.mid)
            #Z_aphids_no_extinct["time"] = Z_aphids_no_extinct["time_bins"].transform(lambda t: t.mid)
            
            Z_aphids_dfs.append(Z_aphids_df)
        
        if TRACK_EXTINCTIONS:
            Z_aphids_df = Z_aphids_dfs[0]
            Z_aphids_no_extinct = Z_aphids_dfs[1]
            
            # plot Z for ABM
            ax1.plot(Z_aphids_df["time"], Z_aphids_df["Z"], label="Agent-based simulation")
            ax2.plot(Z_aphids_no_extinct["time"], Z_aphids_no_extinct["Z"])
            
            # add Z over time from Donnelly model to plot
            don_Z = (don_I/NUM_PLANTS)*a*(1 - omega)*NUM_APHIDS
            ax1.plot(time_steps, don_Z, label = "Donnelly model")
            ax2.plot(time_steps, don_Z)

            # add Z over time from new compartmental model to plot
            ax1.plot(time_steps, new_ode_I_Z["Z"], label="New compartmental model")
            ax2.plot(time_steps, new_ode_I_Z["Z"])
            
            ax1.set_xlabel("Time (arbitrary units)")
            ax1.set_ylabel("Number of infective aphids (Z)")
            ax1.legend()
            ax2.set_xlabel("Time (arbitrary units)")
            
        else:
            Z_aphids_df = Z_aphids_dfs[0]
            
            # plot Z for ABM
            p.plot(Z_aphids_df["time"], Z_aphids_df["Z"], label="Agent-based simulation")
            p.xlabel("Time (arbitrary units)")
            p.ylabel("Number of infective aphids (Z)")
            
            # add Z over time from Don model
            don_Z = (don_I/NUM_PLANTS)*a*(1 - omega)*NUM_APHIDS
            p.plot(time_steps, don_Z, label = "Donnelly model")
            
            # add Z over time from new compartmental model to plot
            p.plot(time_steps, new_ode_I_Z["Z"], label="New compartmental model")
            p.legend()
            
    if NUM_RUNS > 2:
        if DIAGNOSTICS:
            fig1.savefig(f"./simulation_mean_transmissions_{NUM_RUNS}_runs_mean.pdf")
            fig2.savefig(f"./simulation_mean_dispersals_{NUM_RUNS}_runs_mean.pdf")
            
        fig3.savefig(f"./simulation_res_{NUM_RUNS}_runs_mean.pdf")
        fig4.savefig(f"./simulation_ode_comparison_{NUM_RUNS}_runs_mean.pdf")
        
        if TRACK_Z:
            fig5.savefig(f"./simulation_Z_comparison_{NUM_RUNS}_runs_mean.pdf")
    
    p.show()

            
# %%
