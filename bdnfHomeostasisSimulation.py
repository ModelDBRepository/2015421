# This script uses Brian2 to produce network activity that mimics BDNF and 
# glutamate excitotoxicity by modulating connectivity and neuronal cell death, 
# respectively. Instructions for downloading Brian2 can be found at 
# https://brian2.readthedocs.io/en/stable/introduction/install.html 

# Briefly, this code outputs which neurons are spiking and at what time they 
# are spiking, as well as which neurons were connected to which and, for
# excitatory to excitatory connections, the synaptic weight. We later analyzed
# this network activity  in Matlab (bdnfHomeostasisAnalysisMain.m) to generate 
# the data plotted in Figure 9 in O'Neill et al., Time-dependent homeostatic
# mechanisms underlie BDNF action on neural circuitry. Comms Bio, 2023.

# This code was adapted from Masquelier and Deco, PLOS ONE, 2013 by Erin D.
# Anderson and can be accessed at https://www.seas.upenn.edu/~molneuro/

# Updated 11/14/2023
	   
#%% Import packages
from brian2 import *
from itertools import combinations
import networkx as nx
import math
import numpy as np
import numpy.matlib as matlib
from networkx import NetworkXError
import random
import scipy.io as sio
import gc
					 
def datetimeToFileLabelString():
    # Returns a YYYYMMDD_HHMMSS date/time string
    import datetime
    currentdatetime = datetime.datetime.now()
    currentdatetime = str(currentdatetime)
    currentdatetime = currentdatetime.replace("-","")
    currentdatetime = currentdatetime.replace(" ","_")
    currentdatetime = currentdatetime.replace(":","")
    currentdatetime = currentdatetime[0:15]
    return currentdatetime

#%% Defining recurrent network model parameters (almost all of them have been extracted from Masquelier et al. 2013 PlOS ONE)

# Simulation Parameters
simtime = 120*second          # Simulation time for control/injury steps
preinjury_simtime = 120*second # Simulation time for pre-injury step
defaultclock.dt = 0.1*ms      # simulation step length [ms]

condDelay = 3*ms              # experimental estimation of conditional delay of firing

# Pyramidal (excitatory) cells 
Cm_e        = 0.5*nF     # [nF] total capacitance 
gl_e        = 25.0*nS    # [ns] total leak conductance 
El_e        = -70.0*mV   # [mV] leak reversal potential 
Vth_e       = -50.0*mV   # [mV] threshold potential
Vr_e   = -60.0*mV        # [mV] reset potential 
tr_e   = 2.0*ms          # [ms] refractory time 

# Interneuron (inhibitory) cells 
Cm_i        = 0.2*nF     # [nF] total capacitance 
gl_i        = 20.0*nS    # [ns] total leak conductance 
El_i        = -70.0*mV   # [mV] leak reversal potential 
Vth_i       = -50.0*mV   # [mV] threshold potential 
Vr_i   = -60.0*mV        # [mV] reset potential 
tr_i   = 1.0*ms          # [ms] refractory time


# Synaptic weights [a.u.]
w_ee = 8
w_ei = w_ee
w_ie = 1
w_ii = 1

# Noise parameters
sigma_noise  = 6.5*mV        # gaussian white noise (Langevin equation)
					

# AMPA receptor
E_ampa   = 0.0*mV             # synaptic reversal potential (mV)
tau_ampa = 5.0*ms             # Glutamatergic synaptic time constant (AMPA)
g_ampa_e = 0.104*nS
g_ampa_i = 0.081*nS

# GABA receptor
E_gaba     = -70.0*mV        # inhibitory synaptic reversal potential (mV)
tau_gaba   = 10.0*ms         # GABAergic synaptic time constant
g_gaba_e   = 1.250*nS
g_gaba_i   = 0.973*nS

# NMDA receptor 
E_nmda     = 0.0*mV          # synaptic reversal potential (mV)
tau_s_nmda = 100.0*ms        # decay time of NMDA currents 
tau_x_nmda = 2.0*ms          # controls the rise time of NMDAR channels (ms)
a_nmda     = 0.5*kHz         # controls the saturation properties of NMDAR channels (kHz)
g_nmda_e = 0.327*nS
g_nmda_i = 0.258*nS

# STDP parameter for E -> E connections               
dApre = 1e-2
dApost = - 0.55*dApre
tau_post=34.0*ms             #(source: Bi & Poo 2001)
tau_pre=17.0*ms              #(source: Bi & Poo 2001) 
wmax = 2*w_ee

# STP parameters (Masquelier and Deco, 2013)
tau_f = 1600*ms              # Facilitation time scale
tau_d = 800*ms               # Depression time scale
U     = 0.025                

# Adaptation parameters (Masquelier and Deco, 2013); Modeling calcium-activated potassium current for spike-frequency adaptation (Wang et al 1999)
a_Ca = 0.0145
g_AHP = 10*nS
tau_Ca = 4000*ms
E_K = -80*mV

# External Input parameters
mu = 5
sigma = 1
I_on = 200*pA

# Synaptic connection density selected to match in vitro E/I synaptic balance
# Note: have to change manually in lines 178-180 due to a Brian2 quirk
connectionDensityEE = 0.2 
connectionDensityEI = 1
connectionDensityII = 1

EIBalance = 0.88 # selected to match in vitro E/I neuronal balance

#%% Defining neurodynamics equations

# Define equations describing the dynamics of excitatory neurons
exc_eqs = '''
            removed : boolean
            dv/dt=int(not removed)*((sigma_noise*xi*(2*gl_e/Cm_e)**.5) + (-gl_e*(v-El_e) - g_ampa_e*(v-E_ampa)*s_ampa_in - s_nmda_in*g_nmda_e*(v-E_nmda)/(1+exp(-0.062*v/volt)/3.57) - g_gaba_e*(v-E_gaba)*s_gaba - g_AHP*Ca*(v-E_K) + I)/Cm_e): volt (unless refractory)
            ds_ampa/dt = -s_ampa/tau_ampa : 1
            ds_nmda/dt = -(s_nmda/tau_s_nmda) + a_nmda*x*(1-s_nmda) : 1
            ds_gaba/dt = -s_gaba/(tau_gaba) : 1
            dx/dt = -x/tau_x_nmda : 1
            dCa/dt = -Ca/tau_Ca : 1
            I : amp
            s_ampa_in = s_ampa_in1 + s_ampa_in2 + s_ampa_in3: 1
            s_nmda_in = s_nmda_in1 + s_nmda_in2 + s_nmda_in3: 1
            s_ampa_in1 : 1
            s_ampa_in2 : 1
            s_ampa_in3 : 1
            s_nmda_in1 : 1
            s_nmda_in2 : 1
            s_nmda_in3 : 1
            tmp_nmda : 1
            tmp_ampa : 1
            xpos : 1
            ypos : 1
         '''

# Define equations describing the dynamics of inhibitory neurons
inh_eqs = '''
            removed : boolean
            dv/dt=int(not removed)*((sigma_noise*xi*(2*gl_i/Cm_i)**.5) + (-gl_i*(v-El_i) - g_ampa_i*(v-E_ampa)*s_ampa_in - s_nmda_in*g_nmda_i*(v-E_nmda)/(1+exp(-0.062*v/volt)/3.57) - g_gaba_i*(v-E_gaba)*s_gaba - g_AHP*Ca*(v-E_K) + I)/Cm_i): volt (unless refractory)
            ds_ampa/dt = -s_ampa/tau_ampa : 1
            ds_nmda/dt = -(s_nmda/tau_s_nmda) + a_nmda*x*(1-s_nmda) : 1
            ds_gaba/dt = -s_gaba/tau_gaba : 1
            dx/dt = -x/tau_x_nmda : 1
            dCa/dt = -Ca/tau_Ca : 1
            I : amp
            s_ampa_in = s_ampa_in1 + s_ampa_in2 + s_ampa_in3: 1
            s_nmda_in = s_nmda_in1 + s_nmda_in2 + s_nmda_in3: 1
            s_ampa_in1 : 1
            s_ampa_in2 : 1
            s_ampa_in3 : 1
            s_nmda_in1 : 1
            s_nmda_in2 : 1
            s_nmda_in3 : 1
            xpos : 1
            ypos : 1
         '''
    
# Define equations for distance dependence based on connection likelihood; source: Voges et al., 2012
# Note: have to manually define first multiplier for connection density & update in lines 122-124 to output correct values        
ConnectionProbability_ee = '0.2 * 0.8 * exp((-1*((xpos_post-xpos_pre)**2 + (ypos_post - ypos_pre)**2))/(2*(0.1**2)))'
ConnectionProbability_ei = '1 * 0.856 * exp((-1*((xpos_post-xpos_pre)**2 + (ypos_post - ypos_pre)**2))/(2*(0.0875**2)))'
ConnectionProbability_ii = '1 * 0.936 * exp((-1*((xpos_post-xpos_pre)**2 + (ypos_post - ypos_pre)**2))/(2*(0.075**2)))'

#%% Encode distance dependence

# Network Structure
CellDensity = 3500 # [cells/mm^2] cell density on MEA - 3500 cells/mm^2
RecordingWidth = 0.378  # [mm] width of square from which we "record" activity
                        # note: recording from a smaller area than full MEA to reduce the number of neurons we have to simulate to record from the same number of neurons as the MEA does
nReplicates = 6 # how many times to run each condition

for rep in range(nReplicates):
        print(f">>> Starting Sim {rep} of {nReplicates}")
        
        # Network structure
        N = int(round(RecordingWidth**2*CellDensity,0)) # Total # of neurons: recording area * cell density
        NE = math.floor(EIBalance*N)        # Number of Excitatory neurons
        NI = N - NE                         # Number of Inhibitory neurons
        
        # Define neuron XY locations:
        NeuronXPosition = np.random.rand(N,1) * RecordingWidth # [mm]
        NeuronYPosition = np.random.rand(N,1) * RecordingWidth # [mm]        
        
        #%% Define neuronal groups
        # Excitatory neuronal group
        Pe = NeuronGroup(NE, model=exc_eqs, threshold="v>Vth_e",
                        reset="v=Vr_e", refractory=tr_e, 
                        method='euler')
        Pe.xpos = list(NeuronXPosition[0:NE])
        Pe.ypos = list(NeuronYPosition[0:NE])
        
        # Inhibitory neuronal group
        Pi = NeuronGroup(NI, model=inh_eqs, threshold="v>Vth_i",
                        reset="v=Vr_i", refractory=tr_i, 
                        method='euler')
        Pi.xpos = list(NeuronXPosition[NE:])
        Pi.ypos = list(NeuronYPosition[NE:])
        
        #%% Define connections
        # Static Excitatory Connections (autapses)
        self_exc = Synapses(Pe, Pe, model = '''
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            s_nmda_in1_post = int(not removed) * tmp_nmda_pre : 1 (summed)
                                            s_ampa_in1_post = int(not removed) * tmp_ampa_pre : 1 (summed)
                                            w : 1            
                                            ''',
                             delay=condDelay,
                             on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_ampa_post += int(not removed) * w
                                    x_post += int(not removed) * w
                                    tmp_ampa = int(not removed) * clip(s_ampa*r_S*w_ee, 0, 100*r_S*w_ee)
                                    tmp_nmda = int(not removed) * clip(s_nmda*r_S*w_ee, 0, 100*r_S*w_ee)
                                    Ca += int(not removed) * a_Ca
                                    ''')
        
        self_exc.connect(condition='i==j', p=1) 
        self_exc.w = 1
        self_exc.u_f = U  #all synapses have fully-replenished neurotransmitter resources
        self_exc.x_d = 1  #all synapses have fully-replenished neurotransmitter resources
        
        
        # Static Inhibitory Connections (autapses)
        # [Autapses use short term plasticity (STP)]
        self_inh = Synapses(Pi, Pi, model = '''
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            w : 1
                                            ''', 
                             delay=condDelay,
                             on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_gaba_post += int(not removed) * r_S * w
                                    Ca += int(not removed) * a_Ca
                                    ''')
        
        self_inh.connect(condition='i==j', p=1) 
        self_inh.w = 1
        self_inh.u_f = U
        self_inh.x_d = 1
        
        # Remaining connections
        
        # Excitatory - Excitatory connections (with additional STDP)
        con_ee = Synapses(Pe, Pe, model = '''
                                            w_STDP : 1
                                            dApre/dt=-Apre/tau_pre : 1 (event-driven)
                                            dApost/dt=-Apost/tau_post : 1 (event-driven)
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            s_nmda_in2_post = int(not removed) * tmp_nmda_pre : 1 (summed)
                                            s_ampa_in2_post = int(not removed) * tmp_ampa_pre : 1 (summed)
                                            ''',
                          delay=condDelay,
                          on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_ampa_post += int(not removed) * w_STDP
                                    x_post += int(not removed) * w_STDP
                                    tmp_ampa = int(not removed) * clip(s_ampa*r_S*w_ee, 0, 100*r_S*w_ee)
                                    tmp_nmda = int(not removed) * clip(s_nmda*r_S*w_ee, 0, 100*r_S*w_ee)
                                    Apre += dApre
                                    w_STDP = int(not removed) * clip(w_STDP + Apost, 0, wmax)
                                    Ca += int(not removed) * a_Ca
                                    ''',
                         on_post='''
                                     Apost += dApost
                                     w_STDP = int(not removed) * clip(w_STDP + Apre, 0, wmax)
                                ''')
        
        con_ee.connect(condition='i!=j', p=ConnectionProbability_ee)
        con_ee.u_f = U
        con_ee.x_d = 1
        con_ee.w_STDP = w_ee
        
        # Excitatory - Inhibitory connections (STP only)
        con_ei = Synapses(Pe, Pi, model = '''
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            s_nmda_in3_post = int(not removed) * tmp_nmda_pre : 1 (summed)
                                            s_ampa_in3_post = int(not removed) * tmp_ampa_pre : 1 (summed)
                                            w : 1
                                            ''',
                          delay=condDelay,
                          on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_ampa_post += int(not removed) * w
                                    x_post += int(not removed) * w
                                    tmp_ampa = int(not removed) * clip(s_ampa*r_S*w_ei, 0, 100*r_S*w_ei)
                                    tmp_nmda = int(not removed) * clip(s_nmda*r_S*w_ei, 0, 100*r_S*w_ei)
                                    Ca += int(not removed) * a_Ca
                           ''')
        
        con_ei.connect(condition='i!=j', p=ConnectionProbability_ei)
        con_ei.u_f = U
        con_ei.x_d = 1
        con_ei.w = w_ei
        
        # Inhibitory - Excitatory connections (STP only)
        con_ie = Synapses(Pi, Pe, model = '''
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            w : 1
                                            ''', 
                          delay=condDelay,
                          on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_gaba_post += int(not removed) * r_S * w
                                    Ca += int(not removed) * a_Ca
                           ''')
        
        con_ie.connect(condition='i!=j', p=ConnectionProbability_ei)
        con_ie.u_f = U
        con_ie.x_d = 1
        con_ie.w = w_ie
        
        # Inhibitory - Inhibitory connections (STP only)
        con_ii = Synapses(Pi, Pi, model = '''
                                            du_f/dt = (U-u_f)/tau_f : 1 (event-driven)
                                            dx_d/dt = (1-x_d)/tau_d : 1 (event-driven)
                                            w : 1
                                            ''', 
                          delay=condDelay,
                          on_pre='''
                                    u_f += int(not removed) * U * (1 - u_f)
                                    r_S = int(not removed) * u_f * x_d
                                    x_d -= int(not removed) * r_S
                                    s_gaba_post += int(not removed) * r_S * w
                                    Ca += int(not removed) * a_Ca
                           ''')
        
        con_ii.connect(condition='i!=j', p=ConnectionProbability_ii)
        con_ii.u_f = U
        con_ii.x_d = 1
        con_ii.w = w_ii
        
        
        #%% Initial Values
        
        # Excitatory initial values
        Pe.v = 'El_e + rand()*(Vth_e-El_e)'  # random initial voltage
        Pe.s_ampa = 0
        Pe.s_gaba = 0
        Pe.s_nmda = 0
        Pe.Ca = 0
        Pe.I = I_on * np.random.normal(mu, sigma, 1).item()
        
        # Inhibitory initial values
        Pi.v = 'El_i + rand()*(Vth_i-El_i)' 
        Pi.s_gaba = 0
        Pi.Ca = 0
        Pi.s_ampa = 0
        Pi.s_nmda = 0
        Pi.I = I_on * np.random.normal(mu, sigma, 1).item()
        
        
        #%% Set up monitors
        
        exc_spikemon = SpikeMonitor(Pe) # exc neurons' spikes
        inh_spikemon = SpikeMonitor(Pi) # inh neurons' spikes
        
        statemon_STDP = StateMonitor(con_ee,'w_STDP',record=True,dt=10*ms) # w_STDP is the STDP-dependent exc-exc weight
               
        #%% Set up min and max bounds for the synaptic weight variables
        
        scale = ''
        
        for receptor in ['nmda','ampa']:
            for pop in ['Pe', 'Pi']:
                scale += pop+'.s_'+receptor+'_ = clip('+pop+'.s_'+receptor+', 0, 100)'
                scale += '\n'
        
        @network_operation()
        def s():
            exec(scale) # like an eval() statement in Matlab
        
        #%% Run pre-injury
        print("Running pre-injury phase...")    
        run(preinjury_simtime, report='text')
        										
		
        #%% Run control
        print("Running control phase...")
        run(simtime, report='text')
        store('control_phase')
        
        #%% Run injury
           
        inj_simtime = simtime # Simulation time for each injury round
        postinj_simtime = simtime
        
        # BDNF recovers 50% of injured inhibitory neurons at the second timepoint after control phase
        BDNF = [[0, 0.5], # injury then BDNF
                [0, 0], # injury only
                [0, 0]] # contol

        # glutamate injures 30% of exc neurons, 25% of inh neurons at the first timepoint after control phase
        glutamate = [[0.3, 0.25, 0, 0], # injury then BDNF
                     [0.3, 0.25, 0, 0], # injury only
                     [0, 0, 0, 0]] # control
        
        # glutamate also injures 75% of inh synapses at the first timepoint after control phase
        glutamateConnections = [[0.75, 0], # injury then BDNF
                                [0.75, 0], # injury only
                                [0, 0]] # control
        
        for ii in range(len(BDNF)): # for the number of experimental conditions
            restore('control_phase')
              																													 
			# glutamate injures neurons (exc & inh)																																																																					
            exc_injuredNeurons1 = np.random.choice(NE,size = int(floor(NE*glutamate[ii][0])), replace = 0) # apply glutamate at timepoint 1
            inh_injuredNeurons1 = np.random.choice(NI,size = int(floor(NI*glutamate[ii][1])), replace = 0) # apply glutamate at timepoint 1
            for neuron in exc_injuredNeurons1: # exc neurons to injure
                Pe.removed[int(neuron)] = "True" # injure (remove from simulation)
            for neuron in inh_injuredNeurons1: # inh neurons to injure
                Pi.removed[int(neuron)] = "True" # injure (remove from simulation)
			
            # glutamate injures inh-directed connections
            nConnections_ie = len(con_ie.i)
            ie_injuredConnections1 = np.random.choice(nConnections_ie,size = int(floor(nConnections_ie*glutamateConnections[ii][0])), replace = 0) # apply glutamate to connections at timepoint 1
            con_ie.w[ie_injuredConnections1] = np.repeat(0,len(ie_injuredConnections1)) # set the injured connections' weight to 0
                
            nConnections_ii = len(con_ii.i)
            ii_injuredConnections1 = np.random.choice(nConnections_ii,size = int(floor(nConnections_ii*glutamateConnections[ii][0])), replace = 0) # apply glutamate to connections at timepoint 1
            con_ii.w[ii_injuredConnections1] = np.repeat(0,len(ii_injuredConnections1)) # set the injured connections' weight to 0
            
            # Run injury
            print(f"Running timepoint 1 for simulation {ii+1} of {len(glutamate)}")
            run(inj_simtime,report = 'text')
            
            
            # apply BDNF if appropriate
            inh_resurrectedNeurons = np.random.choice(inh_injuredNeurons1, size = int(floor(BDNF[ii][1]*NI*glutamate[ii][1])), replace = 0)
            for neuron in inh_resurrectedNeurons:
                Pi.removed[int(neuron)] = "False" # resurrect (reintroduce to simulation along with pre-injury connections)
                
                
			# Run BDNF							
            print(f"Running timepoint 2 for simulation {ii+1} of {len(glutamate)}")
            run(inj_simtime,report = 'text')
            
            
			# Run Analysis				   
            print(f"Running analysis period")
            run(postinj_simtime,report = 'text')
            
            #%% Save data:
                                       		
            variables_to_save = {
                "TimeStepInSeconds" : float(defaultclock.dt),
                "SimTimeInSeconds" : float(simtime),
                "PreInjurySimTimeInSeconds" : float(preinjury_simtime),
                "InjurySimTimeInSeconds" : float(inj_simtime),
                "PostInjurySimTimeInSeconds" : float(postinj_simtime),
                "CellDensity" : CellDensity,
                "RecordingWidth" : RecordingWidth,
                "ConnectionDensityEE" : connectionDensityEE,
                "ConnectionDensityEI" : connectionDensityEI,
                "ConnectionDensityII" : connectionDensityII,
                "EIBalance" : EIBalance,
                "N" : N,
                "NE" : NE,
                "NI" : NI,

                "exc_spiketimes" : np.fromiter(map(float,exc_spikemon.t),float), # convert from Brian objects to floats
                "inh_spiketimes" : np.fromiter(map(float,inh_spikemon.t),float),
                "exc_spikeindexes" : exc_spikemon.i[:],
                "inh_spikeindexes" : inh_spikemon.i[:],
                
                "w_STDP_statemon" : statemon_STDP.w_STDP[:],
                "w_STDP_timesInSeconds" : statemon_STDP.t[:],
							  
                "w_ee" : w_ee,
                			
                "tau_Ca" : tau_Ca,	
                "g_AHP" : g_AHP,
                "tau_ampa" : tau_ampa,
		   
                "exc_injuredNeurons1" : exc_injuredNeurons1[:],
                "inh_injuredNeurons1" : inh_injuredNeurons1[:],
                
                "ie_injuredConnections1" : ie_injuredConnections1[:],
                "ii_injuredConnections1" : ii_injuredConnections1[:],
                
                "inh_resurrectedNeurons" : inh_resurrectedNeurons[:],
				                										  
                "BDNF" : BDNF[ii],
                "glutamate" : glutamate[ii],
                "glutamateConnections" : glutamateConnections[ii],
                
                "NeuronXPosition" : NeuronXPosition[:],
                "NeuronYPosition" : NeuronYPosition[:],
                
                "nConnections_ee" : len(con_ee.i),
                "nConnections_ei" : len(con_ei.i),
                "nConnections_ie" : len(con_ie.i),                                 
                "nConnections_ii" : len(con_ii.i),

                "con_ee_i" : con_ee.i[:], # ID connection partners for each connection type
                "con_ee_j" : con_ee.j[:],

                "con_ei_i" : con_ei.i[:],
                "con_ei_j" : con_ei.j[:],
                
                "con_ii_i" : con_ii.i[:],
                "con_ii_j" : con_ii.j[:],
                
                "con_ie_i" : con_ie.i[:],
                "con_ie_j" : con_ie.j[:]}
             
            savename = str.format("BDNF Homeostasis - BDNF T1- {}, T2- {}, Glutamate T1- E {} I {} I Conn {}, T2- E {} I {} I Conn {} Saved {}.mat", BDNF[ii][0], BDNF[ii][1], glutamate[ii][0],glutamate[ii][1], glutamateConnections[ii][0], glutamate[ii][2], glutamate[ii][3], glutamateConnections[ii][1], datetimeToFileLabelString())
            sio.savemat(savename,variables_to_save)
            print('Saved')
            
