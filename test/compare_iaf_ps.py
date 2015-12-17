#!/usr/bin/env python
#-*- coding:utf-8 -*-

""" Script to test the precise iaf behaviours """

import argparse
import time
import numpy as np

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib import rcParams

import nest
from nest.raster_plot import from_device
nest.Install("nngt_module")



#-----------------------------------------------------------------------------#
# Parser
#------------------------
#

parser = argparse.ArgumentParser(description="Script to compare the grid-precise and usual models.", usage='%(prog)s [options]')
parser.add_argument("-i", "--indivdual", action="store", default=True,
          help="Compare the individual-neuron dynamics.")
parser.add_argument("-nn", "--no_network", action="store_true",
          help="Compare network dynamics.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-nt", "--notime", action="store_true",
          help="Do not compare runtime.")
group.add_argument("-s", "--size", action="store", default=5000,
          help="Compare for a given network size.")

## parse
args = parser.parse_args()


#-----------------------------------------------------------------------------#
# Individual dynamics
#------------------------
#

nest.ResetKernel()
nest.SetKernelStatus({"local_num_threads": 6})
if args.indivdual:
    r_resolution = 0.01
    nest.SetKernelStatus({"resolution":r_resolution})
    d_step_current = 160.
    r_min_voltage = -70.

    # create AdExp neurons
    di_param = {
        'V_reset': -65.,
        'V_th': -50.,
        'I_e': 0.0,
        'g_L': 12.,
        'E_L': -60.,
        'C_m': 100.,
        'V_m': -60.
    }

    # models
    models = [ "iaf_cond_alpha", "ps_iaf_cond_alpha" ]
    #~ models = [ "iaf_cond_alpha" ]
    lst_neurons = [ nest.Create(model,params=di_param) for model in models ]
    num_neurons = len(lst_neurons)


    step_gen = nest.Create("step_current_generator",1,{"amplitude_times": [50.,1500.], "amplitude_values":[d_step_current,0.]})
    pg = nest.Create("poisson_generator", params={"rate": 100.})
    pn = nest.Create("parrot_neuron")
    nest.Connect(pg,pn)
    multimeter = nest.Create("multimeter",num_neurons)
    nest.SetStatus(multimeter, {"withtime":True, "interval":r_resolution, "record_from":["V_m","g_ex"]})

    for i,neuron in enumerate(lst_neurons):
        nest.Connect(step_gen,neuron)
        nest.Connect(multimeter[i],neuron[0])
        nest.Connect(pn,(neuron[0],))

    nest.Simulate(1600.0)

    # plot

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)

    # get the neuron's membrane potential
    for i in range(num_neurons):
        dmm = nest.GetStatus(multimeter)[i]
        da_voltage = dmm["events"]["V_m"]
        da_adapt = dmm["events"]["g_ex"]
        da_time = dmm["events"]["times"]
        ax1.plot(da_time,da_voltage,c=cm.hot(i/float(num_neurons)), label=models[i])
        ax1.set_ylabel('Voltage (mV)')
        ax2.plot(da_time,da_adapt,c=cm.hot(i/float(num_neurons)), label=models[i])
        ax2.set_xlabel('Time (ms)')
        ax2.set_ylabel('Conductance (nS)')

    plt.legend(loc=4)


#-----------------------------------------------------------------------------#
# Compare network dynamics
#------------------------
#

#~ if not args.no_network:
if False:
    # time the simulations for each neural model and network size
    sim_time = 1000.
    lst_network_sizes = [args.size] if args.notime else np.arange(1000, 6000, 1000)
    num_runs = len(lst_network_sizes)
    lst_times = [ np.zeros(num_runs) for _ in range(len(models)) ]
    lst_spikes = [ np.zeros(num_runs) for _ in range(len(models)) ]
    # fraction of inhibitory neurons
    ifrac = 0.2
    # average degree
    avg_deg = 100

    graph, gids = None, None

    for i,size in enumerate(lst_network_sizes):
        for j,model in enumerate(models):
            nest.ResetKernel()
            nest.SetKernelStatus({"local_num_threads": 6})
            # synaptic weight
            weight = 30. if "exp" in model else 20.
            inhib_start = int(size*(1-ifrac))
            gids = nest.Create(model, size)
            nest.Connect(gids[:inhib_start], gids, conn_spec={'rule': 'fixed_indegree', 'indegree': int(avg_deg/2), 'autapses': False}, syn_spec={'weight':weight})
            nest.Connect(gids[inhib_start:], gids, conn_spec={'rule': 'fixed_indegree', 'indegree': int(avg_deg/2), 'autapses': False}, syn_spec={'weight':-weight})
            # in nest
            dc = nest.Create("dc_generator", params={"amplitude": 800.})
            sd = nest.Create("spike_detector", params={"withtime": True, "withgid": True})
            nest.Connect(dc, gids[:int(inhib_start/3)])
            nest.Connect(gids,sd)

            start = time.time()
            nest.Simulate(sim_time)
            lst_times[j][i] = time.time() - start
            lst_spikes[j][i] = len(nest.GetStatus(sd)[0]["events"]["senders"]) / sim_time

            #~ from_device(sd, title="Raster for {} neurons of type {}".format(size, model))

    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
    for i, model in enumerate(models):
        ax1.plot(lst_network_sizes, lst_times[i],  c=cm.hot(i/float(num_neurons)), label=model)
        ax2.scatter(lst_network_sizes, lst_spikes[i], c=cm.hot(i/float(num_neurons)), label=model)
        ax1.legend(loc=2)
        ax2.legend(loc=2)
        
plt.show()
