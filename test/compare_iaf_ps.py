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
# Parameters
#------------------------
#

num_threads = 6

models = [ "iaf_psc_alpha_canon", "ps_iaf_psc_alpha" ]
num_models = len(models)
g_L = 12.
tau_m = 8.
lst_param = [ {"tau_syn": 0.2, "tau_m": 8., "C_m": tau_m*g_L }, {'g_L': g_L, "C_m": tau_m*g_L}, {'g_L': 12.} ]
lst_record = [ "y2", "I_syn_ex", "g_ex" ]
lst_ylabel = [ "I (pA)", "I (pA)", 'Conductance (nS)' ]

di_param = {
    'V_reset': -65.,
    'V_th': -50.,
    'I_e': 0.0,
    'E_L': -60.,
    'C_m': 100.,
    'V_m': -60.,
}

individual = False
network = True
#~ individual = True
#~ network = False


#-----------------------------------------------------------------------------#
# Individual dynamics
#------------------------
#

nest.ResetKernel()
if individual:
    r_resolution = 0.01
    nest.SetKernelStatus({"resolution":r_resolution})
    d_step_current = 160.
    r_min_voltage = -70.

    # models
    lst_neurons = [ nest.Create(model,params=di_param) for model in models ]

    #~ step_gen = nest.Create("step_current_generator",1,{"amplitude_times": [50.,1500.], "amplitude_values":[d_step_current,0.]})
    dc = nest.Create("dc_generator", params={"amplitude": 200.})
    pg = nest.Create("poisson_generator_ps", params={"rate": 100.})
    pn = nest.Create("parrot_neuron_ps")
    nest.Connect(pg,pn)
    multimeter = nest.Create("multimeter",num_models)

    nest.Connect([lst_neurons[0][0], lst_neurons[1][0]], [lst_neurons[0][0], lst_neurons[1][0]], conn_spec={'rule': 'all_to_all', 'autapses': False}, syn_spec={'weight':30.})

    for i,neuron in enumerate(lst_neurons):
        #~ nest.Connect(step_gen,neuron)
        nest.Connect(dc,neuron)
        nest.SetStatus(neuron, lst_param[i])
        print nest.GetStatus(neuron)
        nest.SetStatus((multimeter[i],), {"withtime":True, "interval":r_resolution, "record_from":["V_m", lst_record[i]]})
        nest.Connect(multimeter[i],neuron[0])
        nest.Connect(pn,(neuron[0],))

    nest.Simulate(1600.0)

    # plot

    plt.close("all")
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)

    # get the neuron's membrane potential
    for i in range(num_models):
        dmm = nest.GetStatus(multimeter)[i]
        da_voltage = dmm["events"]["V_m"]
        da_adapt = dmm["events"][lst_record[i]]
        da_time = dmm["events"]["times"]
        ax1.plot(da_time,da_voltage,c=cm.hot(i/float(num_models)), label=models[i])
        ax1.set_ylabel('Voltage (mV)')
        ax2.plot(da_time,da_adapt,c=cm.hot(i/float(num_models)), label=models[i])
        ax2.set_xlabel('Time (ms)')
        ax2.set_ylabel(lst_ylabel[i])

    plt.legend(loc=4)


#-----------------------------------------------------------------------------#
# Compare network dynamics
#------------------------
#

if network:
    # time the simulations for each neural model and network size
    sim_time = 1000.
    r_resolution = 0.1
    #~ lst_network_sizes = np.arange(1000, 2000, 5000)
    lst_network_sizes = [2]
    num_runs = len(lst_network_sizes)
    lst_times = [ np.zeros(num_runs) for _ in range(len(models)) ]
    lst_spikes = [ np.zeros(num_runs) for _ in range(len(models)) ]

    for i,size in enumerate(lst_network_sizes):
        fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
        for j,model in enumerate(models):
            nest.ResetKernel()
            nest.SetKernelStatus({"resolution":r_resolution})
            nest.SetKernelStatus({"local_num_threads": num_threads})
            # synaptic weight
            weight = 30.
            gids = nest.Create(model, size, params=di_param)
            nest.SetStatus(gids, lst_param[j])
            nest.Connect(gids, gids, conn_spec={'rule': 'all_to_all', 'autapses': False}, syn_spec={'weight':weight})
            multimeter = nest.Create("multimeter")
            nest.SetStatus(multimeter, {"withtime":True, "interval":r_resolution, "record_from":["V_m", lst_record[j]]})
            # in nest
            dc = nest.Create("dc_generator", params={"amplitude": 200.})
            sd = nest.Create("spike_detector", params={"withtime": True, "withgid": True, "precise_times": True})
            nest.Connect(dc, gids)
            nest.Connect(gids,sd)
            nest.Connect(multimeter, [gids[0]])

            start = time.time()
            nest.Simulate(sim_time)
            lst_times[j][i] = time.time() - start
            lst_spikes[j][i] = len(nest.GetStatus(sd)[0]["events"]["senders"]) / sim_time

            dmm = nest.GetStatus(multimeter)[0]
            da_voltage = dmm["events"]["V_m"]
            da_adapt = dmm["events"][lst_record[j]]
            da_time = dmm["events"]["times"]
            ax1.plot(da_time,da_voltage,c=cm.hot(j/float(num_models)), label=models[j])
            ax1.set_ylabel('Voltage (mV)')
            ax2.plot(da_time,da_adapt,c=cm.hot(j/float(num_models)), label=models[j])
            ax2.set_xlabel('Time (ms)')
            ax2.set_ylabel(lst_ylabel[j])

            from_device(sd, title="Raster for {} neurons of type {}".format(size, model))

    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True)
    for i, model in enumerate(models):
        ax1.plot(lst_network_sizes, lst_times[i],  c=cm.hot(i/float(num_models)), label=model)
        ax2.scatter(lst_network_sizes, lst_spikes[i], c=cm.hot(i/float(num_models)), label=model)
        ax2.set_xlabel("Network size")
        ax1.set_ylabel("Simulation time")
        ax2.set_ylabel("Number of spikes")
        ax1.legend(loc=2)
        ax2.legend(loc=2)

plt.show()
