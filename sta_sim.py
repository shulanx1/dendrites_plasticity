#!/usr/bin/env python3
"""
Simulate model with random input and compute synaptic influence before each spike.
"""
#%%
import sys
import os
wd = 'E:\\Code\\dendrites_plasticity' # working directory
sys.path.insert(1, wd)

import numpy as np
import pickle
import sys
import time
import os
from dendrites import comp_model
from dendrites import neuron_model
from dendrites import plot_raster
from dendrites import sequences
from dendrites import training


from neuron import h

seed = int(time.time())
model = 'l5_act_b'
results_file = wd + '\\outputs\\sta\\sta_'+model+'_data\\sta_'+model+'_'+str(seed)
path = wd + '\\outputs\\sta\\sta_'+model+'_data'
if not os.path.exists(path):
    os.makedirs(path)


### Simulation parameters ###
T = 10000            # simulation time (ms)
dt = 0.1            # time step (ms)
v_init = -75        # initial voltage (mV)
reps = 3           # number of reps
w_jitter = 0.5      # perturbation to initial weights
t_window = 150	    # analysis window (ms)
offset = 2		    # spike time offset (ms)
mu = 0              # mean parameter for lognormal rate dist
sigma = 1.5           # sd parameter for lognormal rate dist
np.random.seed(int(seed))

ion_specific = True


def init_rand_sequence(rates_e, rates_i, T):
    """
    build sequences of Poisson presynaptic input to excitatory and inhibitory
    synapses

    Parameters
    ----------
    rates_e, rates_i : list
        excitatory and inhibitory input rates
    T : int
        total simulation time (ms)
    Returns
    -------
    S_e, S_i :  list
        excitatory and inhibitory presynaptic spike times
    """
    S_e = sequences.build_rate_seq(rates_e, 0, T)
    S_i = sequences.build_rate_seq(rates_i, 0, T)
    return S_e, S_i


def init_weights(P):
    """
    Initialise synaptic weights by perturbing initial values.

    Parameters
    ----------
    P : dict
        model parameters

    Returns
    -------
    w_e, w_i :  ndarray
        excitatory and inhibitory weight vectors
    """
    w_e = np.ones(P['N_e'])*[P['g_max_A'] + P['g_max_N']] + w_jitter*(P['g_max_A'] +
            P['g_max_N'])*(np.random.rand(P['N_e']) - 1/2)
    w_i = np.ones(P['N_i'])*[P['g_max_G']] + w_jitter*P['g_max_G']*(
            np.random.rand(P['N_i']) - 1/2)
    return w_e, w_i


def spike_times(dt, v):
    """ Get spike times from voltage trace.

    Parameters
    ----------
    dt : float
        simulation timestep
    v : ndarray
        compartment voltages v=v[compartment, time]
    Returns
    -------
    t_spike : ndarray
        spike times
    """
    thresh_cross = np.where(v[0, :] > 0)[0]
    if thresh_cross.size > 0:
        spikes = np.where(np.diff(thresh_cross) > 1)[0] + 1
        spikes = np.insert(spikes, 0, 0)
        spikes = thresh_cross[spikes]
        t_spike = spikes*dt - offset
    else:
        t_spike = np.array([])
    return t_spike

def spike_times_l5(dt, v):
    """ Get spike times from voltage trace.

    Parameters
    ----------
    dt : float
        simulation timestep
    v : ndarray
        compartment voltages v=v[compartment, time]
    Returns
    -------
    t_spike : ndarray
        spike times
    """
    thresh_cross = np.where(v[1, :] > 0)[0]
    if thresh_cross.size > 0:
        spikes = np.where(np.diff(thresh_cross) > 1)[0] + 1
        spikes = np.insert(spikes, 0, 0)
        spikes = thresh_cross[spikes]
        t_spike = spikes*dt - 2
    else:
        t_spike = np.array([])
    return t_spike

def get_grad(cell, t0, t1, dt, S_e, S_i, soln, stim):
    """ Get gradients (dv_soma/dw) for individual synaptic activations by
    solving variational equations between times t0 and t1 with fast somatic
    conductances set to zero.

    Parameters
    ----------
    cell : dendrites.comp_model.CModel object
        compartmental model instance used for simulation
    t0, t1 :  float
        initial and final times for computing gradients
    dt : float
        simulation time step
    S_e, S_i : array_like
        input spike patterns for E and I synapses
    soln : list
        model states [v, m, h, n, p] (output from cell.simulate)
    stim : list
        synapse indices and states [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

    Returns
    -------
    v_pre : list
        voltage throughout the morphology preceding somatic spike
    F_e, F_i : list
        computed gradients with associated synapse indices and presynaptic
        spike times F_e = [f_e, z_ind_e, z_e]
    """
    IC_sub = cell.set_IC(soln, stim, int(t0/dt))
    g_na_temp, g_k_temp = P['g_na'], P['g_k']
    cell.P['g_na'] = 0
    cell.P['g_k'] = 0
    t_s, soln_s, stim_s = cell.simulate(t0, t1, dt, IC_sub, S_e, S_i)
    Z_e = sequences.subsequence(S_e, t0, t1)
    Z_i = sequences.subsequence(S_i, t0, t1)
    Z_e, z_ind_e = sequences.rate2temp(Z_e)
    Z_i, z_ind_i = sequences.rate2temp(Z_i)
    f_e, f_i = cell.grad_w(soln_s, stim_s, t_s, dt, Z_e, Z_i, z_ind_e, z_ind_i)
    cell.P['g_na'] = g_na_temp
    cell.P['g_k'] = g_k_temp
    v_pre = soln[0][:, int(t1 / dt)]
    f_e = f_e[:, -1]
    f_i = f_i[:, -1]
    z_e = Z_e - t1
    z_i = Z_i - t1
    return [v_pre], [f_e, z_ind_e, z_e], [f_i, z_ind_i, z_i]

def get_grad_l5(cell, t0, t1, dt, S_e, S_i, soln, stim):
    """ Get gradients (dv_soma/dw) for individual synaptic activations by
    solving variational equations between times t0 and t1 with fast somatic
    conductances set to zero.

    Parameters
    ----------
    cell : dendrites.comp_model.CModel object
        compartmental model instance used for simulation
    t0, t1 :  float
        initial and final times for computing gradients
    dt : float
        simulation time step
    S_e, S_i : array_like
        input spike patterns for E and I synapses
    soln : list
        model states [v, m, h, n, p, hcn] (output from cell.simulate)
    stim : list
        synapse indices and states [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

    Returns
    -------
    F_e, F_i : list
        computed gradients with associated synapse indices and presynaptic
        spike times F_e = [f_e, z_ind_e, z_e]
    """
    IC_sub = cell.set_IC(soln, stim, int(t0 / dt))
    # g_na_temp, g_k_temp = P['g_na'], P['g_k']
    # cell.P['g_na'] = 0
    # cell.P['g_k'] = 0
    # g_ion_temp = cell.g_ion
    # cell.g_ion[0][1] = 0.0
    # cell.g_ion[1][1] = 0.0
    t_s, soln_s, stim_s = cell.simulate_L5(t0, t1, dt, IC_sub, S_e, S_i)
    Z_e = sequences.subsequence(S_e, t0, t1)
    Z_i = sequences.subsequence(S_i, t0, t1)
    Z_e, z_ind_e = sequences.rate2temp(Z_e)
    Z_i, z_ind_i = sequences.rate2temp(Z_i)
    f_e, f_i = cell.grad_w_l5(soln_s, stim_s, t_s, dt, Z_e, Z_i, z_ind_e, z_ind_i)
    # cell.g_ion = g_ion_temp

    v_pre = soln[0][:, int(t1 / dt)]
    f_e = f_e[:, -1]
    f_i = f_i[:, -1]
    z_e = Z_e - t1
    z_i = Z_i - t1
    return [v_pre], [f_e, z_ind_e, z_e], [f_i, z_ind_i, z_i]

def get_grad_l5_channels(cell, t0, t1, dt, S_e, S_i, soln, stim):
    """ Get gradients (dv_soma/dw) and channel gradient for individual synaptic activations by
    solving variational equations between times t0 and t1 with fast somatic
    conductances set to zero.

    Parameters
    ----------
    cell : dendrites.comp_model.CModel object
        compartmental model instance used for simulation
    t0, t1 :  float
        initial and final times for computing gradients
    dt : float
        simulation time step
    S_e, S_i : array_like
        input spike patterns for E and I synapses
    soln : list
        model states [v, m, h, n, p, hcn] (output from cell.simulate)
    stim : list
        synapse indices and states [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

    Returns
    -------
    F_e, F_i : list
        computed gradients with associated synapse indices and presynaptic
        spike times F_e = [f_e, z_ind_e, z_e]
    """
    IC_sub = cell.set_IC(soln, stim, int(t0 / dt))
    t_s, soln_s, stim_s = cell.simulate_L5(t0, t1, dt, IC_sub, S_e, S_i)
    Z_e = sequences.subsequence(S_e, t0, t1)
    Z_i = sequences.subsequence(S_i, t0, t1)
    Z_e, z_ind_e = sequences.rate2temp(Z_e)
    Z_i, z_ind_i = sequences.rate2temp(Z_i)
    f_e, f_i, c_e, c_i,  = cell.grad_w_l5_channels(soln_s, stim_s, t_s, dt, Z_e, Z_i, z_ind_e, z_ind_i)

    v_pre = soln[0][:, int(t1 / dt)]
    f_e = f_e[[0, 1, 800], :, -1]
    f_i = f_i[[0, 1, 800], :, -1]
    c_e = c_e[:, [0, 1, 800], :, -1]
    c_i = c_i[:, [0, 1, 800], :, -1]
    z_e = Z_e - t1
    z_i = Z_i - t1
    return [v_pre], [f_e, c_e, z_ind_e, z_e], [f_i, c_i, z_ind_i, z_i]

# P = params.init_params()
from dendrites import parametersL5_Hay
P = parametersL5_Hay.init_params(wd)

h('forall pop_section()')
h('forall delete_section()')

cell = comp_model.CModel(P, verbool = False)
# t, soln, stim = cell.simulate(0, T, dt, v_init, S_e, S_i)
# np.random.seed(int(seed))
print("initialized")
F_e = []
F_i = []
V = []
W_e = []
W_i = []

#%%
for rep in range(reps):
    print('starting rep %d' % (rep+1))
    rates_e, rates_i = sequences.lognormal_rates(1, P['N_e'], P['N_i'], mu, sigma)
    w_e, w_i = init_weights(P)
    cell.set_weights(w_e, w_i)
    S_e, S_i = init_rand_sequence(rates_e[0], rates_i[0], T)
    print('simulating rep %d' % (rep+1))
    t, soln, stim = cell.simulate_L5(0, T, dt, v_init, S_e, S_i)
    v = soln[0]
    t_spikes = spike_times_l5(dt, v)
    num_spikes = len(t_spikes)
    print('solve the gradients for rep %d' % (rep+1))
    if ion_specific:
        t_window = 100
        if num_spikes > 0:
            incl = np.where(np.diff(np.insert(t_spikes, 0, 0)) > t_window)[0]
            t1 = t_spikes[incl]
            t0 = t1 - t_window
            num_test_spikes = len(t1)
            print('%d spikes in total' % num_test_spikes)
            for k, spike in enumerate(t1):
                v_pre, f_e, f_i = get_grad_l5_channels(cell, t0[k], t1[k], dt, S_e, S_i, soln, stim)
                V.append(v_pre)
                F_e.append(f_e)
                F_i.append(f_i)
                W_e.append(np.array(w_e[f_e[2]]))
                W_i.append(np.array(w_i[f_i[2]]))
                print('gradients solved for spike %d in rep %d' % (k + 1, rep + 1))
    
        if seed == '0':
            pickle.dump([cell, V, F_e, F_i, W_e, W_i], open(results_file + '_channels', 'wb'))
            pickle.dump(soln[0], open(results_file + '_channels_v', 'wb'))
        else:
            pickle.dump([cell, V, F_e, F_i, W_e, W_i], open(results_file + '_channels', 'wb'))
            pickle.dump(soln[0], open(results_file + '_channels_v', 'wb'))

        del V
        del F_e
        del F_i
        del soln
        import gc
        gc.collect()
    else:
        if num_spikes > 0:
            incl = np.where(np.diff(np.insert(t_spikes, 0, 0)) > t_window)[0]
            t1 = t_spikes[incl]
            t0 = t1 - t_window
            num_test_spikes = len(t1)
            print('%d spikes in total' % num_test_spikes)
            for k, spike in enumerate(t1):
                v_pre, f_e, f_i = get_grad_l5(cell, t0[k], t1[k], dt, S_e, S_i, soln, stim)
                # f_e, f_i = get_k_grad(t0[spike], t1[spike], dt, S_e, S_i, v,
                #                       kernel_params)
                V.append(v_pre)
                F_e.append(f_e)
                F_i.append(f_i)
                W_e.append(np.array(w_e[f_e[1]]))
                W_i.append(np.array(w_i[f_i[1]]))
                print('gradients solved for spike %d in rep %d' % (k + 1, rep + 1))
    
        if seed == '0':
            pickle.dump([cell, V, F_e, F_i, W_e, W_i], open(results_file, 'wb'))
        else:
            pickle.dump([cell, V, F_e, F_i, W_e, W_i], open(results_file, 'wb'))

        del V
        del F_e
        del F_i
        del soln
        import gc
        gc.collect()
