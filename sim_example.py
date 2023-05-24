#!/usr/bin/env python3
"""
Example script to demonstrate simulation of the model and calculation of
dv_soma/dw.
"""
#%%
import sys
import os
wd = 'E:\\Code\\dendrites_pasticity' # working directory
sys.path.insert(1, wd)

import numpy as np
import pickle

from matplotlib import pyplot
import matplotlib.pyplot as plt
from neuron import h
import neuron

from dendrites import comp_model
from dendrites import neuron_model
from dendrites import parameters1
from dendrites import plot_raster
from dendrites import sequences
from dendrites import training
from dendrites import parametersL5_Hay
P = parametersL5_Hay.init_params(wd)


T = 200        # simulation time (ms)
dt = 0.1        # time step (ms)
v_init = -75    # initial voltage (mV)
seed = 1        # random seed
stim_dur = 300							# stimulus duration
stim_on = 100							# stimulus on
stim_off = stim_on + stim_dur           # stimulus off
t_on = 0								# background on
t_off = stim_on							# background off
r_0 = 1.25								# background rate
dt = 0.1            					# time step
r_mean = 2.5
num_patterns = 4
input = 'opt'
param_sets = {'rate':[40., 0, 0.], 'temp':[2.5, 1, 1.], 'opt':[20., 1, 1.]}
r_max, num_t, s = param_sets[input]
mu = 0
sigma = 2
w_jitter = 0.5      # perturbation to initial weights


kernel_fit = wd + "\\input\\kernel_fit_act1"  # fitted plasticity kernel
# P = parameters1.init_params(wd)            # stores model parameters in dict P


c_model = True # True for custom compartmental model with explicit gradient
                 # calculations, False for Neuron model using fitted approximation
num_spikes = 3   # max number of somatic spikes for which to compute gradients
np.random.seed(seed)


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
    z_e = Z_e - t1
    z_i = Z_i - t1
    F_e = [f_e[:, -1], z_ind_e, z_e]
    F_i = [f_i[:, -1], z_ind_i, z_i]
    return F_e, F_i


def get_k_grad(t0, t1, dt, S_e, S_i, v, kernel_params):
    """ Get gradients (dv_soma/dw) for individual synaptic activations from
    fitted approximations.

    Parameters
    ----------
    t0, t1 :  float
        initial and final times for computing gradients
    dt : float
        simulation time step
    S_e, S_i : array_like
        input spike patterns for E and I synapses
    v : ndarray
        compartment voltages from simulation
    kernel_params : list
        parameters for fitted kernels (see dendrites.kernels)

    Returns
    -------
    F_e, F_i : list
        approximated gradients with associated synapse indices and presynaptic
        spike times F_e = [f_e, z_ind_e, z_e]
    """
    Z_e = sequences.subsequence(S_e, t0, t1)
    Z_i = sequences.subsequence(S_i, t0, t1)
    Z_e, z_ind_e = sequences.rate2temp(Z_e)
    Z_i, z_ind_i = sequences.rate2temp(Z_i)
    f_e, f_i = training.kernel_grad_ss(int(t0 / dt), int(t1 / dt), S_e, S_i,
                        cell.seg_e, cell.seg_i, cell.b_type_e, cell.b_type_i, v,
                        kernel_params, dt)
    F_e = [f_e, z_ind_e, Z_e - t1]
    F_i = [f_i, z_ind_i, Z_i - t1]
    return F_e, F_i


def init_input(P, num_patterns, stim_on, stim_off, r_mean, r_max, num_t, s):
    """
    Initialise input rates and spike time sequences for feature-binding task.

    Parameters
    ----------
    P : dict
        model parameters
    num_patterns : int
        number of input patterns to be classified
    stim_on, stim_off : int
        time of stimulus onset and termination (ms)
    r_mean : float
        average presynaptic population rate (Hz)
    r_max : float
        time averaged input rate to active synapses
    num_t : int
        number of precisely timed events per active synapse
    s : float
        interpolates between rate (s=0) and temporal (s=1) input signals (mostly
        unused parameter -- to be removed)

    Returns
    -------
    rates_e, rates_i : list
        excitatory and inhibitory input rates for all patterns
    S_E, S_I : list
        times of precisely timed events for all patterns
    """
    N_e, N_i = P['N_e'], P['N_i']
    ind_e = np.arange(N_e)
    ind_i = np.arange(N_i)
    np.random.shuffle(ind_e)
    np.random.shuffle(ind_i)
    rates_e, rates_i = sequences.assoc_rates(num_patterns, N_e, N_i, r_mean,
                                             r_max)
    rates_e = [r[ind_e] for r in rates_e]
    rates_i = [r[ind_i] for r in rates_i]
    if s > 0:
        S_E, S_I = sequences.assoc_seqs(num_patterns, N_e, N_i, stim_on, stim_off,
                                        num_t)
        S_E = [s[ind_e] for s in S_E]
        S_I = [s[ind_i] for s in S_I]
        for s_e, r_e in zip(S_E, rates_e):
            s_e[r_e == 0] = np.inf
        for s_i, r_i in zip(S_I, rates_i):
            s_i[r_i == 0] = np.inf
    else:
        S_E, S_I = sequences.build_seqs(num_patterns, N_e, N_i, stim_on, stim_off,
                                        0)
    return rates_e, rates_i, S_E, S_I

def pad_S(S0):
    l = np.max([len(s) for s in S0])
    S = np.full((len(S0), l), np.inf)
    for k, s in enumerate(S0):
        S[k, :len(s)] = s
    return S

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

#%%
h('forall pop_section()')
h('forall delete_section()')
rates_e, rates_i = sequences.lognormal_rates(1, P['N_e'], P['N_i'], mu, sigma)
w_e, w_i = init_weights(P)
S_e, S_i = init_rand_sequence(rates_e[0], rates_i[0], T)
if c_model:
    cell_comp = comp_model.CModel(P, verbool = True)
    t, soln, stim = cell_comp.simulate_L5(0, T, dt, v_init, S_e, S_i)
    v = soln[0]
    gates = soln[1]
    plt.plot(t, v[1])
else:
    cell = neuron_model.NModel(P, verbool = True)
    # cell.set_deficite_channels('SK_E2', sec_name = 'all',  percentage = 0)
    # cell.set_deficite_channels('CaDynamics_E2', sec_name = 'all',  percentage = 0)
    cell.set_weights(w_e, w_i)
    _, kernel_params = pickle.load(open(kernel_fit, 'rb'))
    cell.kernel = kernel_params
    t, v = cell.simulate(T, dt, v_init, S_e, S_i)
    plt.plot(t, v[1])

#%%### Compute Gradients ###
t_window = 100  # synaptic plasticity window (fixed parameter)
t1 = spike_times(dt, v)
if len(t1) > 0:
    t0 = t1 - t_window
    E_data = []
    I_data = []
    for spike in range(min(num_spikes, len(t1))):
        if c_model:
            F_e, F_i = get_grad(cell, t0[spike], t1[spike], dt, S_e,
                                S_i, soln, stim)
        else:
            F_e, F_i = get_k_grad(t0[spike], t1[spike], dt, S_e, S_i, v,
                                  kernel_params)
        E_data.append(F_e)
        I_data.append(F_i)

### Plot Results ###
fig, ax = pyplot.subplots(figsize=(8, 2.5))
ax.plot(t, v[0, :], 'k')
ax.set_xlim([0, T])
ax.set_ylim([-80, 40])
ax.set_yticks(np.arange(-75, 50, 25))
ax.set_xlabel('time (ms)', fontsize=14)
ax.set_ylabel('V' + r'$_{soma}$' + ' (mV)', fontsize=14)
ax.tick_params(labelsize=12)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
pyplot.tight_layout()

if len(t1) > 0:
    i_positions, index, boundaries = plot_raster.raster_params(cell)
    for e_dat, i_dat in zip(E_data, I_data):
        f_e, e_ind, s_e = e_dat
        f_i, i_ind, s_i = i_dat
        ff_e = np.array(np.abs(f_e)) / np.max(np.abs(f_e))
        ff_i = np.array(f_i)
        if np.min(ff_i) < 0:
            ff_i = -ff_i / np.min(ff_i)
        plot_raster.plot_grad_example(ff_e, ff_i, e_ind, i_ind, s_e, s_i,
                                     t_window, i_positions, index, boundaries)

pyplot.show()
