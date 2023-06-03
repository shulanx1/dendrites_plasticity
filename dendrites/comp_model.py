#!/usr/bin/env python3
"""
Compartmental model class and functions for simulation and training.
"""
import numpy as np
import numba as nb
import json

from dendrites import nak, Ih
from dendrites import morphology
from dendrites import neuron_model
from dendrites.l5_biophys import *

class CModel:
    """Compartmental model object for simulation of dendritic computation and
    learning.

    Parameters
    ----------
        P : dict
        model and simulation parameters

    Attributes
    ----------
        P : dict
            model and simulation paramaters
        L_s, a_s : ndarray
            segment lengths and radii
        C : ndarray
            adjacency matrix for segments
        Q : ndarray
            axial conductance matrix
        sec_e, sec_i : ndarray
            synapse locations
        seg_e, seg_i : ndarray
            synapse segment numbers
        seg2sec : ndarray
            sections and positions for each segment
        b_type_e, b_type_i : list
            branch types (basal/apical/soma) for all synapses
        H_e, H_i : ndarray
            projection matrices (inputs->compartments)
        w_e, w_i : ndarray
            synaptic weights
        g_ops : ndarray
            list of Gaussian elimination steps for solve
        f_ops : ndarray
            list of forward sub steps for solve
    """

    def __init__(self, P, verbool = False):
        self.P = P
        self.verbool = verbool
        A, nseg, L, a, self.sec_e, self.sec_i, self.b_type_e, self.b_type_i = \
        self.define_morphology()
        self.L_s, self.a_s, self.area = morphology.seg_geometry(L, a, nseg)
        if 'param_file' in P: #l5
            self.area[0] = 629.147/1e8
            self.L_s[0] = self.area[0]/(2*np.pi*self.a_s[0])
            [cm, gpas, gna, gkv, gcah, gcal, gih, gim, gnap, gkt, gkp, gsk,gamma, decay, dend, axon, apic, soma]  = self.insert_biophysical_L5()
            self.g_ion = np.vstack((gna, gkv, gcah, gcal, gih, gim, gnap, gkt, gkp, gsk))*1e3  #mS/cm^2
            self.dend = dend
            self.axon = axon
            self.apic = apic
            self.soma = soma
            self.cm = cm
            self.gamma = gamma
            self.decay = decay
        self.C, self.seg2sec = self.build_comp_matrix(A, nseg)
        self.Q = self.build_axial_mat()
        self.H_e, self.seg_e = self.syn_mat(self.sec_e)
        self.H_i, self.seg_i = self.syn_mat(self.sec_i)
        self.w_e = np.array(self.sec_e.shape[1]*[P['g_max_A'] + P['g_max_N']])
        self.w_i = np.array(self.sec_i.shape[1]*[P['g_max_G']])
        self.g_ops = self.gauss_ops()
        self.f_ops = self.fsub_ops()

    def define_morphology(self):
        """ Create adjacency matrix, lists of segment numbers and dimensions,
        and synapse locations.

        Returns
        -------
        A : ndarray
            adjacency matrix for dendritic sections
        nseg : ndarray
            number of segments in each section
        L : ndarray
            section lengths
        a : ndarray
            section radii
        sec_points : list
            section coordinates
        sec_e, sec_i :  ndarray
            locations of excitatory and inhibitory synapses
        b_type_e, b_type_i : list
            branch types (basal/apical/soma) for all synapses
        """
        P = self.P
        N_e, N_i, l_seg, locs_e_p, locs_i_p, branch_ids = \
        (P['N_e'], P['N_i'], P['l_seg'], P['locs_e'], P['locs_i'],
        [P['basal'], P['oblique'], P['apical']])
        A, L, a, sec_points, secs, basal, apical, trunk, axon = morphology.reconstruction(P['tree'])

        branch_ids = [basal, apical, trunk]
        nseg = np.array(len(L)*[1])
        dseg = L[secs]//(l_seg*1e-4)+1
        dseg[dseg == 1] = 2
        nseg[secs] = dseg
        # locs_e = np.array(
        #     basal + trunk + apical)  # location of excitatory synapses
        # locs_i = np.array(
        #     basal + trunk + apical)  # location of inhibitory synapses
        # locs_e = np.asarray(basal)
        # locs_i = np.asarray(basal)
        locs_e = []
        locs_i = []
        self.basal = basal
        self.apical = apical
        self.axon = axon
        if 'basal' in locs_e_p:
            locs_e = locs_e + basal
        if 'trunk' in locs_e_p:
            locs_e = locs_e + trunk
        if 'apical' in locs_e_p:
            locs_e = locs_e + apical
        if 'basal' in locs_i_p:
            locs_i = locs_i + basal
        if 'trunk' in locs_i_p:
            locs_i = locs_i + trunk
        if 'apical' in locs_i_p:
            locs_i = locs_i + apical
        locs_e = np.asarray(locs_e)
        locs_i = np.asarray(locs_i)
        sec_e = morphology.synapse_locations_rand(locs_e, N_e, nseg[locs_e], 0)
        sec_i = morphology.synapse_locations_rand(locs_i, N_i, nseg[locs_i], 0)
        b_type_e = morphology.branch_type(sec_e, branch_ids)
        b_type_i = morphology.branch_type(sec_i, branch_ids)
        return A, nseg, L, a, sec_e, sec_i, b_type_e, b_type_i

    def build_comp_matrix(self, A, nseg):
        """ Expand adjacency matric for sections into adjacency matrix for all
        segments. Uses Wye-delta transformation at branch points.

        Parameters
        ----------
        A : ndarray
            adjacency matrix for dendritic sections
        nseg : ndarray
            number of segments in each section

        Returns
        -------
        C : ndarray
            adjacency matrix for dendritic segments
        seg2sec : ndarray
            array mapping the new locations [sec, sec_location].
        """
        n_secs = A.shape[0]
        secs = [nseg[k]*[k] for k in range(n_secs)]
        secs = [s for sub in secs for s in sub]
        sec_locs = [np.arange(1/(2*nseg[k]), 1, 1/nseg[k]) for k in
                    range(n_secs)]
        sec_locs = [s for sub in sec_locs for s in sub]
        seg2sec = np.vstack((secs, sec_locs)).T
        A = np.array(A + A.T)
        branch_points = np.array([k for k in range(n_secs) if sum(A[k, k:]) >=
                                2])
        terminal_points = np.array([k for k in range(n_secs-1) if sum(A[k, k:])
                                == 0])
        M = sum(nseg)
        C = np.diag(np.ones(M - 1), 1)
        for tp in terminal_points:
            tp_new = np.where(seg2sec[:, 0] == tp)[0][-1]
            C[tp_new, tp_new + 1] = 0
        for bp in branch_points:
            bp_new = np.where(seg2sec[:, 0] == bp)[0][-1]
            C[bp_new, bp_new + 1] = 0
        for bp in branch_points:
            bp_new = np.where(seg2sec[:, 0] == bp)[0][-1]
            daughters = np.where(A[bp, :])[0]
            daughters = daughters[daughters > bp]
            d_new = [np.where(seg2sec[:, 0] == d)[0][0] for d in daughters]
            C[bp_new, d_new] = 1
            for n, d_i in enumerate(d_new):
                for d_j in d_new[n+1:]:
                    C[d_i, d_j] = 1
        C = C + C.T
        return C, seg2sec

    def g_axial(self, a_i, a_j, L_i, L_j, R_a):
        """Axial conductance from compartment j to i in unbranched section."""
        return (a_i*a_j**2)/(R_a*L_i*(L_j*a_i**2 + L_i*a_j**2))

    def g_axial_b(self, a_i, a_j, L_i, L_j, a_k, L_k, R_a):
        """Axial conductance from compartment j to i through branch point."""
        return ((a_i*a_j**2)/(L_i**2*L_j))/(R_a*(a_i**2/L_i + a_j**2/L_j +
                sum(a_k**2 / L_k)))

    def build_axial_mat(self):
        """Build and return axial conductance matrix Q."""
        R_a = self.P['R_a']
        L_s = self.L_s
        a_s = self.a_s
        C = self.C
        Q = np.zeros(C.shape)
        ind_i, ind_j = np.where(C)
        for i, j in zip(ind_i, ind_j):
                adj = np.intersect1d(np.where(C[i, :])[0], np.where(C[j, :])[0])
                a_k = a_s[adj]
                L_k = L_s[adj]
                Q[i, j] = self.g_axial_b(a_s[i], a_s[j], L_s[i], L_s[j],
                                        a_k, L_k, R_a)
        Q = Q + np.diag(-np.sum(Q, axis=1))
        return Q

    def set_weights(self, w_e, w_i):
        """Set synaptic weights."""
        self.w_e = w_e
        self.w_i = w_i

    def syn_mat(self, syns):
        """Matrix to project conductances to specified compartments with
        conversion to current in units of mS/area.

        Parameters
        ----------
        syns : ndarray
            locations of synapses

        Returns
        -------
        H : ndarray
            synapse->segment projection matrix
        syn_segs : ndarray
            segment locations of synapses

        """
        H = np.zeros((self.C.shape[0], syns.shape[1]))
        seg_soma = [0 for c in syns.T if c[0] == 0]
        seg_dend = [int(np.where(np.sum(np.abs(self.seg2sec - c), axis=1) <
                                1e-5)[0][0])
                    for c in syns.T if c[0] > 0]
        for k, s in enumerate(seg_soma + seg_dend):
            H[s, k] = 1/(1e3*self.area[s])
        syn_segs = np.array(seg_soma + seg_dend)
        return H, syn_segs

    def gauss_ops(self):
        """Returns sequence of pivots and targets for Gaussian elimination in
        solve.
        """
        targets = [np.where(self.C[:k, k])[0] for k in range(self.C.shape[0])]
        g_ops = []
        for k in range(1, self.C.shape[0]):
            for target in targets[k]:
                g_ops.append([k, target])
        g_ops = np.array(g_ops[::-1])
        return g_ops

    def fsub_ops(self):
        """Return array of non-zero elements for forward substitution in solve.
        """
        Q = self.C + np.diag(np.arange(self.C.shape[0]))
        row_reduce(Q, self.g_ops)
        Q[np.abs(Q) < 1e-10] = 0
        np.fill_diagonal(Q, 0)
        f_ops = np.vstack((np.where(Q)[0], np.where(Q)[1])).T
        return f_ops

    def init_IC(self, v_init):
        """ Inititialse voltage, gating and synaptic variables.

        Parameters
        ----------
        v_init : int
            initial voltage

        Returns
        -------
        v0 : list
            initial voltage in all compartments
        gate0 : list
            initial states of gating variables
        syn0 : list
            initial states of synapse kinetics states
        """
        N_e, N_i = self.P['N_e'], self.P['N_i']
        v0 = len(self.C)*[v_init]
        gate0 = []
        syn0 = [np.zeros((2, N_e)), np.zeros((2, N_e)), np.zeros((2, N_i))]
        return v0, gate0, syn0

    def init_IC_L5(self, v_init):
        """ Inititialse voltage, gating and synaptic variables.

        Parameters
        ----------
        v_init : int
            initial voltage

        Returns
        -------
        v0 : list
            initial voltage in all compartments
        gate0 : list
            initial states of gating variables
        syn0 : list
            initial states of synapse kinetics states
        """
        N_e, N_i = self.P['N_e'], self.P['N_i']
        v0 = len(self.C)*[v_init]
        gate0 = [nak.m_inf(v_init), nak.h_inf(v_init), nak.n_inf(v_init),
                 nak.p_inf(v_init), Ih.m_inf(v_init)]
        syn0 = [np.zeros((2, N_e)), np.zeros((2, N_e)), np.zeros((2, N_i))]
        return v0, gate0, syn0

    def set_IC(self, soln, stim, t0):
        """ Set conditions from specific time point in previous simulation.

        Parameters
        ----------
        soln :  list
            solution returned by `simulate`
        stim : list
            conductance states returned by `simulate`
        t0 : int
            time index to extract model states

        Returns
        -------
        v0 : list
            initial voltage in all compartments
        gate0 : list
            initial states of gating variables
        syn0 : list
            initial states of synapse kinetics states
        """
        ind_e, ind_i = stim[0], stim[1]
        v0 = soln[0][:, t0]
        gate0 = [] # [soln[1][:, t0], soln[2][:, t0], soln[3][:, t0], soln[4][:, t0], soln[5][:, t0]]
        syn0 = [np.zeros((2, self.P['N_e'])), np.zeros((2, self.P['N_e'])),
            np.zeros((2,self.P['N_i']))]
        syn0[0][:, ind_e] = np.vstack((stim[2][:, t0], stim[3][:, t0]))
        syn0[1][:, ind_e] = np.vstack((stim[4][:, t0], stim[5][:, t0]))
        syn0[2][:, ind_i] = np.vstack((stim[6][:, t0], stim[7][:, t0]))
        return v0, gate0, syn0

    def simulate(self, t_0, t_1, dt, IC, S_e, S_i, I_inj=0, break_flag=False):
        """Simulate instance of CModel using input sequences S_e and S_i from
        initial conditions IC. Records detailed synaptic state variables.

        Parameters
        ----------
        t_0, t_1 : int
            start and end times of simulation
        dt : float
            timestep
        IC : array_like
            initial conditions for all state variables (v0, gate0, syn0)
        S_e, S_i : array_like
            presynaptic spike times for each E and I synapse
        I_inj : int
            injected current at soma (default 0)
        break_flag : bool
            interupt simulation at time of first spike (default False)

        Returns
        -------
        t : ndarray
            time vector
        soln : list
            arrays of model states (voltage and gating variables) [v, m, h, n, p, hcn]
        stim :  list
            arrays of synaptic conductance and kinetic states and associated
            indices [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

         """
        P = self.P
        E_r, E_e, E_i, E_na, E_k, E_hcn, g_na, g_k, g_km, g_Ih, g_na_d, g_k_d, \
        g_km_d, g_Ih_d, r_na, tau_m, cm_s, cm_d, tauA, tauN, tauG, active_d, \
        active_n = (P['E_r'], P['E_e'], P['E_i'], P['E_na'], P['E_k'], P['E_hcn'],
                    P['g_na'], P['g_k'],P['g_km'], P['g_Ih'], P['g_na_d'],
                    P['g_k_d'], P['g_km_d'], P['g_Ih_d'], P['r_na'], P['tau_m'],
                    P['c_m'], P['c_m_d'], P['tauA'], P['tauN'], P['tauG'],
                    P['active_d'], P['active_n'])

        t = np.arange(t_0, t_1+dt, dt)
        if isinstance(IC, (int, float)):
            v_0, gate_0, syn_0 = self.init_IC(IC)
        else:
            v_0, gate_0, syn_0 = IC

        M = self.Q.shape[0]
        cm = np.hstack((cm_s, (M-1)*[cm_d]))

        Id = np.eye(M)
        d_inds = np.diag_indices(M)
        ind_e = np.where(S_e[:, 0] < t_1)[0]
        ind_i = np.where(S_i[:, 0] < t_1)[0]
        w_e = self.w_e[ind_e]
        w_i = self.w_i[ind_i]
        H_e = self.H_e[:, ind_e]
        H_i = self.H_i[:, ind_i]
        A_r, A_d, N_r, N_d, G_r, G_d = build_stim2(t, dt, syn_0[0][:, ind_e],
                                    syn_0[1][:, ind_e], syn_0[2][:, ind_i],
                                    S_e[ind_e], S_i[ind_i], tauA, tauN, tauG)
        I_inj *= 1/(self.area[0]*1e3)

        if active_d:
            a_inds = np.arange(M)
        else:
            a_inds = [0]
        M_active = len(a_inds)

        v = np.zeros((M, len(t)))
        n = np.zeros((M_active, len(t)))
        m = np.zeros((M_active, len(t)))
        h = np.zeros((M_active, len(t)))
        p = np.zeros((M_active, len(t)))
        hcn = np.zeros((M_active, len(t)))
        v[:, 0] = v_0
        m[:, 0], h[:, 0], n[:, 0], p[:, 0], hcn[:, 0] = gate_0
        g_na = np.hstack((g_na, (M_active-1)*[g_na_d]))
        g_k = np.hstack((g_k, (M_active-1)*[g_k_d]))
        g_km = np.hstack((g_km, (M_active-1)*[g_km_d]))
        g_Ih = np.hstack((g_Ih, (M_active-1)*[g_Ih_d]))

        J = dt*(self.Q.T*1/cm).T
        q = self.Q[d_inds]
        g_a = H_e@(w_e/(1 + r_na)*(A_d - A_r).T).T
        g_n = H_e@(w_e*r_na/(1 + r_na)*(N_d - N_r).T).T
        g_g = H_i@(w_i*(G_d - G_r).T).T

        if active_n:
            update_J = update_jacobian
            rhs = dvdt
        else:
            update_J = update_jacobian_pas
            rhs = dvdt_pas

        for k in range(1, len(t)):
            m[:, k] = m[:, k-1] + (1 - np.exp(-dt/nak.tau_m(v[a_inds, k - 1])))*(
                    nak.m_inf(v[a_inds, k - 1]) - m[:, k - 1])
            h[:, k] = h[:, k-1] + (1 - np.exp(-dt/nak.tau_h(v[a_inds, k - 1])))*(
                    nak.h_inf(v[a_inds, k - 1]) - h[:, k - 1])
            n[:, k] = n[:, k-1] + (1 - np.exp(-dt/nak.tau_n(v[a_inds, k - 1])))*(
                    nak.n_inf(v[a_inds, k - 1]) - n[:, k - 1])
            p[:, k] = p[:, k-1] + (1 - np.exp(-dt/nak.tau_p(v[a_inds, k - 1])))*(
                    nak.p_inf(v[a_inds, k - 1]) - p[:, k - 1])
            hcn[:, k] = hcn[:, k-1] + (1 - np.exp(-dt/Ih.tau_m(v[a_inds, k - 1])))*(
                    Ih.m_inf(v[a_inds, k - 1]) - hcn[:, k - 1])

            update_J(J, q, v[:, k-1], g_a[:, k], g_n[:, k], g_g[:, k],
                            E_e, tau_m, cm, dt, d_inds)
            f = rhs(v[:, k-1], g_a[:, k], g_n[:, k], g_g[:, k], self.Q, E_r,
                    E_e, E_i, tau_m, cm)
            f[0] += I_inj
            f *= dt/cm
            a = Id - J	 # note to future self: J multiplied by dt in update step
            a[a_inds, a_inds] += dt/cm[a_inds]*(g_na*m[:, k]**3*h[:, k] + g_k*n[:, k]**4 +
            g_km*p[:, k] + g_Ih*hcn[:, k])
            b = v[:, k-1] + f - J@v[:, k-1]
            b[a_inds] += dt/cm[a_inds]*(g_na*m[:, k]**3*h[:, k]*E_na + g_k*n[:, k]**4*E_k +
            g_km*p[:, k]*E_k + g_Ih*hcn[:, k]*E_hcn)
            v[:, k] = solve(a, b, self.g_ops, self.f_ops)
            if v[0, k] > 0 and break_flag:
                break
        soln = [v, m, h, n, p, hcn]
        stim = [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]
        return t, soln, stim

    def insert_biophysical_L5(self):
        P = self.P
        E_r, E_e, E_i, E_na, E_k, E_hcn, tauA, tauN, tauG, active_d, \
        active_n = (P['E_r'], P['E_e'], P['E_i'], P['E_na'], P['E_k'], P['E_hcn'], \
                    P['tauA'], P['tauN'], P['tauG'], P['active_d'], P['active_n'])
        from neuron import h
        h('forall pop_section()')
        h('forall delete_section()')
        cell = neuron_model.NModel(P, verbool = False)
        cm = np.ones(cell.totnsegs)
        gpas = 3e-5*np.ones(cell.totnsegs)
        gna = np.zeros(cell.totnsegs)     # NasT2, NaTa_t for axon
        gkv = np.zeros(cell.totnsegs)   # SKv3.1, kv for basal
        gsk = np.zeros(cell.totnsegs)   # SK_E2
        gcah = np.zeros(cell.totnsegs)   # CaHva, ca for basal
        gcal = np.zeros(cell.totnsegs)   # CaLva, it for basal
        decay = np.zeros(cell.totnsegs)   # Cadynamic decay
        gamma = np.zeros(cell.totnsegs)   # casynamic gamma
        gih = np.zeros(cell.totnsegs)   # Ih
        gim = np.zeros(cell.totnsegs)   # Im
        gnap = np.zeros(cell.totnsegs)   # Nap_Et2 (axon)
        gkt = np.zeros(cell.totnsegs)   # K_Tst (axon), kad (basal)
        gkp = np.zeros(cell.totnsegs)   # K_Pst (axon), kap (basal)
        dend = []  # index for basal dendrites
        axon = []     # index for axon
        apic = []
        soma = np.asarray([0])

        count = 0
        from neuron import h
        for k, sec in enumerate(h.allsec()):
            if k in cell.basal: # active basal dendrites
                dend.append(cell.get_idx(sec.name()))
            if k in cell.apical: # active basal dendrites
                apic.append(cell.get_idx(sec.name()))
            if k in cell.axon: # axon
                axon.append(cell.get_idx(sec.name()))
            for seg in sec:
                cm[count] = seg.cm
                gpas[count] = seg.g_pas
                if k in cell.axon: # axon
                    if hasattr(seg, 'gNaTa_tbar_NaTa_t'):
                        gna[count] = seg.gNaTa_tbar_NaTa_t
                    elif hasattr(seg, 'gNaTa_tbar_NaTa_t_2F'):
                        gna[count] = seg.gNaTa_tbar_NaTa_t_2F
                else:
                    if hasattr(seg, 'gNaTs2_tbar_NaTs2_t'):
                        gna[count] = seg.gNaTs2_tbar_NaTs2_t
                    elif hasattr(set, 'gNaTs2_tbar_NaTs2_t_2F'):
                        gna[count] = seg.gNaTs2_tbar_NaTs2_t_2F
                if hasattr(seg, 'gSKv3_1bar_SKv3_1'):
                    gkv[count] = seg.gSKv3_1bar_SKv3_1
                if hasattr(seg, 'gSK_E2bar_SK_E2'):
                    gsk[count] = seg.gSK_E2bar_SK_E2
                if hasattr(seg, 'gCa_HVAbar_Ca_HVA'):
                    gcah[count] = seg.gCa_HVAbar_Ca_HVA
                if hasattr(seg, 'gCa_LVAbar_Ca_LVA'):
                    gcal[count] = seg.gCa_LVAbar_Ca_LVA
                if hasattr(seg, 'decay_CaDynamics_E2'):
                    decay[count] = seg.decay_CaDynamics_E2
                    gamma[count] = seg.gamma_CaDynamics_E2
                if hasattr(seg, 'gIhbar_Ih'):
                    gih[count] = seg.gIhbar_Ih
                if hasattr(seg, 'gImbar_Im'):
                    gim[count] = seg.gImbar_Im
                if hasattr(seg, 'gNap_Et2bar_Nap_Et2'):
                    gnap[count] = seg.gNap_Et2bar_Nap_Et2
                elif hasattr(seg, 'gNap_Et2bar_Nap_Et2_2F'):
                    gnap[count] = seg.gNap_Et2bar_Nap_Et2_2F
                if hasattr(seg, 'gK_Tstbar_K_Tst'):
                    gkt[count] = seg.gK_Tstbar_K_Tst
                elif hasattr(seg, 'gK_Tstbar_K_Tst_2F'):
                    gkt[count] = seg.gK_Tstbar_K_Tst_2F
                if hasattr(seg, 'gK_Pstbar_K_Pst'):
                    gkp[count] = seg.gK_Pstbar_K_Pst
                elif hasattr(seg, 'gK_Pstbar_K_Pst_2F'):
                    gkp[count] = seg.gK_Pstbar_K_Pst_2F
                if active_d:
                    if k in cell.basal: # active basal dendrites
                        if hasattr(seg, 'gbar_na'):
                            gna[count] = seg.gbar_na*1e-4 #pS/um^2 -> S/cm^2
                        elif hasattr(seg, 'gbar_na_2F'):
                            gna[count] = seg.gbar_na_2F*1e-4
                        if hasattr(seg, 'gbar_kv'):
                            gkv[count] = seg.gbar_kv*1e-4
                        if hasattr(seg, 'gkabar_kad'):
                            gkt[count] = seg.gkabar_kad
                        if hasattr(seg, 'gkabar_kap'):
                            gkp[count] = seg.gkabar_kap
                        if hasattr(seg, 'gbar_ca'):
                            gcah[count] = seg.gbar_ca*1e-4
                        if hasattr(seg, 'gbar_it'):
                            gcal[count] = seg.gbar_it
                        if hasattr(seg, 'gpeak_kBK'):
                            gsk[count] = seg.gpeak_kBK/100

                count = count + 1
        if count != cell.totnsegs:
            raise ValueError("something wrong with the loop")
        h('forall pop_section()')
        h('forall delete_section()')
        axon = np.hstack(np.asarray(axon))
        apic = np.hstack(np.asarray(apic))
        dend = np.hstack(np.asarray(dend))
        return [cm, gpas, gna, gkv, gcah, gcal, gih, gim, gnap, gkt, gkp, gsk,gamma, decay, dend, axon, apic, soma]

    def simulate_L5(self, t_0, t_1, dt, IC, S_e, S_i, I_inj=0, break_flag=False):
        """Simulate instance of CModel using input sequences S_e and S_i from
        initial conditions IC. Records detailed synaptic state variables.

        Parameters
        ----------
        t_0, t_1 : int
            start and end times of simulation
        dt : float
            timestep
        IC : array_like
            initial conditions for all state variables (v0, gate0, syn0)
        S_e, S_i : array_like
            presynaptic spike times for each E and I synapse
        I_inj : int
            injected current at soma (default 0)
        break_flag : bool
            interupt simulation at time of first spike (default False)

        Returns
        -------
        t : ndarray
            time vector
        soln : list
            arrays of model states (voltage and gating variables) [v, m, h, n, p, hcn]
        stim :  list
            arrays of synaptic conductance and kinetic states and associated
            indices [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

         """
        P = self.P
        E_r, E_e, E_i, E_na, E_k, E_hcn, tauA, tauN, tauG, active_d,tau_m, \
        active_n, r_na = (P['E_r'], P['E_e'], P['E_i'], P['E_na'], P['E_k'], P['E_hcn'], \
                    P['tauA'], P['tauN'], P['tauG'], P['active_d'],P['tau_m'],  P['active_n'], P['r_na'])

        g_ion = self.g_ion
        dend = self.dend
        axon = self.axon
        apic = self.apic
        soma = self.soma
        cm = self.cm
        gamma = self.gamma
        decay = self.decay

        gcah = self.g_ion[2,:]/1e3
        gcal = self.g_ion[3,:]/1e3

        t = np.arange(t_0, t_1+dt, dt)
        if isinstance(IC, (int, float)):
            v_0, gate_0, syn_0 = self.init_IC(IC)
        else:
            v_0, gate_0, syn_0 = IC

        M = self.Q.shape[0]
        Id = np.eye(M)
        d_inds = np.diag_indices(M)
        ind_e = np.where(S_e[:, 0] < t_1)[0]
        ind_i = np.where(S_i[:, 0] < t_1)[0]
        w_e = self.w_e[ind_e]
        w_i = self.w_i[ind_i]
        H_e = self.H_e[:, ind_e]
        H_i = self.H_i[:, ind_i]
        A_r, A_d, N_r, N_d, G_r, G_d = build_stim2(t, dt, syn_0[0][:, ind_e],
                                    syn_0[1][:, ind_e], syn_0[2][:, ind_i],
                                    S_e[ind_e], S_i[ind_i], tauA, tauN, tauG)
        I_inj *= 1/(self.area[0]*1e3)


        a_inds = np.arange(M)
        M_active = len(a_inds)

        v = np.zeros((M, len(t)))
        gates = []

        n_dend = np.asarray([i for i in a_inds if i not in dend])
        v[:, 0] = v_0
        gates = np.zeros((17, M, len(t)))
        channels = [NaTs2_t(v_0[0]), NaTa_t(v_0[0]), na(v_0[0]), SKv3_1(v_0[0]), kv(v_0[0]), Ca_HVA(v_0[0]), ca(v_0[0]), Ca_LVAst(v_0[0]), it(v_0[0]), Ih(v_0[0]), Im(v_0[0]), Nap_Et2(v_0[0]), K_Tst(v_0[0]), kad(v_0[0]), K_Pst(v_0[0]), kap(v_0[0]), SK_E2(), kBK(v_0[0]), CaDynamics_E2()]
        gates[0] = np.zeros((M, len(t)))  # na_m
        gates[1] = np.zeros((M, len(t)))  # na_h
        gates[0, n_dend, 0] = channels[0].m
        gates[1, n_dend, 0] = channels[0].h
        gates[0, axon, 0] = channels[1].m
        gates[1, axon, 0] = channels[1].h
        gates[0, dend, 0] = channels[2].m
        gates[1, dend, 0] = channels[2].h
        gates[2] = np.zeros((M, len(t)))  # kv_m
        gates[2,n_dend, 0] = channels[3].m
        gates[2, dend, 0] = channels[4].m
        gates[3] = np.zeros((M, len(t)))  # cahva_m
        gates[4] = np.zeros((M, len(t)))  # cahva_h
        gates[3, n_dend, 0] = channels[5].m
        gates[4, n_dend, 0] = channels[5].h
        gates[3, dend, 0] = channels[6].m
        gates[4, dend, 0] = channels[6].h
        gates[5] = np.zeros((M, len(t)))  # calva_m
        gates[6] = np.zeros((M, len(t)))  # calva_h
        gates[5, n_dend, 0] = channels[7].m
        gates[6, n_dend, 0] = channels[7].h
        gates[5, dend, 0] = channels[8].m
        gates[6, dend, 0] = channels[8].h
        gates[7] = np.zeros((M, len(t)))  # Ih_m
        gates[7, :, 0] = channels[9].m
        gates[8] = np.zeros((M, len(t)))  # Im_m
        gates[8, :, 0] = channels[10].m
        gates[9] = np.zeros((M, len(t)))  # nap_m
        gates[10] = np.zeros((M, len(t)))  # nap_h
        gates[9, axon, 0] = channels[11].m
        gates[10, axon, 0] = channels[11].h
        gates[11] = np.zeros((M, len(t)))  # k_Tst/kad_m
        gates[12] = np.zeros((M, len(t)))  # k_Tst/kad_h
        gates[11, axon, 0] = channels[12].m
        gates[12, axon, 0] = channels[12].h
        gates[11, dend, 0] = channels[13].m
        gates[12, dend, 0] = channels[13].h
        gates[13] = np.zeros((M, len(t)))  # k_Pst/kap_m
        gates[14] = np.zeros((M, len(t)))  # k_Pst/kap_h
        gates[13, axon, 0] = channels[14].m
        gates[14, axon, 0] = channels[14].h
        gates[13, dend, 0] = channels[15].m
        gates[14, dend, 0] = channels[15].h
        gates[15] = np.zeros((M, len(t)))  #sk.bk
        gates[15, n_dend, 0] = channels[16].m
        gates[15, dend, 0] = channels[17].m
        gates[16] = np.zeros((M, len(t)))   # cai
        gates[16, :, 0] = channels[18].ca


        J = dt*(self.Q.T*1/cm).T
        q = self.Q[d_inds]
        g_a = H_e@(w_e/(1 + r_na)*(A_d - A_r).T).T
        g_n = H_e@(w_e*r_na/(1 + r_na)*(N_d - N_r).T).T
        g_g = H_i@(w_i*(G_d - G_r).T).T

        ica = np.zeros((M, len(t)))
        if active_n:
            update_J = update_jacobian
            rhs = dvdt
        else:
            update_J = update_jacobian_pas
            rhs = dvdt_pas

        for k in range(1, len(t)):
            gates[0, n_dend, k] = gates[0, n_dend, k-1] + (1 - np.exp(-dt/channels[0].mtau(v[n_dend, k - 1])))*(channels[0].minf(v[n_dend, k - 1]) - gates[0, n_dend, k-1])
            gates[1, n_dend, k] = gates[1, n_dend, k-1] + (1 - np.exp(-dt/channels[0].htau(v[n_dend, k - 1])))*(channels[0].hinf(v[n_dend, k - 1]) - gates[1,  n_dend, k-1])
            gates[0, axon, k] = gates[0, axon, k-1] + (1 - np.exp(-dt/channels[1].mtau(v[axon, k - 1])))*(channels[1].minf(v[axon, k - 1]) - gates[0, axon, k-1])
            gates[1, axon, k] = gates[1, axon, k-1] + (1 - np.exp(-dt/channels[1].htau(v[axon, k - 1])))*(channels[1].hinf(v[axon, k - 1]) - gates[1, axon, k-1])
            if active_d:
                gates[0, dend, k] = gates[0, dend, k-1] + (1 - np.exp(-dt/channels[2].mtau(v[dend, k - 1])))*(channels[2].minf(v[dend, k - 1]) - gates[0, dend, k-1])
                gates[1, dend, k] = gates[1, dend, k-1] + (1 - np.exp(-dt/channels[2].htau(v[dend, k - 1])))*(channels[2].hinf(v[dend, k - 1]) - gates[1, dend, k-1])
            gates[2,n_dend, k] = gates[2, n_dend, k-1] + (1 - np.exp(-dt/channels[3].mtau(v[n_dend, k - 1])))*(channels[3].minf(v[n_dend, k - 1]) - gates[2, n_dend, k-1])
            if active_d:
                gates[2, dend, k] = gates[2, dend, k-1] + (1 - np.exp(-dt/channels[4].mtau(v[dend, k - 1])))*(channels[4].minf(v[dend, k - 1]) - gates[2, dend, k-1])
            gates[3, n_dend, k] = gates[3, n_dend, k-1] + (1 - np.exp(-dt/channels[5].mtau(v[n_dend, k - 1])))*(channels[5].minf(v[n_dend, k - 1]) - gates[3, n_dend, k-1])
            gates[4, n_dend, k] = gates[4, n_dend, k-1] + (1 - np.exp(-dt/channels[5].htau(v[n_dend, k - 1])))*(channels[5].hinf(v[n_dend, k - 1]) - gates[4, n_dend, k-1])
            if active_d:
                gates[3, dend, k] = gates[3, dend, k-1] + (1 - np.exp(-dt/channels[6].mtau(v[dend, k - 1])))*(channels[6].minf(v[dend, k - 1]) - gates[3, dend, k-1])
                gates[4, dend, k] = gates[4, dend, k-1] + (1 - np.exp(-dt/channels[6].htau(v[dend, k - 1])))*(channels[6].hinf(v[dend, k - 1]) - gates[4, dend, k-1])
            gates[5, n_dend, k] = gates[5, n_dend, k-1] + (1 - np.exp(-dt/channels[7].mtau(v[n_dend, k - 1])))*(channels[7].minf(v[n_dend, k - 1]) - gates[5, n_dend, k-1])
            gates[6, n_dend, k] = gates[6, n_dend, k-1] + (1 - np.exp(-dt/channels[7].htau(v[n_dend, k - 1])))*(channels[7].hinf(v[n_dend, k - 1]) - gates[6, n_dend, k-1])
            if active_d:
                gates[5, dend, k] = gates[5, dend, k-1] + (1 - np.exp(-dt/channels[8].mtau(v[dend, k - 1])))*(channels[8].minf(v[dend, k - 1]) - gates[5, dend, k-1])
                gates[6, dend, k] = gates[6, dend, k-1] + (1 - np.exp(-dt/channels[8].htau(v[dend, k - 1])))*(channels[8].hinf(v[dend, k - 1]) - gates[6, dend, k-1])
            gates[7, :, k] = gates[7, :, k-1] + (1 - np.exp(-dt/channels[9].mtau(v[:, k - 1])))*(channels[9].minf(v[:, k - 1]) - gates[7, :, k-1])
            gates[8, :, k] = gates[8, :, k-1] + (1 - np.exp(-dt/channels[10].mtau(v[:, k - 1])))*(channels[10].minf(v[:, k - 1]) - gates[8, :, k-1])
            gates[9, axon, k] = gates[9, axon, k-1] + (1 - np.exp(-dt/channels[11].mtau(v[axon, k - 1])))*(channels[11].minf(v[axon, k - 1]) - gates[9, axon, k-1])
            gates[10, axon, k] = gates[10, axon, k-1] + (1 - np.exp(-dt/channels[11].htau(v[axon, k - 1])))*(channels[11].hinf(v[axon, k - 1]) - gates[10, axon, k-1])
            gates[11, axon, k] = gates[11, axon, k-1] + (1 - np.exp(-dt/channels[12].mtau(v[axon, k - 1])))*(channels[12].minf(v[axon, k - 1]) - gates[11, axon, k-1])
            gates[12, axon, k] = gates[12, axon, k-1] + (1 - np.exp(-dt/channels[12].htau(v[axon, k - 1])))*(channels[12].hinf(v[axon, k - 1]) - gates[12, axon, k-1])
            if active_d:
                gates[11, dend, k] = gates[11, dend, k-1] + (1 - np.exp(-dt/channels[13].mtau(v[dend, k - 1])))*(channels[13].minf(v[dend, k - 1]) - gates[11, dend, k-1])
                gates[12, dend, k] = gates[12, dend, k-1] + (1 - np.exp(-dt/channels[13].htau(v[dend, k - 1])))*(channels[13].hinf(v[dend, k - 1]) - gates[12, dend, k-1])
            gates[13, axon, k] = gates[13, axon, k-1] + (1 - np.exp(-dt/channels[14].mtau(v[axon, k - 1])))*(channels[14].minf(v[axon, k - 1]) - gates[13, axon, k-1])
            gates[14, axon, k] = gates[14, axon, k-1] + (1 - np.exp(-dt/channels[14].htau(v[axon, k - 1])))*(channels[14].hinf(v[axon, k - 1]) - gates[14, axon, k-1])
            if active_d:
                gates[13, dend, k] = gates[13, dend, k-1] + (1 - np.exp(-dt/channels[15].mtau(v[dend, k - 1])))*(channels[15].minf(v[dend, k - 1]) - gates[13, dend, k-1])
                gates[14, dend, k] = gates[14, dend, k-1] + (1 - np.exp(-dt/channels[15].htau(v[dend, k - 1])))*(channels[15].hinf(v[dend, k - 1]) - gates[14, dend, k-1])
            gates[15, n_dend, k] = gates[15, n_dend, k-1] + (1 - np.exp(-dt/channels[16].mtau(gates[16, n_dend, k - 1])))*(channels[16].minf(gates[16, n_dend, k - 1]) - gates[15, n_dend, k-1])
            # if active_d:
            #     gates[15, dend, k] = gates[15, dend, k-1] + (1 - np.exp(-dt/channels[17].mtau(v[dend, k - 1], gates[16, dend, k - 1])))*(channels[17].minf(v[dend, k - 1], gates[16, dend, k - 1]) - gates[15, dend, k-1])
            for kk in n_dend:
                ica[kk, k] = gcah[kk]*channels[5].g_s(gates[3, kk, k-1], gates[4, kk, k-1])*(v[kk, k-1] - channels[5].E) + gcal[kk]*channels[7].g_s(gates[5, kk, k-1], gates[6, kk, k-1])*(v[kk, k-1] - channels[7].E)
            # if active_d:
            #     for kk in dend:
            #         ica[kk, k] = gcah[kk]*channels[6].g_s(gates[3, kk, k-1], gates[4, kk, k-1])*(v[kk, k-1] - channels[6].E) + gcal[kk]*channels[8].g_s(gates[5, kk, k-1], gates[6, kk, k-1])*(v[kk, k-1] - channels[8].E)
            gates[16,:,k] = channels[18].update(ica[:,k], gates[16, :, k-1], gamma, decay, dt)


            update_J(J, q, v[:, k-1], g_a[:, k], g_n[:, k], g_g[:, k],
                            E_e, tau_m, cm, dt, d_inds)
            f = rhs(v[:, k-1], g_a[:, k], g_n[:, k], g_g[:, k], self.Q, E_r,
                    E_e, E_i, tau_m, cm)
            f[0] += I_inj
            f *= dt/cm
            a = Id - J	 # note to future self: J multiplied by dt in update step
            gates_current = gates[:,:,k].reshape((gates.shape[0], M))

            g_a_ion = np.zeros((g_ion.shape[0], M))
            g_b_ion = np.zeros((g_ion.shape[0], M))

            g_a_ion[0,n_dend] = channels[0].g_s(gates_current[0, n_dend], gates_current[1, n_dend])      # na
            g_b_ion[0,n_dend] = channels[0].g_s(gates_current[0, n_dend], gates_current[1, n_dend])*channels[0].E
            g_a_ion[0,axon] = channels[1].g_s(gates_current[0, axon], gates_current[1, axon])
            g_b_ion[0,axon] = channels[1].g_s(gates_current[0, axon], gates_current[1, axon])*channels[1].E
            if active_d:
                g_a_ion[0,dend] = channels[2].g_s(gates_current[0, dend], gates_current[1, dend])
                g_b_ion[0,dend] = channels[2].g_s(gates_current[0, dend], gates_current[1, dend])*channels[2].E

            g_a_ion[1,n_dend] = channels[3].g_s(gates_current[2, n_dend])       #kv
            g_b_ion[1,n_dend] = channels[3].g_s(gates_current[2, n_dend])*channels[3].E
            if active_d:
                g_a_ion[1,dend] = channels[4].g_s(gates_current[2, dend])
                g_b_ion[1,dend] = channels[4].g_s(gates_current[2, dend])*channels[4].E

            g_a_ion[2,n_dend] = channels[5].g_s(gates_current[3, n_dend], gates_current[4, n_dend])     #cah
            g_b_ion[2,n_dend] = channels[5].g_s(gates_current[3, n_dend], gates_current[4, n_dend])*channels[5].E
            if active_d:
                g_a_ion[2,dend] = channels[6].g_s(gates_current[3, dend], gates_current[4, dend])
                g_b_ion[2,dend] = channels[6].g_s(gates_current[3, dend], gates_current[4, dend])*channels[6].E

            g_a_ion[3,n_dend] = channels[7].g_s(gates_current[5, n_dend], gates_current[6, n_dend])     #cal
            g_b_ion[3,n_dend] = channels[7].g_s(gates_current[5, n_dend], gates_current[6, n_dend])*channels[7].E
            if active_d:
                g_a_ion[3,dend] = channels[8].g_s(gates_current[5, dend], gates_current[6, dend])
                g_b_ion[3,dend] = channels[8].g_s(gates_current[5, dend], gates_current[6, dend])*channels[8].E

            g_a_ion[4,:] = channels[9].g_s(gates_current[7, :])       #ih
            g_b_ion[4,:] = channels[9].g_s(gates_current[7, :])*channels[9].E

            g_a_ion[5,:] = channels[10].g_s(gates_current[8, :])       #im
            g_b_ion[5,:] = channels[10].g_s(gates_current[8, :])*channels[10].E

            g_a_ion[6,axon] = channels[11].g_s(gates_current[9, axon], gates_current[10, axon])  #nap
            g_b_ion[6,axon] = channels[11].g_s(gates_current[9, axon], gates_current[10, axon])*channels[11].E

            g_a_ion[7,n_dend] = channels[12].g_s(gates_current[11, n_dend], gates_current[12, n_dend])     #kad
            g_b_ion[7,n_dend] = channels[12].g_s(gates_current[11, n_dend], gates_current[12, n_dend])*channels[12].E
            if active_d:
                g_a_ion[7,dend] = channels[13].g_s(gates_current[11, dend], gates_current[12, dend])
                g_b_ion[7,dend] = channels[13].g_s(gates_current[11, dend], gates_current[12, dend])*channels[13].E

            g_a_ion[8,n_dend] = channels[14].g_s(gates_current[13, n_dend], gates_current[14, n_dend])     #kap
            g_b_ion[8,n_dend] = channels[14].g_s(gates_current[13, n_dend], gates_current[14, n_dend])*channels[14].E
            if active_d:
                g_a_ion[8,dend] = channels[15].g_s(gates_current[13, dend], gates_current[14, dend])
                g_b_ion[8,dend] = channels[15].g_s(gates_current[13, dend], gates_current[14, dend])*channels[15].E

            g_a_ion[9,n_dend] = channels[16].g_s(gates_current[15, n_dend])     #sk
            g_b_ion[9,n_dend] = channels[16].g_s(gates_current[15, n_dend])*channels[16].E
            if active_d:
                g_a_ion[9,dend] = channels[17].g_s(gates_current[15, dend])
                g_b_ion[9,dend] = channels[17].g_s(gates_current[15, dend])*channels[17].E


            a[a_inds, a_inds] += dt/cm[a_inds]*np.sum(g_ion*g_a_ion, axis = 0)
            b = v[:, k-1] + f - J@v[:, k-1]
            b[a_inds] += dt/cm[a_inds]*np.sum(g_ion*g_b_ion, axis = 0)
            v[:, k] = solve(a, b, self.g_ops, self.f_ops)
            if self.verbool:
                print('%d th time point solved'% k)
            if v[0, k] > 0 and break_flag:
                break
        soln = [v, gates, ica]
        stim = [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]

        return t, soln, stim

    def grad_w(self, soln, stim, t, dt, Z_e, Z_i, z_ind_e, z_ind_i):
        """Compute gradients associated with individual input
        spikes using solution from `simulate2'. Z_e and Z_i are expanded copies
        of the input pattern between those times (one spike per synapse; see
        `sequences.rate2temp`).

        Parameters
        ----------
        soln : list
            arrays of model states (voltage and gating variables)
            [v, m, h, n, p]. See `simulate2`.
        stim : list
            arrays of synaptic conductance states and associated indices
            [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]
        t : ndarray
            simulation time vector
        dt : float
            timestep from original forward simulation
        Z_e, Z_i : array_like
            presynaptic spike time for dummy copies of E and I synapses
        z_ind_e, z_ind_i :
            original indices of dummy synapses

        Returns
        -------
        f_e, f_i : ndarray
            dv_soma/dw for E and I synapses as a function of time
        """

        P = self.P
        E_e, E_i, tau_m, E_k, E_na, E_hcn, g_k, g_na, g_km, g_Ih, g_na_d, g_k_d, g_km_d, g_Ih_d, \
        r_na, cm_s, cm_d, tauA, tauN, tauG, active_d, active_n = (P['E_e'],
            P['E_i'], P['tau_m'], P['E_k'], P['E_na'], P['E_hcn'], P['g_k'],
            P['g_na'], P['g_km'], P['g_Ih'], P['g_na_d'], P['g_k_d'], P['g_km_d'],
            P['g_Ih_d'], P['r_na'], P['c_m'], P['c_m_d'], P['tauA'], P['tauN'],
            P['tauG'], P['active_d'], P['active_n'])
        ind_e, ind_i = stim[0], stim[1]
        M = self.Q.shape[0]
        cm = np.hstack((cm_s, (M - 1) * [cm_d]))
        if active_d:
            a_inds = np.arange(M)
        else:
            a_inds = [0]
            g_Ih = 0
        M_active = len(a_inds)
        ZA, ZN, ZG = build_stim(t, dt, Z_e, Z_i, tauA, tauN, tauG)

        N_e = len(Z_e)
        N_i = len(Z_i)
        Hz_e = np.zeros((self.H_e.shape[0], len(z_ind_e)))
        Hz_i = np.zeros((self.H_i.shape[0], len(z_ind_i)))
        Hz_e[np.where(self.H_e)[0][z_ind_e], np.arange(len(z_ind_e))] = self.H_e[
        np.where(self.H_e)[0][z_ind_e], np.where(self.H_e)[1][z_ind_e]]
        Hz_i[np.where(self.H_i)[0][z_ind_i], np.arange(len(z_ind_i))] = self.H_i[
        np.where(self.H_i)[0][z_ind_i], np.where(self.H_i)[1][z_ind_i]]
        he_inds = (np.where(Hz_e)[0], np.where(Hz_e)[1])
        hi_inds = (np.where(Hz_i)[0], np.where(Hz_i)[1] + N_e)

        g_na = np.hstack((g_na, (M_active-1)*[g_na_d]))
        g_k = np.hstack((g_k, (M_active-1)*[g_k_d]))
        g_km = np.hstack((g_km, (M_active-1)*[g_km_d]))
        g_Ih = np.hstack((g_Ih, (M_active-1)*[g_Ih_d]))

        v, m, h, n, p, hcn = soln
        GA = stim[3] - stim[2]
        GN = stim[5] - stim[4]
        GG = stim[7] - stim[6]

        w_e = self.w_e[ind_e]
        w_i = self.w_i[ind_i]
        H_e = self.H_e[:, ind_e]
        H_i = self.H_i[:, ind_i]
        dhQ = dt*(self.Q.T*1/cm).T

        Y_m = np.zeros((N_e + N_i, M_active))
        Y_h = np.zeros((N_e + N_i, M_active))
        Y_n = np.zeros((N_e + N_i, M_active))
        Y_p = np.zeros((N_e + N_i, M_active))
        Y_hcn = np.zeros((N_e + N_i, M_active))
        B = np.zeros((M, N_e + N_i))
        f_soma = B[0, :]
        f_e = np.zeros((N_e, v.shape[1]))
        f_i = np.zeros((N_i, v.shape[1]))

        a_m = nak.d_alpha_m(v[a_inds, :])*(1 - m) - nak.d_beta_m(v[a_inds, :])*m
        a_h = nak.d_alpha_h(v[a_inds, :])*(1 - h) - nak.d_beta_h(v[a_inds, :])*h
        a_n = nak.d_alpha_n(v[a_inds, :])*(1 - n) - nak.d_beta_n(v[a_inds, :])*n
        a_p = nak.d_alpha_p(v[a_inds, :])*(1 - p) - nak.d_beta_p(v[a_inds, :])*p
        a_hcn = Ih.d_alpha_m(v[a_inds, :])*(1 - hcn) - Ih.d_beta_m(v[a_inds, :])*hcn

        b_m = nak.alpha_m(v[a_inds, :]) + nak.beta_m(v[a_inds, :])
        b_h = nak.alpha_h(v[a_inds, :]) + nak.beta_h(v[a_inds, :])
        b_n = nak.alpha_n(v[a_inds, :]) + nak.beta_n(v[a_inds, :])
        b_p = nak.alpha_p(v[a_inds, :]) + nak.beta_p(v[a_inds, :])
        b_hcn = Ih.alpha_m(v[a_inds, :]) + Ih.beta_m(v[a_inds, :])

        c_m = (g_na*3*(m**2*h).T*(E_na - v[a_inds, :].T)).T
        c_h = (g_na*(m**3).T*(E_na - v[a_inds, :].T)).T
        c_n = (g_k*4*(n**3).T*(E_k - v[a_inds, :].T)).T
        c_p = (g_km*(E_k - v[a_inds, :].T)).T
        c_hcn = (g_Ih*(E_hcn - v[a_inds, :].T)).T

        if active_n:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T*sigma(v) -
                H_e@(w_e*r_na/(1 + r_na)*GN.T).T*d_sigma(v)*(E_e - v) +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            (E_e - v)*sigma(v)))*ZN
        else:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            E_e - v))*ZN

        g_s[a_inds, :] += (g_k*(n**4).T + g_na*(m**3*h).T + g_km*p.T + g_Ih*hcn.T).T
        gw_i = (Hz_i.T@(E_i - v))*ZG

        for k in range(1, v.shape[1]):
            Y_m += (a_m[:, k-1]/b_m[:, k-1]*B[a_inds, :].T - Y_m)*(
                    1 - np.exp(-dt*b_m[:, k-1]))
            Y_h += (a_h[:, k-1]/b_h[:, k-1]*B[a_inds, :].T - Y_h)*(
                    1 - np.exp(-dt*b_h[:, k-1]))
            Y_n += (a_n[:, k-1]/b_n[:, k-1]*B[a_inds, :].T - Y_n)*(
                    1 - np.exp(-dt*b_n[:, k-1]))
            Y_p += (a_p[:, k-1]/b_p[:, k-1]*B[a_inds, :].T - Y_p)*(
                    1 - np.exp(-dt*b_p[:, k-1]))
            Y_hcn += (a_hcn[:, k-1]/b_hcn[:, k-1]*B[a_inds, :].T - Y_hcn)*(
                    1 - np.exp(-dt*b_hcn[:, k-1]))
            A = np.diag(1 + dt/cm*g_s[:, k]) - dhQ
            B[he_inds] += dt/cm[self.seg_e[z_ind_e]]*gw_e[:, k]
            B[hi_inds] += dt/cm[self.seg_i[z_ind_i]]*gw_i[:, k]
            B[a_inds, :] += (dt/cm[a_inds]*(c_m[:, k]*Y_m + c_h[:, k]*Y_h +
                            c_n[:, k]*Y_n + c_p[:, k]*Y_p + c_hcn[:, k]*Y_hcn)).T
            solve_grad(A, B, self.g_ops, self.f_ops)
            f_e[:, k] = f_soma[:N_e]
            f_i[:, k] = f_soma[N_e:]
        return f_e, f_i

    def grad_w_l5(self, soln, stim, t, dt, Z_e, Z_i, z_ind_e, z_ind_i):
        """Compute gradients associated with individual input
        spikes using solution from `simulate2'. Z_e and Z_i are expanded copies
        of the input pattern between those times (one spike per synapse; see
        `sequences.rate2temp`).

        Parameters
        ----------
        soln : list
            arrays of model states (voltage and gating variables)
            [v, m, h, n, p]. See `simulate2`.
        stim : list
            arrays of synaptic conductance states and associated indices
            [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]
        t : ndarray
            simulation time vector
        dt : float
            timestep from original forward simulation
        Z_e, Z_i : array_like
            presynaptic spike time for dummy copies of E and I synapses
        z_ind_e, z_ind_i :
            original indices of dummy synapses

        Returns
        -------
        f_e, f_i : ndarray
            dv_soma/dw for E and I synapses as a function of time
        """

        P = self.P
        E_r, E_e, E_i, E_na, E_k, E_hcn, tauA, tauN, tauG, active_d,tau_m, \
        active_n, r_na = (P['E_r'], P['E_e'], P['E_i'], P['E_na'], P['E_k'], P['E_hcn'], \
                    P['tauA'], P['tauN'], P['tauG'], P['active_d'],P['tau_m'],  P['active_n'], P['r_na'])
        # [cm, gpas, gna, gkv, gcah, gcal, gih, gim, gnap, gkt, gkp, gsk,gamma, decay, dend, axon, apic, soma]  = self.insert_biophysical_L5()
        g_ion = self.g_ion[:-1,:]
        dend = self.dend
        axon = self.axon
        apic = self.apic
        soma = self.soma
        cm = self.cm

        ind_e, ind_i = stim[0], stim[1]
        M = self.Q.shape[0]
        a_inds = np.arange(M)

        M_active = len(a_inds)
        ZA, ZN, ZG = build_stim(t, dt, Z_e, Z_i, tauA, tauN, tauG)

        N_e = len(Z_e)
        N_i = len(Z_i)
        Hz_e = np.zeros((self.H_e.shape[0], len(z_ind_e)))
        Hz_i = np.zeros((self.H_i.shape[0], len(z_ind_i)))
        Hz_e[np.where(self.H_e)[0][z_ind_e], np.arange(len(z_ind_e))] = self.H_e[
        np.where(self.H_e)[0][z_ind_e], np.where(self.H_e)[1][z_ind_e]]
        Hz_i[np.where(self.H_i)[0][z_ind_i], np.arange(len(z_ind_i))] = self.H_i[
        np.where(self.H_i)[0][z_ind_i], np.where(self.H_i)[1][z_ind_i]]
        he_inds = (np.where(Hz_e)[0], np.where(Hz_e)[1])
        hi_inds = (np.where(Hz_i)[0], np.where(Hz_i)[1] + N_e)

        a_inds = np.arange(M)
        M_active = len(a_inds)

        v = np.zeros((M, len(t)))
        gates = []

        n_dend = np.asarray([i for i in a_inds if i not in dend])
        v = soln[0]
        gates = soln[1]
        v_0 = v[:,0]

        channels = [NaTs2_t(v_0[0]), NaTa_t(v_0[0]), na(v_0[0]), SKv3_1(v_0[0]), kv(v_0[0]), Ca_HVA(v_0[0]), ca(v_0[0]), Ca_LVAst(v_0[0]), it(v_0[0]), Ih(v_0[0]), Im(v_0[0]), Nap_Et2(v_0[0]), K_Tst(v_0[0]), kad(v_0[0]), K_Pst(v_0[0]), kap(v_0[0]), SK_E2(), kBK(v_0[0]), CaDynamics_E2()]

        GA = stim[3] - stim[2]
        GN = stim[5] - stim[4]
        GG = stim[7] - stim[6]

        w_e = self.w_e[ind_e]
        w_i = self.w_i[ind_i]
        H_e = self.H_e[:, ind_e]
        H_i = self.H_i[:, ind_i]
        dhQ = dt*(self.Q.T*1/cm).T


        gates_a = np.zeros(gates.shape)
        gates_a[0, n_dend, :] = channels[0].m_a(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_a[1, n_dend, :] = channels[0].h_a(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_a[0, axon, :] = channels[1].m_a(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_a[1, axon, :] = channels[1].h_a(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_a[0, dend, :] = channels[2].m_a(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_a[1, dend, :] = channels[2].h_a(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_a[2,n_dend, :] = channels[3].m_a(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_a[2, dend, :] = channels[4].m_a(v[dend,:], gates[2,dend, :])
        gates_a[3, n_dend, :] = channels[5].m_a(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_a[4, n_dend, :] = channels[5].h_a(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_a[3, dend, :] = channels[6].m_a(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_a[4, dend, :] = channels[6].h_a(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_a[5, n_dend, :] = channels[7].m_a(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_a[6, n_dend, :] = channels[7].h_a(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_a[5, dend, :] = channels[8].m_a(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_a[6, dend, :] = channels[8].h_a(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_a[7, :, :] = channels[9].m_a(v, gates[7, :, :])
        gates_a[8, :, :] = channels[10].m_a(v, gates[8, :, :])
        gates_a[9, axon, :] = channels[11].m_a(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_a[10, axon, :] = channels[11].h_a(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_a[11, axon, :] = channels[12].m_a(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_a[12, axon, :] = channels[12].h_a(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_a[11, dend, :] = channels[13].m_a(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_a[12, dend, :] = channels[13].h_a(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_a[13, axon, :] = channels[14].m_a(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_a[14, axon, :] = channels[14].h_a(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_a[13, dend, :] = channels[15].m_a(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_a[14, dend, :] = channels[15].h_a(v[dend,:], gates[13, dend, :], gates[14, dend, :])

        gates_b = np.zeros(gates.shape)
        gates_b[0, n_dend, :] = channels[0].m_b(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_b[1, n_dend, :] = channels[0].h_b(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_b[0, axon, :] = channels[1].m_b(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_b[1, axon, :] = channels[1].h_b(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_b[0, dend, :] = channels[2].m_b(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_b[1, dend, :] = channels[2].h_b(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_b[2,n_dend, :] = channels[3].m_b(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_b[2, dend, :] = channels[4].m_b(v[dend,:], gates[2,dend, :])
        gates_b[3, n_dend, :] = channels[5].m_b(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_b[4, n_dend, :] = channels[5].h_b(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_b[3, dend, :] = channels[6].m_b(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_b[4, dend, :] = channels[6].h_b(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_b[5, n_dend, :] = channels[7].m_b(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_b[6, n_dend, :] = channels[7].h_b(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_b[5, dend, :] = channels[8].m_b(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_b[6, dend, :] = channels[8].h_b(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_b[7, :, :] = channels[9].m_b(v, gates[7, :, :])
        gates_b[8, :, :] = channels[10].m_b(v, gates[8, :, :])
        gates_b[9, axon, :] = channels[11].m_b(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_b[10, axon, :] = channels[11].h_b(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_b[11, axon, :] = channels[12].m_b(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_b[12, axon, :] = channels[12].h_b(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_b[11, dend, :] = channels[13].m_b(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_b[12, dend, :] = channels[13].h_b(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_b[13, axon, :] = channels[14].m_b(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_b[14, axon, :] = channels[14].h_b(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_b[13, dend, :] = channels[15].m_b(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_b[14, dend, :] = channels[15].h_b(v[dend,:], gates[13, dend, :], gates[14, dend, :])

        gates_c = np.zeros(gates.shape)
        gates_c[0, n_dend, :] = channels[0].m_c(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_c[1, n_dend, :] = channels[0].h_c(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_c[0, axon, :] = channels[1].m_c(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_c[1, axon, :] = channels[1].h_c(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_c[0, dend, :] = channels[2].m_c(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_c[1, dend, :] = channels[2].h_c(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_c[2,n_dend, :] = channels[3].m_c(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_c[2, dend, :] = channels[4].m_c(v[dend,:], gates[2,dend, :])
        gates_c[3, n_dend, :] = channels[5].m_c(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_c[4, n_dend, :] = channels[5].h_c(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_c[3, dend, :] = channels[6].m_c(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_c[4, dend, :] = channels[6].h_c(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_c[5, n_dend, :] = channels[7].m_c(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_c[6, n_dend, :] = channels[7].h_c(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_c[5, dend, :] = channels[8].m_c(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_c[6, dend, :] = channels[8].h_c(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_c[7, :, :] = channels[9].m_c(v, gates[7, :, :])
        gates_c[8, :, :] = channels[10].m_c(v, gates[8, :, :])
        gates_c[9, axon, :] = channels[11].m_c(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_c[10, axon, :] = channels[11].h_c(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_c[11, axon, :] = channels[12].m_c(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_c[12, axon, :] = channels[12].h_c(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_c[11, dend, :] = channels[13].m_c(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_c[12, dend, :] = channels[13].h_c(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_c[13, axon, :] = channels[14].m_c(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_c[14, axon, :] = channels[14].h_c(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_c[13, dend, :] = channels[15].m_c(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_c[14, dend, :] = channels[15].h_c(v[dend,:], gates[13, dend, :], gates[14, dend, :])


        g_a_ion = np.zeros((g_ion.shape[0], M, gates.shape[2]))

        g_a_ion[0,n_dend,:] = channels[0].g_s(gates[0, n_dend,:], gates[1, n_dend,:])      # na
        g_a_ion[0,axon,:] = channels[1].g_s(gates[0, axon,:], gates[1, axon,:])
        if active_d:
            g_a_ion[0,dend,:] = channels[2].g_s(gates[0, dend,:], gates[1, dend,:])

        g_a_ion[1,n_dend,:] = channels[3].g_s(gates[2, n_dend,:])       #kv
        if active_d:
            g_a_ion[1,dend,:] = channels[4].g_s(gates[2, dend,:])

        g_a_ion[2,n_dend,:] = channels[5].g_s(gates[3, n_dend,:], gates[4, n_dend,:])     #cah
        if active_d:
            g_a_ion[2,dend,:] = channels[6].g_s(gates[3, dend,:], gates[4, dend,:])

        g_a_ion[3,n_dend,:] = channels[7].g_s(gates[5, n_dend,:], gates[6, n_dend,:])     #cal
        if active_d:
            g_a_ion[3,dend,:] = channels[8].g_s(gates[5, dend,:], gates[6, dend,:])

        g_a_ion[4,:,:] = channels[9].g_s(gates[7, :,:])       #ih

        g_a_ion[5,:,:] = channels[10].g_s(gates[8, :,:])       #im

        g_a_ion[6,axon,:] = channels[11].g_s(gates[9, axon,:], gates[10, axon,:])  #nap

        g_a_ion[7,n_dend,:] = channels[12].g_s(gates[11, n_dend,:], gates[12, n_dend,:])     #kad
        if active_d:
            g_a_ion[7,dend,:] = channels[13].g_s(gates[11, dend,:], gates[12, dend,:])

        g_a_ion[8,n_dend,:] = channels[14].g_s(gates[13, n_dend,:], gates[14, n_dend,:])     #kap
        if active_d:
            g_a_ion[8,dend,:] = channels[15].g_s(gates[13, dend,:], gates[14, dend,:])



        if active_n:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T*sigma(v) -
                H_e@(w_e*r_na/(1 + r_na)*GN.T).T*d_sigma(v)*(E_e - v) +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            (E_e - v)*sigma(v)))*ZN
        else:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            E_e - v))*ZN

        for k in range(g_a_ion.shape[0]):
            g_s += (g_ion[k]*g_a_ion[k,:,:].T).T

        gw_i = (Hz_i.T@(E_i - v))*ZG

        gates_Y = np.zeros((gates.shape[0], N_e + N_i, M))
        B = np.zeros((M, N_e + N_i))
        f_soma = B[0, :]
        f_e = np.zeros((N_e, v.shape[1]))
        f_i = np.zeros((N_i, v.shape[1]))
        gates_b1 = gates_b
        gates_b1[np.where(gates_b1==0)] = np.inf
        for k in range(1, v.shape[1]):
            # Y_m += (a_m[:, k-1]/b_m[:, k-1]*B[a_inds, :].T - Y_m)*(
            #         1 - np.exp(-dt*b_m[:, k-1]))
            # Y_h += (a_h[:, k-1]/b_h[:, k-1]*B[a_inds, :].T - Y_h)*(
            #         1 - np.exp(-dt*b_h[:, k-1]))
            # Y_n += (a_n[:, k-1]/b_n[:, k-1]*B[a_inds, :].T - Y_n)*(
            #         1 - np.exp(-dt*b_n[:, k-1]))
            # Y_p += (a_p[:, k-1]/b_p[:, k-1]*B[a_inds, :].T - Y_p)*(
            #         1 - np.exp(-dt*b_p[:, k-1]))
            # Y_hcn += (a_hcn[:, k-1]/b_hcn[:, k-1]*B[a_inds, :].T - Y_hcn)*(
            #         1 - np.exp(-dt*b_hcn[:, k-1]))
            for kk in range(gates.shape[0]):
                gates_Y[kk] += (gates_a[kk, :, k-1]/gates_b[kk, :, k-1]*B.T - gates_Y[kk])*(1 - np.exp(-dt*gates_b[kk, :, k-1]))
            A = np.diag(1 + dt/cm*g_s[:, k]) - dhQ
            B[he_inds] += dt/cm[self.seg_e[z_ind_e]]*gw_e[:, k]
            B[hi_inds] += dt/cm[self.seg_i[z_ind_i]]*gw_i[:, k]
            Y_temp = np.zeros(B.shape)
            for kk in range(gates.shape[0]):
                Y_temp += (gates_c[kk,:,k]*gates_Y[kk,:,:]).T
            B[a_inds, :] += (dt/cm[a_inds]*Y_temp.T).T
            solve_grad(A, B, self.g_ops, self.f_ops)
            f_e[:, k] = f_soma[:N_e]
            f_i[:, k] = f_soma[N_e:]
            if self.verbool:
                print('%d th time point solved' % k)
        return f_e, f_i

    def grad_w_l5_channels(self, soln, stim, t, dt, Z_e, Z_i, z_ind_e, z_ind_i):
        """Compute gradients associated with individual input
        spikes using solution from `simulate2'. Z_e and Z_i are expanded copies
        of the input pattern between those times (one spike per synapse; see
        `sequences.rate2temp`).

        Parameters
        ----------
        soln : list
            arrays of model states (voltage and gating variables)
            [v, m, h, n, p]. See `simulate2`.
        stim : list
            arrays of synaptic conductance states and associated indices
            [ind_e, ind_i, A_r, A_d, N_r, N_d, G_r, G_d]
        t : ndarray
            simulation time vector
        dt : float
            timestep from original forward simulation
        Z_e, Z_i : array_like
            presynaptic spike time for dummy copies of E and I synapses
        z_ind_e, z_ind_i :
            original indices of dummy synapses

        Returns
        -------
        f_e, f_i : ndarray
            dv_soma/dw for E and I synapses as a function of time
        gates_Y_e, gates_Y_i: dgate/dw for E and I synapses as a function of time
        """

        P = self.P
        E_r, E_e, E_i, E_na, E_k, E_hcn, tauA, tauN, tauG, active_d,tau_m, \
        active_n, r_na = (P['E_r'], P['E_e'], P['E_i'], P['E_na'], P['E_k'], P['E_hcn'], \
                    P['tauA'], P['tauN'], P['tauG'], P['active_d'],P['tau_m'],  P['active_n'], P['r_na'])
        # [cm, gpas, gna, gkv, gcah, gcal, gih, gim, gnap, gkt, gkp, gsk,gamma, decay, dend, axon, apic, soma]  = self.insert_biophysical_L5()
        g_ion = self.g_ion[:-1,:]
        dend = self.dend
        axon = self.axon
        apic = self.apic
        soma = self.soma
        cm = self.cm

        ind_e, ind_i = stim[0], stim[1]
        M = self.Q.shape[0]
        a_inds = np.arange(M)

        M_active = len(a_inds)
        ZA, ZN, ZG = build_stim(t, dt, Z_e, Z_i, tauA, tauN, tauG)

        N_e = len(Z_e)
        N_i = len(Z_i)
        Hz_e = np.zeros((self.H_e.shape[0], len(z_ind_e)))
        Hz_i = np.zeros((self.H_i.shape[0], len(z_ind_i)))
        Hz_e[np.where(self.H_e)[0][z_ind_e], np.arange(len(z_ind_e))] = self.H_e[
        np.where(self.H_e)[0][z_ind_e], np.where(self.H_e)[1][z_ind_e]]
        Hz_i[np.where(self.H_i)[0][z_ind_i], np.arange(len(z_ind_i))] = self.H_i[
        np.where(self.H_i)[0][z_ind_i], np.where(self.H_i)[1][z_ind_i]]
        he_inds = (np.where(Hz_e)[0], np.where(Hz_e)[1])
        hi_inds = (np.where(Hz_i)[0], np.where(Hz_i)[1] + N_e)

        a_inds = np.arange(M)
        M_active = len(a_inds)

        v = np.zeros((M, len(t)))
        gates = []

        n_dend = np.asarray([i for i in a_inds if i not in dend])
        v = soln[0]
        gates = soln[1]
        v_0 = v[:,0]

        channels = [NaTs2_t(v_0[0]), NaTa_t(v_0[0]), na(v_0[0]), SKv3_1(v_0[0]), kv(v_0[0]), Ca_HVA(v_0[0]), ca(v_0[0]), Ca_LVAst(v_0[0]), it(v_0[0]), Ih(v_0[0]), Im(v_0[0]), Nap_Et2(v_0[0]), K_Tst(v_0[0]), kad(v_0[0]), K_Pst(v_0[0]), kap(v_0[0]), SK_E2(), kBK(v_0[0]), CaDynamics_E2()]

        GA = stim[3] - stim[2]
        GN = stim[5] - stim[4]
        GG = stim[7] - stim[6]

        w_e = self.w_e[ind_e]
        w_i = self.w_i[ind_i]
        H_e = self.H_e[:, ind_e]
        H_i = self.H_i[:, ind_i]
        dhQ = dt*(self.Q.T*1/cm).T


        gates_a = np.zeros(gates.shape)
        gates_a[0, n_dend, :] = channels[0].m_a(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_a[1, n_dend, :] = channels[0].h_a(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_a[0, axon, :] = channels[1].m_a(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_a[1, axon, :] = channels[1].h_a(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_a[0, dend, :] = channels[2].m_a(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_a[1, dend, :] = channels[2].h_a(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_a[2,n_dend, :] = channels[3].m_a(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_a[2, dend, :] = channels[4].m_a(v[dend,:], gates[2,dend, :])
        gates_a[3, n_dend, :] = channels[5].m_a(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_a[4, n_dend, :] = channels[5].h_a(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_a[3, dend, :] = channels[6].m_a(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_a[4, dend, :] = channels[6].h_a(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_a[5, n_dend, :] = channels[7].m_a(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_a[6, n_dend, :] = channels[7].h_a(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_a[5, dend, :] = channels[8].m_a(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_a[6, dend, :] = channels[8].h_a(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_a[7, :, :] = channels[9].m_a(v, gates[7, :, :])
        gates_a[8, :, :] = channels[10].m_a(v, gates[8, :, :])
        gates_a[9, axon, :] = channels[11].m_a(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_a[10, axon, :] = channels[11].h_a(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_a[11, axon, :] = channels[12].m_a(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_a[12, axon, :] = channels[12].h_a(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_a[11, dend, :] = channels[13].m_a(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_a[12, dend, :] = channels[13].h_a(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_a[13, axon, :] = channels[14].m_a(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_a[14, axon, :] = channels[14].h_a(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_a[13, dend, :] = channels[15].m_a(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_a[14, dend, :] = channels[15].h_a(v[dend,:], gates[13, dend, :], gates[14, dend, :])

        gates_b = np.zeros(gates.shape)
        gates_b[0, n_dend, :] = channels[0].m_b(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_b[1, n_dend, :] = channels[0].h_b(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_b[0, axon, :] = channels[1].m_b(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_b[1, axon, :] = channels[1].h_b(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_b[0, dend, :] = channels[2].m_b(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_b[1, dend, :] = channels[2].h_b(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_b[2,n_dend, :] = channels[3].m_b(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_b[2, dend, :] = channels[4].m_b(v[dend,:], gates[2,dend, :])
        gates_b[3, n_dend, :] = channels[5].m_b(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_b[4, n_dend, :] = channels[5].h_b(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_b[3, dend, :] = channels[6].m_b(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_b[4, dend, :] = channels[6].h_b(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_b[5, n_dend, :] = channels[7].m_b(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_b[6, n_dend, :] = channels[7].h_b(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_b[5, dend, :] = channels[8].m_b(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_b[6, dend, :] = channels[8].h_b(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_b[7, :, :] = channels[9].m_b(v, gates[7, :, :])
        gates_b[8, :, :] = channels[10].m_b(v, gates[8, :, :])
        gates_b[9, axon, :] = channels[11].m_b(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_b[10, axon, :] = channels[11].h_b(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_b[11, axon, :] = channels[12].m_b(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_b[12, axon, :] = channels[12].h_b(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_b[11, dend, :] = channels[13].m_b(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_b[12, dend, :] = channels[13].h_b(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_b[13, axon, :] = channels[14].m_b(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_b[14, axon, :] = channels[14].h_b(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_b[13, dend, :] = channels[15].m_b(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_b[14, dend, :] = channels[15].h_b(v[dend,:], gates[13, dend, :], gates[14, dend, :])

        gates_c = np.zeros(gates.shape)
        gates_c[0, n_dend, :] = channels[0].m_c(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_c[1, n_dend, :] = channels[0].h_c(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_c[0, axon, :] = channels[1].m_c(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_c[1, axon, :] = channels[1].h_c(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_c[0, dend, :] = channels[2].m_c(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_c[1, dend, :] = channels[2].h_c(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_c[2,n_dend, :] = channels[3].m_c(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_c[2, dend, :] = channels[4].m_c(v[dend,:], gates[2,dend, :])
        gates_c[3, n_dend, :] = channels[5].m_c(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_c[4, n_dend, :] = channels[5].h_c(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_c[3, dend, :] = channels[6].m_c(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_c[4, dend, :] = channels[6].h_c(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_c[5, n_dend, :] = channels[7].m_c(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_c[6, n_dend, :] = channels[7].h_c(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_c[5, dend, :] = channels[8].m_c(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_c[6, dend, :] = channels[8].h_c(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_c[7, :, :] = channels[9].m_c(v, gates[7, :, :])
        gates_c[8, :, :] = channels[10].m_c(v, gates[8, :, :])
        gates_c[9, axon, :] = channels[11].m_c(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_c[10, axon, :] = channels[11].h_c(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_c[11, axon, :] = channels[12].m_c(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_c[12, axon, :] = channels[12].h_c(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_c[11, dend, :] = channels[13].m_c(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_c[12, dend, :] = channels[13].h_c(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_c[13, axon, :] = channels[14].m_c(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_c[14, axon, :] = channels[14].h_c(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_c[13, dend, :] = channels[15].m_c(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_c[14, dend, :] = channels[15].h_c(v[dend,:], gates[13, dend, :], gates[14, dend, :])

        gates_d = np.zeros(gates.shape)
        gates_d[0, n_dend, :] = channels[0].m_d(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_d[1, n_dend, :] = channels[0].h_d(v[n_dend, :], gates[0, n_dend, :], gates[1, n_dend, :])
        gates_d[0, axon, :] = channels[1].m_d(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        gates_d[1, axon, :] = channels[1].h_d(v[axon, :], gates[0, axon, :], gates[1, axon, :])
        if active_d:
            gates_d[0, dend, :] = channels[2].m_d(v[dend, :], gates[0, dend, :], gates[1, dend, :])
            gates_d[1, dend, :] = channels[2].h_d(v[dend, :], gates[0, dend, :], gates[1, dend, :])
        gates_d[2,n_dend, :] = channels[3].m_d(v[n_dend,:], gates[2,n_dend, :])
        if active_d:
            gates_d[2, dend, :] = channels[4].m_d(v[dend,:], gates[2,dend, :])
        gates_d[3, n_dend, :] = channels[5].m_d(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        gates_d[4, n_dend, :] = channels[5].h_d(v[n_dend,:], gates[3, n_dend, :], gates[4, n_dend, :])
        if active_d:
            gates_d[3, dend, :] = channels[6].m_d(v[dend,:], gates[3, dend, :], gates[4, dend, :])
            gates_d[4, dend, :] = channels[6].h_d(v[dend,:], gates[3, dend, :], gates[4, dend, :])
        gates_d[5, n_dend, :] = channels[7].m_d(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        gates_d[6, n_dend, :] = channels[7].h_d(v[n_dend,:], gates[5, n_dend, :], gates[6, n_dend, :])
        if active_d:
            gates_d[5, dend, :] = channels[8].m_d(v[dend,:], gates[5, dend, :], gates[6, dend, :])
            gates_d[6, dend, :] = channels[8].h_d(v[dend,:], gates[5, dend, :], gates[6, dend, :])
        gates_d[7, :, :] = channels[9].m_d(v, gates[7, :, :])
        gates_d[8, :, :] = channels[10].m_d(v, gates[8, :, :])
        gates_d[9, axon, :] = channels[11].m_d(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_d[10, axon, :] = channels[11].h_d(v[axon,:], gates[9, axon, :], gates[10, axon, :])
        gates_d[11, axon, :] = channels[12].m_d(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        gates_d[12, axon, :] = channels[12].h_d(v[axon,:], gates[11, axon, :], gates[12, axon, :])
        if active_d:
            gates_d[11, dend, :] = channels[13].m_d(v[dend,:], gates[11, dend, :], gates[12, dend, :])
            gates_d[12, dend, :] = channels[13].h_d(v[dend,:], gates[11, dend, :], gates[12, dend, :])
        gates_d[13, axon, :] = channels[14].m_d(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        gates_d[14, axon, :] = channels[14].h_d(v[axon,:], gates[13, axon, :], gates[14, axon, :])
        if active_d:
            gates_d[13, dend, :] = channels[15].m_d(v[dend,:], gates[13, dend, :], gates[14, dend, :])
            gates_d[14, dend, :] = channels[15].h_d(v[dend,:], gates[13, dend, :], gates[14, dend, :])


        g_a_ion = np.zeros((g_ion.shape[0], M, gates.shape[2]))

        g_a_ion[0,n_dend,:] = channels[0].g_s(gates[0, n_dend,:], gates[1, n_dend,:])      # na
        g_a_ion[0,axon,:] = channels[1].g_s(gates[0, axon,:], gates[1, axon,:])
        if active_d:
            g_a_ion[0,dend,:] = channels[2].g_s(gates[0, dend,:], gates[1, dend,:])

        g_a_ion[1,n_dend,:] = channels[3].g_s(gates[2, n_dend,:])       #kv
        if active_d:
            g_a_ion[1,dend,:] = channels[4].g_s(gates[2, dend,:])

        g_a_ion[2,n_dend,:] = channels[5].g_s(gates[3, n_dend,:], gates[4, n_dend,:])     #cah
        if active_d:
            g_a_ion[2,dend,:] = channels[6].g_s(gates[3, dend,:], gates[4, dend,:])

        g_a_ion[3,n_dend,:] = channels[7].g_s(gates[5, n_dend,:], gates[6, n_dend,:])     #cal
        if active_d:
            g_a_ion[3,dend,:] = channels[8].g_s(gates[5, dend,:], gates[6, dend,:])

        g_a_ion[4,:,:] = channels[9].g_s(gates[7, :,:])       #ih

        g_a_ion[5,:,:] = channels[10].g_s(gates[8, :,:])       #im

        g_a_ion[6,axon,:] = channels[11].g_s(gates[9, axon,:], gates[10, axon,:])  #nap

        g_a_ion[7,n_dend,:] = channels[12].g_s(gates[11, n_dend,:], gates[12, n_dend,:])     #kad
        if active_d:
            g_a_ion[7,dend,:] = channels[13].g_s(gates[11, dend,:], gates[12, dend,:])

        g_a_ion[8,n_dend,:] = channels[14].g_s(gates[13, n_dend,:], gates[14, n_dend,:])     #kap
        if active_d:
            g_a_ion[8,dend,:] = channels[15].g_s(gates[13, dend,:], gates[14, dend,:])



        if active_n:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T*sigma(v) -
                H_e@(w_e*r_na/(1 + r_na)*GN.T).T*d_sigma(v)*(E_e - v) +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            (E_e - v)*sigma(v)))*ZN
        else:
            g_s = (H_e@(w_e/(1 + r_na)*GA.T).T + H_e@(w_e*r_na/(1 + r_na)*GN.T).T +
                H_i@(w_i*GG.T).T)
            g_s = (g_s.T + cm/tau_m).T
            gw_e = 1/(1 + r_na)*(Hz_e.T@(E_e - v))*ZA + r_na/(1 + r_na)*(Hz_e.T@(
            E_e - v))*ZN

        for k in range(g_a_ion.shape[0]):
            g_s += (g_ion[k]*g_a_ion[k,:,:].T).T

        gw_i = (Hz_i.T@(E_i - v))*ZG

        gates_Y = np.zeros((gates.shape[0], N_e + N_i, M))
        c = np.zeros((g_ion.shape[0], M,N_e + N_i, v.shape[1]))
        B = np.zeros((M, N_e + N_i))
        f_soma = B[0, :]
        f_e = np.zeros((M,N_e,v.shape[1]))
        f_i = np.zeros((M,N_i,v.shape[1]))
        gates_b1 = gates_b
        gates_b1[np.where(gates_b1==0)] = np.inf
        for k in range(1, v.shape[1]):
            # Y_m += (a_m[:, k-1]/b_m[:, k-1]*B[a_inds, :].T - Y_m)*(
            #         1 - np.exp(-dt*b_m[:, k-1]))
            # Y_h += (a_h[:, k-1]/b_h[:, k-1]*B[a_inds, :].T - Y_h)*(
            #         1 - np.exp(-dt*b_h[:, k-1]))
            # Y_n += (a_n[:, k-1]/b_n[:, k-1]*B[a_inds, :].T - Y_n)*(
            #         1 - np.exp(-dt*b_n[:, k-1]))
            # Y_p += (a_p[:, k-1]/b_p[:, k-1]*B[a_inds, :].T - Y_p)*(
            #         1 - np.exp(-dt*b_p[:, k-1]))
            # Y_hcn += (a_hcn[:, k-1]/b_hcn[:, k-1]*B[a_inds, :].T - Y_hcn)*(
            #         1 - np.exp(-dt*b_hcn[:, k-1]))
            for kk in range(gates.shape[0]):
                gates_Y[kk] += (gates_a[kk, :, k-1]/gates_b[kk, :, k-1]*B.T - gates_Y[kk])*(1 - np.exp(-dt*gates_b[kk, :, k-1]))
            A = np.diag(1 + dt/cm*g_s[:, k]) - dhQ
            B[he_inds] += dt/cm[self.seg_e[z_ind_e]]*gw_e[:, k]
            B[hi_inds] += dt/cm[self.seg_i[z_ind_i]]*gw_i[:, k]
            Y_temp = np.zeros(B.shape)
            for kk in range(gates.shape[0]):
                Y_temp += (gates_c[kk,:,k]*gates_Y[kk,:,:]).T
            B[a_inds, :] += (dt/cm[a_inds]*Y_temp.T).T
            solve_grad(A, B, self.g_ops, self.f_ops)
            f_e[:,:, k] = B[:,:N_e]
            f_i[:,:, k] = B[:,N_e:]
            c1_temp = np.zeros((gates.shape[0], M, N_e + N_i))
            for kk in range(gates.shape[0]):
                c1_temp[kk,:,:] = (gates_d[kk,:,k]*gates_Y[kk,:,:]).T
            c_temp = np.zeros((g_ion.shape[0], M, N_e + N_i))
            c_temp[0,:,:] = c1_temp[0, :,:] + c1_temp[1, :,:]      # na
            c_temp[1,:,:] = c1_temp[2, :,:]       #kv
            c_temp[2,:,:] = c1_temp[3, :,:] +  c1_temp[4, :,:]     #cah
            c_temp[3,:,:] = c1_temp[5, :,:] + c1_temp[6, :,:]     #cal
            c_temp[4,:,:] = c1_temp[7, :,:]       #ih
            c_temp[5,:,:] = c1_temp[8, :,:]       #im
            c_temp[6,:,:] = c1_temp[9, :,:] + c1_temp[10, :,:]  #nap
            c_temp[7,:,:] = c1_temp[11, :,:]+ c1_temp[12, :,:]    #kad
            c_temp[8,:,:] = c1_temp[13, :,:] +  c1_temp[14, :,:]     #kap
            c[:,:,:,k] = c_temp

            if self.verbool:
                print('%d th time point solved' % k)
        c_e = c[:, :,:N_e, :]
        c_i = c[:, :,:N_i, :]
        return f_e, f_i, c_e, c_i

@nb.jit(nopython=True, cache=True)
def kernel(t, tau_r, tau_d):
    """Returns value of double-exponential synaptic kernel.

    Parameters
    ---------
    t : float
        time
    tau_r, tau_d : float
        rise and decay time constants
    """
    t_peak = tau_r*tau_d/(tau_d - tau_r)*np.log(tau_d/tau_r)
    Z = np.exp(-t_peak/tau_d) - np.exp(-t_peak/tau_r)
    return 1/Z*(np.exp(-t/tau_d) - np.exp(-t/tau_r))


@nb.jit(nopython=True, cache=True)
def g(t, S, tau_r, tau_d):
    """Returns vector of synaptic conductances for sequence S at time t.

    Parameters
    ----------
    t : float
        time
    S : ndarray
        presynaptic spike times
    tau_r, tau_d : float
        rise and decay time constants
    """
    s_vec = (t - S)
    for i in range(s_vec.shape[0]):
            for j in range(s_vec.shape[1]):
                if ~(s_vec[i, j] > 0):
                    s_vec[i, j] = 0
    return np.sum(kernel(s_vec, tau_r, tau_d), axis=1)


def sigma(v):
    """NMDA voltage nonlinearity.

    Parameters
    ----------
    v : array_like
        voltage (mV)
    """
    return 1/(1 + 1/3.75*np.exp(-0.062*v))


def d_sigma(v):
    """Derivative of NMDA nonlinearity with respect to v.

    Parameters
    ----------
    v : array_like
        voltage (mV)
    """
    return 0.062*sigma(v)*(1 - sigma(v))


@nb.jit(nopython=True, cache=True)
def build_stim(t, dt, S_e, S_i, tauA, tauN, tauG):
    """AMPA, NMDA and GABA conductance time series.

    Parameters
    ----------
    t : ndarray
        time vector
    dt : float
        timestep
    S_e, S_i : ndarray
        presynaptic spike times for each E and I synapse
    tauA, tauN, tauG : list
        rise and decay time constants for AMPA, NMDA, GABA receptors

    Returns
    -------
    GA, GN, GG : ndarray
        conductance states as a function of time for AMPA, NMDA, GABA receptors
    """
    GA = build_G(t, dt, S_e, tauA[0], tauA[1])
    GN = build_G(t, dt, S_e, tauN[0], tauN[1])
    GG = build_G(t, dt, S_i, tauG[0], tauG[1])
    return GA, GN, GG


@nb.jit(nopython=True, cache=True)
def build_G(t, dt, S, tau_r, tau_d):
    """Build synaptic conductance time series using two-state scheme

    Parameters
    ----------
    t : ndarray
        time vector
    dt : float
        timestep
    S : ndarray
        presynaptic spike times for set of synapses
    tau_r, tau_d : float
        rise and decay time constants

    Returns
    -------
    G : ndarray
        conductance states as a function of time
    """
    G = np.zeros((len(S), len(t)))
    r = np.zeros(len(S))
    d = np.zeros(len(S))
    alpha_r = np.exp(-dt/tau_r)
    alpha_d = np.exp(-dt/tau_d)
    t_peak = tau_r*tau_d/(tau_d - tau_r)*np.log(tau_d/tau_r)
    Z = np.exp(-t_peak/tau_d) - np.exp(-t_peak/tau_r)
    for k, t_k in enumerate(t):
        r *= alpha_r
        d *= alpha_d
        dt_k = S - t_k
        ind = np.where((dt_k > 0) & (dt_k < dt))
        for j, i in enumerate(ind[0]):
            r[i] += 1/Z*np.exp(-dt_k[i, ind[1][j]]/tau_r)
            d[i] += 1/Z*np.exp(-dt_k[i, ind[1][j]]/tau_d)
        G[:, k] = d - r
    return G


@nb.jit(nopython=True, cache=True)
def build_stim2(t, dt, a_init, n_init, g_init, S_e, S_i, tauA, tauN, tauG):
    """AMPA, NMDA and GABA conductance time series with detailed kinetic states.

    Parameters
    ----------
    t : ndarray
        time vector
    dt : float
        timestep
    a_init, n_init, g_init : list
        initial conditions for AMPA, NMDA, GABA receptors
    S_e, S_i : ndarray
        presynaptic spike times for each E and I synapse
    tauA, tauN, tauG : list
        rise and decay time constants for AMPA, NMDA, GABA receptors

    Returns
    -------
    A_r, A_d, N_r, N_d, G_r, G_d : ndarray
        kinetic states as a function of time for AMPA, NMDA, GABA receptors
    """
    A_r, A_d = build_G2(t, dt, S_e, tauA[0], tauA[1], a_init[0], a_init[1])
    N_r, N_d = build_G2(t, dt, S_e, tauN[0], tauN[1], n_init[0], n_init[1])
    G_r, G_d = build_G2(t, dt, S_i, tauG[0], tauG[1], g_init[0], g_init[1])
    return A_r, A_d, N_r, N_d, G_r, G_d


@nb.jit(nopython=True, cache=True)
def build_G2(t, dt, S, tau_r, tau_d, r_init, d_init):
    """Build synaptic conductance time series using two-state scheme, and return
    kinetic states

    Parameters
    ----------
    t : ndarray
        time vector
    dt : float
        timestep
    S : ndarray
        presynaptic spike times for set of synapses
    tau_r, tau_d : float
        rise and decay time constants
    r_init, d_init : ndarray
        kinetics state initial conditions

    Returns
    -------
    R, D : ndarray
        kinetic states as a function of time
    """
    R = np.zeros((len(S), len(t)))
    D = np.zeros((len(S), len(t)))
    r = r_init
    d = d_init
    R[:, 0] = r_init
    D[:, 0] = d_init
    alpha_r = np.exp(-dt/tau_r)
    alpha_d = np.exp(-dt/tau_d)
    t_peak = tau_r*tau_d/(tau_d - tau_r)*np.log(tau_d/tau_r)
    Z = np.exp(-t_peak/tau_d) - np.exp(-t_peak/tau_r)
    for k, t_k in enumerate(t[1:]):
        r *= alpha_r
        d *= alpha_d
        dt_k = S - t_k
        ind = np.where((dt_k > 0) & (dt_k < dt))
        for j, i in enumerate(ind[0]):
            r[i] += 1/Z*np.exp(-dt_k[i, ind[1][j]]/tau_r)
            d[i] += 1/Z*np.exp(-dt_k[i, ind[1][j]]/tau_d)

        R[:, k+1] = r
        D[:, k+1] = d
    return R, D


def dvdt(v, g_a, g_n, g_g, Q, E_r, E_e, E_i, tau_m, cm):
    """ Returns right-hand side of ODE system.

    Parameters
    ----------
    v : ndarray
        voltage in all compartments
    g_a, g_n, gg : ndarray
        AMPA, NMDA GABA conductance states
    Q : ndarray
        axial conductance matrix
    E_r, E_e, E_i : int
        resting, E, and I reversal potentials
    tau_m : float
        membrane time constant
    I_inj : float
        injected current at soma
    """
    return (g_a*(E_e - v) + g_n*sigma(v)*(E_e - v) + g_g*(E_i - v) + cm/tau_m*(E_r - v) + Q@v)


def update_jacobian(J, q, v, g_a, g_n, g_g, E_e, tau_m, cm, dt, d_inds):
    """Update ODE Jacobian matrix.

    Parameters
    ----------
    J : ndarry
        Jacobian matrix
    q : ndarray
        diagonal of axial conductance matrix
    v : ndarray
        voltage in all compartments
    g_a, g_n, gg : ndarray
        AMPA, NMDA GABA conductance states
    E_e : int
        resting potential
    tau_m : float
        membrane time constant
    dt : float
        time step
    d_inds : tuple
        diagonal indices of J
    """
    J[d_inds] = dt/cm*(-g_a - g_n*sigma(v) + g_n*d_sigma(v)*(E_e - v) - g_g - cm/tau_m + q)


def dvdt_pas(v, g_a, g_n, g_g, Q, E_r, E_e, E_i, tau_m, cm):
    """ Returns right-hand side of ODE system with Ohmic NMDA receptors

    Parameters
    ----------
    v : ndarray
        voltage in all compartments
    g_a, g_n, gg : ndarray
        AMPA, NMDA GABA conductance states
    Q : ndarray
        axial conductance matrix
    E_r, E_e, E_i : int
        resting, E, and I reversal potentials
    tau_m : float
        membrane time constant
    I_inj : float
        injected current at soma
    """
    return (g_a*(E_e - v) + g_n*(E_e - v) + g_g*(E_i - v)) + (cm/tau_m*(E_r - v) +
            Q@v)


def update_jacobian_pas(J, q, v, g_a, g_n, g_g, E_e, tau_m, cm, dt, d_inds):
    """Update ODE Jacobian matrix with Ohmic NMDa receptors

    Parameters
    ----------
    J : ndarry
        Jacobian matrix
    q : ndarray
        diagonal of axial conductance matrix
    v : ndarray
        voltage in all compartments (unused)
    g_a, g_n, gg : ndarray
        AMPA, NMDA GABA conductance states
    E_e : int
        resting potential (unused)
    tau_m : float
        membrane time constant
    dt : float
        time step
    d_inds : tuple
        diagonal indices of J
    """
    J[d_inds] = dt/cm*(-g_a - g_n - g_g - cm/tau_m + q)


def solve(Q, b, g_ops, f_ops):
    """Solve linear system of equations Qx=b with Gaussian elimination
    (using v[0] as soma requires clearing upper triangle first).

    Parameters
    ----------
    Q : ndarray
        coefficient matrix
    b : ndarray
        right-hand side
    g_ops : ndarray
        sequence of operations for Gaussian elimination
    g_ops : ndarray
        sequence of operations for forward substitution

    Returns
    -------
    x : ndarray
        solution
    """
    gauss_elim(Q, b, g_ops)
    x = b
    forward_sub(Q, x, f_ops)
    return x


# def solve_grad_l5(Q, B, g_ops, f_ops):
#     """Solve linear system of matrix equations QX=B with Gaussian elimination
#     (using v[0] as soma requires clearing upper triangle first). Note: modifies
#     B in place to produce solution X.

#     Parameters
#     ----------
#     Q : ndarray
#         coefficient matrix
#     B : ndarray
#         right-hand side
#     g_ops : ndarray
#         sequence of operations for Gaussian elimination
#     f_ops : ndarray
#         sequence of operations for forward substitution

#     """
#     gauss_elim_mat(Q, B, g_ops)
#     X = B
#     forward_sub_mat(Q, X, f_ops)
#     return X

def solve_grad(Q, B, g_ops, f_ops):
    """Solve linear system of matrix equations QX=B with Gaussian elimination
    (using v[0] as soma requires clearing upper triangle first). Note: modifies
    B in place to produce solution X.

    Parameters
    ----------
    Q : ndarray
        coefficient matrix
    B : ndarray
        right-hand side
    g_ops : ndarray
        sequence of operations for Gaussian elimination
    f_ops : ndarray
        sequence of operations for forward substitution

    """
    gauss_elim_mat(Q, B, g_ops)
    X = B
    forward_sub_mat(Q, X, f_ops)


@nb.jit(nopython=True, cache=True)
def gauss_elim(Q, b, g_ops):
    """Gaussian elimination (upper triangle cleared)

    Parameters
    ----------
    Q : ndarray
        coefficient matrix
    b : ndarray
        right-hand side
    g_ops : ndarray
        sequence of operations for Gaussian elimination
    """
    for k in range(g_ops.shape[0]):
        p = g_ops[k][0]
        t = g_ops[k][1]
        if t != p-1:
            b[t] -= Q[t, p]/Q[p, p]*b[p]
            Q[t, :] = Q[t, :] - Q[t, p]/Q[p, p]*Q[p, :]
        else:
            b[t] -= Q[t, p]/Q[p, p]*b[p]
            Q[t, t] = Q[t, t] - Q[t, p]/Q[p, p]*Q[p, t]
            Q[t, p] = 0


@nb.jit(nopython=True, cache=True)
def gauss_elim_mat(Q, B, g_ops):
    """Gaussian elimination (upper triangle cleared) for matrix system

    Parameters
    ----------
    Q : ndarray
        coefficient matrix
    B : ndarray
        right-hand side
    g_ops : ndarray
        sequence of operations for Gaussian elimination
    """
    for k in range(g_ops.shape[0]):
        p = g_ops[k][0]
        t = g_ops[k][1]
        if t != p-1:
            B[t, :] = B[t, :] - Q[t, p]/Q[p, p]*B[p, :]
            Q[t, :] = Q[t, :] - Q[t, p]/Q[p, p]*Q[p, :]
        else:
            B[t, :] = B[t, :] - Q[t, p]/Q[p, p]*B[p, :]
            Q[t, t] = Q[t, t] - Q[t, p]/Q[p, p]*Q[p, t]
            Q[t, p] = 0


@nb.jit(nopython=True, cache=True)
def row_reduce(Q, g_ops):
    """Row reduction of Q to precompute forward subtitution operations.

    Parameters
    ----------
    Q : ndarray
        matrix to be reduced
    g_ops : ndarray
        sequence of operations for Gaussian elimination of Q
    """
    for k in range(g_ops.shape[0]):
        p = g_ops[k][0]
        t = g_ops[k][1]
        Q[t, :] = Q[t, :] - Q[t, p]/Q[p, p]*Q[p, :]


@nb.jit(nopython=True, cache=True)
def forward_sub(Q, x, f_ops):
    """Forward substitution after gauss_elim.

    Parameters
    ----------
    Q : ndarray
        row-reduced matrix
    x : ndarray
        view to rhs b.
    f_ops : ndarray
        sequence of operations for forward substitution
    """
    x /= np.diag(Q)
    for k in range(f_ops.shape[0]):
        r = f_ops[k][0]
        c = f_ops[k][1]
        x[r] -= Q[r, c]/Q[r, r]*x[c]


@nb.jit(nopython=True, cache=True)
def forward_sub_mat(Q, X, f_ops):
    """Forward substitution for matrix system after gauss_elim_mat

    Parameters
    ----------
    Q : ndarray
        row-reduced matrix
    X : ndarray
        view to rhs B.
    f_ops : ndarray
        sequence of operations for forward substitution
    """
    q = np.expand_dims(np.diag(Q), 1)
    X /= q
    for k in range(f_ops.shape[0]):
        r = f_ops[k][0]
        c = f_ops[k][1]
        X[r, :] -= Q[r, c]/Q[r, r]*X[c, :]
