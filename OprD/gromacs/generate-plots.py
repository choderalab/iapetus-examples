#!/bin/env python

import matplotlib
matplotlib.use('agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
import os
from netCDF4 import Dataset
import mdtraj as md
from pymbar import timeseries, MBAR
import pymbar
from scipy.misc import logsumexp

#from yank.multistate import MultiStateReporter

def estimate_f_k(ncfile):
    """
    Estimate the free energy of each thermodynamic state
    This is equivalent to estimating the PMF using a Gaussian kernel density estimate

    Parameters
    ----------
    ncfile : netCDF4.Dataset
        The NetCDF file

    Returns
    -------
    f_k : np.array with shape [nstates-1]
        f_k[k] is the relative free energy (in kT) of thermodynamic state k
        The unbiased state is omitted
        The lowest free energy state has free energy 0
    df_k : np.array with shape [nstates-1]
        Uncertainty in f_k relative to lowest free energy state

    """
    niterations, nreplicas, nstates = ncfile.variables['energies'].shape
    niterations -= 1
    print('niterations = {}'.format(niterations))
    print('nreplicas = {}'.format(nreplicas))
    print('nstates = {}'.format(nstates))

    # Use the state correlation time for the statistical statisticalInefficiency
    # TODO: Use u_n = u(x_n; s_n) + g_{s_n}  timeseries instead
    # TODO: Generalize to support multiple replicas?
    replica_index = 0
    states_t = ncfile.variables['states'][:,replica_index]
    [t0, g, Neff_max] = timeseries.detectEquilibration(states_t, fast=True, nskip=10)
    print(t0, g, Neff_max)

    # Use log deviance instead
    u_n = np.zeros([niterations], np.float32)
    g_k_final = - np.array(ncfile.groups['online_analysis'].variables['logZ'][:]) # last (best) g_k estimate
    states = np.array(ncfile.variables['states'][:,:], np.int32)
    energies = np.array(ncfile.variables['energies'][:,:])
    logZ = np.array(ncfile.groups['online_analysis'].variables['logZ_history'][:,:])
    for iteration in range(niterations):
        g_k = - logZ[iteration,:]
        state_index = states[iteration,replica_index]
        u_n[iteration] = energies[iteration,replica_index,state_index] - g_k[state_index] + logsumexp(-g_k_final + g_k)
    [t0, g, Neff_max] = timeseries.detectEquilibration(u_n, fast=True, nskip=10)
    print(t0, g, Neff_max)

    # DEBUG
    #print('Resetting t0 and g')
    #t0 = 0
    #g = 5

    # Subsample energies
    indices = np.array(timeseries.subsampleCorrelatedData(states_t[t0:niterations], g=g)) + t0
    u_kn = np.array(ncfile.variables['energies'][indices,0,:].transpose(), np.float64)
    states_n = ncfile.variables['states'][indices,0]
    N_k = np.zeros([nstates], np.int32)
    for state in states_n:
        N_k[state] += 1
    print(N_k)
    u_kn -= u_kn.min()
    print(u_kn)
    print(np.isnan(u_kn).sum())

    # Estimate f_k
    print('Estimating MBAR...')
    mbar = MBAR(u_kn, N_k, relative_tolerance=1.0e-10)
    delta_f_ij, ddelta_f_ij, theta = mbar.getFreeEnergyDifferences(uncertainty_method='svd-ew')

    # Return free energy with respect to minimum free energy state
    reference_state_index = np.argmin(delta_f_ij[0,:])
    f_k = delta_f_ij[reference_state_index,:]
    df_k = ddelta_f_ij[reference_state_index,:]

    return f_k, df_k, u_n

def generate_plots(filename, refpdb_filename, prefix):
    print(filename)
    if not os.path.exists(filename):
        return
    with Dataset(filename, 'r') as ncfile:
        logZ = np.array(ncfile.groups['online_analysis'].variables['logZ_history'])
        logZ = logZ[:-1,:]
        n_iterations, n_replicas, n_states = ncfile.variables['energies'].shape
        #n_iterations, n_states = logZ.shape
        print(n_iterations, n_states)

        states = ncfile.variables['states']

        gamma = ncfile.groups['online_analysis'].variables['gamma_history']

        pdf_filename = prefix + '.pdf'
        with PdfPages(pdf_filename) as pdf:
            fig = plt.figure(figsize=(10,5));
            plt.plot(logZ[:,:],'-')
            plt.xlabel('iteration');
            plt.ylabel('logZ / kT');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            fig = plt.figure(figsize=(10,5));
            plt.plot(logZ[:,0] - logZ[:,-1],'.');
            plt.xlabel('iteration');
            plt.ylabel('$\Delta G_\mathrm{complex}$ / kT');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            fig = plt.figure(figsize=(10,5));
            plt.plot(states,'.');
            plt.xlabel('iteration');
            plt.ylabel('thermodynamic state');
            plt.axis([0, n_iterations, 0, n_states]);
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            # invert y-axis
            fig = plt.figure(figsize=(10,5));
            plt.plot(states,'.');
            plt.xlabel('iteration');
            plt.ylabel('thermodynamic state');
            plt.axis([0, n_iterations, 0, n_states]);
            plt.gca().invert_yaxis()
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            fig = plt.figure(figsize=(10,5));
            plt.plot(gamma,'.');
            plt.xlabel('iteration');
            plt.ylabel('gamma');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            fig = plt.figure(figsize=(10,5));
            plt.plot(-logZ[-1,:],'.');
            plt.xlabel('thermodynamic state');
            plt.ylabel('-logZ');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            # TODO: Estimate PMF with MBAR
            f_k, df_k, u_n = estimate_f_k(ncfile)

            fig = plt.figure(figsize=(10,5));
            plt.plot(u_n,'.');
            plt.xlabel('iteration');
            plt.ylabel('u_n');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

            x = np.arange(n_states-1)
            fig = plt.figure(figsize=(10,5));
            plt.fill_between(x, f_k - 2*df_k, f_k + 2*df_k, color=[0.5, 0.5, 0.5])
            plt.plot(f_k, '.');
            plt.xlabel('pore distance');
            plt.ylabel('free energy / kT');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        # Extract trajectory
        print('Extracting positions...')
        mdtraj_refpdb = md.load(refpdb_filename)
        solute_indices = mdtraj_refpdb.topology.select('not water')
        solute_topology = mdtraj_refpdb.topology.subset(solute_indices)
        mdtraj_solute_trajectory = md.Trajectory(ncfile.variables['positions'][:-2,0,solute_indices,:], solute_topology)
        mdtraj_solute_trajectory[0].save(prefix + '.pdb')
        mdtraj_solute_trajectory.save(prefix + '.xtc')

#generate_plots('output.nc', 'comp7_nowat/3sy7_lig_nowat_GMX.pdb')
jobnames = ['arg', 'comp7', 'comp7_nowat', 'comp8', 'comp8_nowat', 'glu', 'his_neut', 'his_posit', 'imi', 'mero']
#jobnames = ['imi']
for prefix in jobnames:
    print(prefix)
    try:
        generate_plots(prefix + '.nc', prefix + '/3sy7_lig_solv_GMX.pdb', prefix)
    except pymbar.utils.ParameterError as e:
        print(e)
    except OSError as e:
        print(e)
