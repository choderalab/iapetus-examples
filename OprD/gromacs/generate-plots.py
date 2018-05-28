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

#from yank.multistate import MultiStateReporter

def generate_plots(filename, refpdb_filename):
    print(filename)
    if not os.path.exists(filename):
        return
    with Dataset(filename, 'r') as ncfile:
        logZ = np.array(ncfile.groups['online_analysis'].variables['logZ_history'])
        logZ = logZ[:-1,:]
        n_iterations, n_states = logZ.shape
        print(n_iterations, n_states)

        states = ncfile.variables['states']

        gamma = ncfile.groups['online_analysis'].variables['gamma_history']

        pdf_filename = 'output.pdf'
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

            fig = plt.figure(figsize=(10,5));
            plt.plot(gamma,'.');
            plt.xlabel('iteration');
            plt.ylabel('gamma');
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

        # Extract trajectory
        print('Extracting positions...')
        mdtraj_refpdb = md.load(refpdb_filename)
        solute_indices = mdtraj_refpdb.topology.select('not water')
        solute_topology = mdtraj_refpdb.topology.subset(solute_indices)
        mdtraj_solute_trajectory = md.Trajectory(ncfile.variables['positions'][:-2,0,solute_indices,:], solute_topology)
        mdtraj_solute_trajectory[0].save('output.pdb')
        mdtraj_solute_trajectory.save('output.xtc')

#generate_plots('output.nc', 'comp7_nowat/3sy7_lig_nowat_GMX.pdb')
generate_plots('comp7.nc', 'comp7/3sy7_lig_solv_GMX.pdb')
