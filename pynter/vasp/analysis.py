
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.plotter import BSPlotter,BSDOSPlotter,DosPlotter
from pymatgen.electronic_structure.bandstructure import BandStructure
from pymatgen.analysis.eos import EOS

matplotlib.rcParams.update({'font.size': 15})


class JobAnalysis:
    
    def __init__(self,job):
        """
        Class to analyse a VaspJob object

        Parameters
        ----------
        job : 
            Job object
        """
        
        self.job = job
        
    def plot_dos(self,xlim=(-3,3),stack=False):
        """
        Plot DOS from data in vasprun.xml with Pymatgen
        """
        job = self.job
        wdir = os.getcwd()
        os.chdir(job.path)
        if job.is_converged:   
            complete_dos = job.complete_dos
            partial_dos = complete_dos.get_spd_dos()        
            dos_plotter = DosPlotter(stack=stack)
            dos_plotter.add_dos('Total',complete_dos)
            for orbital in partial_dos:
                dos_plotter.add_dos(orbital,partial_dos[orbital])
            eg = job.energy_gap
            plt = dos_plotter.get_plot(xlim=(xlim[0],eg+xlim[1]))
        else:
            raise ValueError(f'Job %s is not converged' %self.job.name)
        os.chdir(wdir)
        return plt


    def plot_dos_bs(self):
        """
        Plot DOS and BS from data in vasprun.xml with Pymatgen
        """
        job = self.job
        wdir = os.getcwd()
        os.chdir(job.path)
        if job.is_converged:
            bs = job.get_band_structure(line_mode=True)
            dos = job.complete_dos
            plt = BSDOSPlotter(bs_projection=None,dos_projection=None).get_plot(bs,dos)           
        else:
            raise ValueError(f'Job %s is not converged' %self.name)
        os.chdir(wdir)
        return plt


class DatasetAnalysis:
    
    def __init__(self,jobs):
        """
        Class to analyse multiple Jobs

        Parameters
        ----------
        jobs : 
            List of VaspJob objects
        """
        self.jobs = jobs


    def birch_murnaghan(self,title=None):
        """
        Fit E,V data with Birch-Murnaghan EOS

        Parameters
        ----------
        title : (str)
            Plot title

        Returns
        -------
        plt :
            Matplotlib object.
        eos_fit : 
            Pymatgen EosBase object.
        """
        V,E = [],[]
        for j in self.jobs:
            V.append(j.initial_structure.lattice.volume)
            E.append(j.final_energy)
        eos = EOS(eos_name='birch_murnaghan')
        eos_fit = eos.fit(V, E)
        eos_fit.plot(width = 10, height = 10 , text = '', markersize = 15,  label= 'Birch-Murnaghan fit')  
        plt.legend(loc=2, prop={'size': 20})
        if title:
            plt.title(title,size=25)
        plt.tight_layout()
        
        return plt, eos_fit


    def plot_convergence_encut(self,title=None,delta=1e-03):
        """
        Plot total energy convergence with respect to ENCUT

        Parameters
        ----------
        title : (str), optional
            Plot title. The default is None.
        delta : (float), optional
            Value of delta E for energy convergence. The default is 1e-03.

        Returns
        -------
        plt : 
            Matplotlib object.
        """
        plt.figure(figsize=(10,8))
        cutoff_max = 0
        energies = {}
        for j in self.jobs:
            cutoff = j.incar['ENCUT']
            if cutoff > cutoff_max:
                cutoff_max = cutoff
            energies[cutoff] = j.final_energy/len(j.initial_structure)
        
        E0 = energies[cutoff_max]
        cutoffs = [c for c in energies.keys()]
        deltas = [E0 - e for e in energies.values()]

        plt.plot(cutoffs,deltas,'o--',markersize=10)
        plt.axhline(y=delta,ls='-',color='k')
        plt.axhline(y=-1*delta,ls='-',color='k')
        plt.xlabel('PW cutoff (eV)')
        plt.ylabel('Energy difference (eV/atom)')
        plt.grid()

        return plt         

        
    def plot_fractional_charge(self,reference='electrons',name='',new_figure=True,legend_out=False):
        """
        Plot fractional charge study
        
        Parameters
        ----------
        reference: (str)
            Reference for fractional charge sign. 'electrons' or 'holes'. Default is 'electrons'.
        name: (str)
            Name to put in legend
        new_figure: (bool)
            If true initialize new matplotlib figure
        legend_out: (bool)
            If True place legend outside of the plot

        Returns
        -------
        plt : 
            Matplotlib object
        """
        if new_figure:
            plt.figure(figsize=(8,6))            
        jobs = self.jobs
        energy_dict = {}
        
        for j in jobs:
            if reference=='electrons':
                n = -1*np.around(j.charge,decimals=1) # express in terms of occupation
            if reference=='holes':
                n = 1 - np.around(j.charge,decimals=1) #shift to 0 charge of +1
            energy_dict[n] = j.final_energy            
        energy_dict = {k: v for k, v in sorted(energy_dict.items(), key=lambda item: item[0])} #order by charge

        def linear(x,E0,E1):
            return (E1-E0)*x + E0
        
        e_rescaled = {}
        e0 = float(energy_dict[0])
        e1 = float(energy_dict[1])
        for n in energy_dict:
            e = float(energy_dict[n])
            e_resc = e - linear(n, e0, e1)
            e_rescaled[n] = e_resc

        charges = list(e_rescaled.keys())
        energies = list(e_rescaled.values())

        ax = plt.gca()
        if not ax.lines:
            plt.hlines(0,xmin=0,xmax=1,linestyles='dashed',label='Exact',color='k')
        plt.plot(charges,energies,'o--',label=name)
        if new_figure:
            width = max([abs(e) for e in energies])
            plt.ylim(-3*width,+3*width)
        plt.xlabel('Fractional charge')
        plt.ylabel('Energy (eV)')
        if legend_out:
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            plt.legend()
        plt.grid()
            
        return plt
            
    
    def plot_U_dependency(self,element,lattice_parameter='a'):
        """
        Plot dependency of lattice parameter and energy gap on U parameter assigned
        to a specific element

        Parameters
        ----------
        element : (str)
            String of element to which U parameter is applied.
        lattice_parameter : (str), optional
            Lattice parameter whose dependency is evaluated. The default is 'a'.

        Returns
        -------
        plt :
            Matplotlib object
        """
        el = Element(element)
        plt.figure(figsize=(12,8))
        U_list = []
        lattice_params = []
        band_gaps = []
        for j in self.jobs:
            U_list.append(j.hubbard[el])
            lattice = j.final_structure.lattice
            lattice_params.append(getattr(lattice,lattice_parameter))
            band_gaps.append(j.energy_gap)
        
        plt.subplot(1,2,1)
        plt.xlabel('U on %s'%element)
        plt.ylabel(lattice_parameter+' ($\AA$)')
        plt.plot(U_list,lattice_params,'o--')
        plt.grid()
        
        plt.subplot(1,2,2)
        plt.xlabel('U on %s'%element)
        plt.ylabel('Energy gap (eV)')
        plt.plot(U_list,band_gaps,'o--')
        plt.grid()
        
        return plt
            
                
        
        
        
        
        
        