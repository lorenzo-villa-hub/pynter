
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


class Vaspout: # in progress, will probably be eliminated because inefficient
    
    def __init__(self,**kwargs):
        
        for key,value in kwargs.items():
            setattr(self, key, value)
 
           
    def __str__(self):
        return self.__dict__.__str__()
    
    def __repr__(self):
        return self.__str__()
    
    
    def as_dict(self):
        d = {}
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
    #    d['structures'] = [s.as_dict() for s in self.structures]
        d['epsilon_static'] = self.epsilon_static
        d['epsilon_static_wolfe'] = self.epsilon_static_wolfe
        d['epsilon_ionic'] = self.epsilon_ionic
 #       d['dielectric'] = self.dielectric if hasattr(self,'dielectric') else None
 #       d['optical_absorption_coeff'] = self.optical_absorption_coeff if hasattr(self,'optical_absorption_coeff') else None
        d['converged_electronic'] = self.converged_electronic
        d['converged_ionic'] = self.converged_ionic
        d['converged'] = self.converged
        d['final_energy'] = self.final_energy
        d['complete_dos'] = self.complete_dos.as_dict()
        d['hubbards'] = self.hubbards
        d['is_hubbard'] = self.is_hubbard
        d['is_spin'] = self.is_spin
        d['bandstructure'] = self.bandstructure.as_dict()
        d['eigenvalue_band_properties'] = self.eigenvalue_band_properties
        
        return d
        
        
    @staticmethod
    def from_dict(d):
        k = {}
        k['structures'] = [Structure.from_dict(s) for s in d['structures']]
        pass
        

    
    @staticmethod
    def from_Vasprun(vasprun):
        v = vasprun
        k = {}
        k['structures'] = v.structures
        k['epsilon_static'] = v.epsilon_static
        k['epsilon_static_wolfe'] = v.epsilon_static_wolfe
        k['epsilon_ionic'] = v.epsilon_ionic
        if v.dielectric_data:
            k['dielectric'] = v.dielectric
            k['optical_absorption_coeff'] = v.optical_absorption_coeff
        k['converged_electronic'] = v.converged_electronic
        k['converged_ionic'] = v.converged_ionic
        k['converged'] = v.converged
        k['final_energy'] = v.final_energy
        k['complete_dos'] = v.complete_dos
        k['hubbards'] = v.hubbards
        k['is_hubbard'] = v.is_hubbard
        k['is_spin'] = v.is_spin
        k['bandstructure'] = v.get_band_structure(force_hybrid_mode=True)
        k['eigenvalue_band_properties'] = v.eigenvalue_band_properties
                
        kwargs = k
            
        return Vaspout(**kwargs)
        
        

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
            vasprun = job.outputs['Vasprun']       
            complete_dos = vasprun.complete_dos
            partial_dos = complete_dos.get_spd_dos()        
            dos_plotter = DosPlotter(stack=stack)
            dos_plotter.add_dos('Total',complete_dos)
            for orbital in partial_dos:
                dos_plotter.add_dos(orbital,partial_dos[orbital])
            eg = job.energy_gap()
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
            vasprun = job.outputs['Vasprun']
            bs = vasprun.get_band_structure(line_mode=True)
            dos = vasprun.complete_dos
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
        V,E = [],[]
        for j in self.jobs:
            V.append(j.initial_structure.lattice.volume)
            E.append(j.final_energy())
        eos = EOS(eos_name='birch_murnaghan')
        eos_fit = eos.fit(V, E)
        eos_fit.plot(width = 10, height = 10 , text = '', markersize = 15,  label= 'Birch-Murnaghan fit')  
        plt.legend(loc=2, prop={'size': 20})
        if title:
            plt.title(title,size=25)
        plt.tight_layout()
        
        return plt, eos_fit

        
    def plot_fractional_charge(self,name='',new_figure=True,legend_out=False):
        """
        Plot fractional charge study
        
        Parameters
        ----------
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
            n = -1*np.around(j.charge(),decimals=1) # express in terms of occupation
            energy_dict[n] = j.final_energy()            
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
            U_list.append(j.hubbard()[el])
            lattice = j.final_structure().lattice
            lattice_params.append(getattr(lattice,lattice_parameter))
            band_gaps.append(j.energy_gap())
        
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
            
                
        
        
        
        
        
        