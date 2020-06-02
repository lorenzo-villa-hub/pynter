
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 15})

class JobsAnalysis:
    
    def __init__(self,dataset):
        """
        Class to analyse multiple Jobs groupped in a Dataset

        Parameters
        ----------
        dataset : (Datset object)
        """
        self.ds = dataset

        
    def plot_fractional_charge(self):
        """
        Plot fractional charge study from data in Dataset

        Returns
        -------
        plt : 
            Matplotlib object
        """
        plt.figure(figsize=(8,6))     
        charges = []
        energies = []        
        jobs = self.ds.jobs
        energy_dict = {}
        
        for j in jobs:
            n = -1*np.around(j.charge(),decimals=1) # express in terms of occupation
            energy_dict[n] = j.final_energy()            
        energy_dict = {k: v for k, v in sorted(energy_dict.items(), key=lambda item: item[0])} #order by charge

        def linear(x,E0,E1):
            return (E1-E0)*x + E0
        
        for n in energy_dict:
            e = float(energy_dict[n])
            e0 = float(energy_dict[0])
            e1 = float(energy_dict[1])
            charges.append(n)
            e_rescaled = e - linear(n, e0, e1) 
            energies.append(e_rescaled)

        plt.plot(charges,energies,'o--',label=self.ds.name)
        plt.hlines(0,xmin=0,xmax=1,linestyles='dashed',label='Exact',color='k')
        width = max([abs(e) for e in energies])
        plt.ylim(-3*width,+3*width)
        plt.xlabel('Fractional charge')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.grid()
            
        return plt
            