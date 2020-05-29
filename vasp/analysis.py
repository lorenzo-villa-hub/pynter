
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
        
        for n in energy_dict:
            E = energy_dict[n]
            if n == 0:
                E0 = E
            if n == 1:
                E1 = E
            charges.append(n)
            energies.append(E0-E)
        
        plt.plot(charges,energies,label=self.ds.name)
        plt.plot([0,1],[0,E0-E1],label='Ideal')
        plt.xlabel('Fractional charge')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.grid()
            
        return plt
            