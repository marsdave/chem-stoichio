# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 13:19:00 2021

@author: vogt_dv
"""

from chempy import mass_fractions
from chempy import Substance
import numpy as np

from functions import comp_from_mix
from functions import mixture_stoichiometry
from functions import mixture_weights

# EXAMPLE 1: Calculate stoichiometry of a mixture of compounds given in mass fractions
print("EXAMPLE 1")
# 1) Set up a Python dictionairy where the entries have the names of chemical compounds, and the values are their weights.
mixture_dict = {'MgO': 0.3, 'SiO2': 0.3, 'H2O': 0.4}
# 2) The output of the function mixture_stoichiometry is a dictionairy containing the atomic fractions of the elements in the mixture.
stoichio_dict = mixture_stoichiometry(mixture_dict)
print(stoichio_dict)
# 3) The entries can be accessed by giving the name of the element as the key
print(stoichio_dict['H'])
# 4) The opposite direction is available with the mass_fractions function by ChemPy
stoichio_water = {'H': 2, 'O': 1}
mass_frac = mass_fractions(stoichio_water)
print(mass_frac)


# EXAMPLE 2: Calculate the mixture in mass fractions for a given set of atomic concentrations
print("EXAMPLE 2")
# 1) Define the desired atomic fractions of certain elements as a dictionairy
atomic_fractions_dict = {'Ca': 0.015, 'F': 0.03, 'Cl': 0.015}
# 2) Define the substances in the mixture
substances_list = ['CaF2', 'K2SO4', 'KCl']
# 3) Minimize a fit function with mixture_weights to find the best mixture of the substances
mix_wt = mixture_weights(atomic_fractions_dict, substances_list)
print(mix_wt)
# 4) Check to see if the desired atomic fractions have been found (small deviations are possible)
mix_at = mixture_stoichiometry(mix_wt)
print({key: mix_at[key] for key in atomic_fractions_dict.keys()})


# EXAMPLE 3: Calculate stoichiometry of a mixture of Mars regolith simulant JEZ-1 and olivine
print("EXAMPLE 3")
# 1) Define the major oxides as ChemPy substances
major_oxides_list = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'SO3']
major_oxides_substance_dict = {s: Substance.from_formula(s) for s in major_oxides_list}

# 2) Define Mars Regolith Simulant JEZ-1 (Jezero crater soil simulant) by the weight percentages of its major oxides
mix_jezero = {
    'SiO2': 44.2,
    'TiO2': 0.2,
    'Al2O3': 11.3,
    'Cr2O3': 0.3,
    'FeO': 9.5,
    'MnO': 0.1,
    'MgO': 25.9,
    'CaO': 3.5,
    'Na2O': 1.9,
    'K2O': 0.3,
    'P2O5': 0.6,
    'SO3': 2.1,
    }

# 3) Define the mineral olivine, (Mg,Fe)2SiO4, as a mixture containing different pseudo chemical formulas
mix_olivine = {
    'Mg2SiO4': 0.5,
    'Fe2SiO4': 0.5
    }

# 4) Create ChemPy substances from the mixtures and create a dictionairy of all substances within those two substances
subst_jezero = Substance(name='Exolith JEZ-1', composition=comp_from_mix(mix_jezero))
subst_olivine = Substance(name='Olivine', composition=comp_from_mix(mix_olivine))
substances = {s.name: s for s in [subst_jezero, subst_olivine]}

# 5) Calculate the stoichiometry of the mixture of the two substances for different weight ratios
wt_olivine = np.linspace(0, 1, num=21)  # weight percentage of olivine in the mixture
for i in range(len(wt_olivine)):
    mix_wt = {'Olivine': wt_olivine[i], 'Exolith JEZ-1': 1-wt_olivine[i]}
    mix_stoi = mixture_stoichiometry(mix_wt, substances=substances)
    wt_str = str(int(wt_olivine[i] * 100)) + ' wt%'
    print(wt_str, np.round((mix_stoi['Fe'] + mix_stoi['Mg'])/mix_stoi['Si'], decimals=3), sep='\t')

# substances_list = ['NaCl', 'CaSO4(H2O)2']
# atperc_Cl_arr = np.array([0, 0.25, 0.75, 1.26, 2.59, 5.48])
# for atperc in atperc_Cl_arr:
#     atomic_fractions_dict = {'Cl': atperc/100}
#     mix_wt = mixture_weights(atomic_fractions_dict, substances_list)
#     mix_at = mixture_stoichiometry(mix_wt)
#     print('NaCl: ' + str(np.round(mix_wt['NaCl']*100, decimals=2)) + ' wt%, Cl: ' + str(np.round(mix_at['Cl']*100, decimals=2)) + ' at%')




# element_dict = atomicnumber_to_element_dict()
# for key in element_dict:
#     subst = Substance.from_formula(element_dict[key])
#     print(key, subst.name, subst.mass)

