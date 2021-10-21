# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 18:31:24 2021

@author: vogt_dv
"""

from chempy import mass_fractions
from chempy import Substance
from chempy import balance_stoichiometry
import numpy as np
from scipy.optimize import least_squares

def atomicnumber_to_element_dict(reverse=False):
    element_list = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 
                    'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 
                    'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 
                    'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 
                    'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 
                    'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 
                    'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
                    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 
                    'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 
                    'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 
                    'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 
                    'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
    
    if reverse == False:
        output_dict = {}
        for i in range(len(element_list)):
            output_dict[i + 1] = element_list[i]
    else:
        output_dict = {}
        for i in range(len(element_list)):
            output_dict[element_list[i]] = i + 1
    
    return output_dict

def elements_to_numbers(stoichio_dict):
    element_dict_reverse = atomicnumber_to_element_dict(reverse=True)
    output_dict = {element_dict_reverse[key]: stoichio_dict[key] for key in stoichio_dict.keys()}
    return output_dict

def comp_from_mix(mixture_dict):
    stoichio_mix = mixture_stoichiometry(mixture_dict)
    comp = elements_to_numbers(stoichio_mix)
    return comp

def mixture_stoichiometry(mixture_dict, normalize=True, substances=None):
    element_dict = atomicnumber_to_element_dict()
    
    if substances == None:
        substances = {key: Substance.from_formula(key) for key in mixture_dict.keys()}
    
    # extract all keys for elements present in the mixture, create dictionary for them
    substances_list = [substances[key] for key in substances.keys()]
    element_keys = Substance.composition_keys(substances_list)
    stoichio_dict = {}
    for key in element_keys:
        stoichio_dict[element_dict[key]] = 0
    
    for key in mixture_dict.keys():
        molar_mass = substances[key].mass   # molar mass in g/mol
        weight = mixture_dict[key]          # weight of substance in g
        mol = weight/molar_mass             # mol of substance
        # for each element, add the number of atoms in the formula times the mol number of the substance to the stoichiometry
        for composition_key in substances[key].composition.keys():
            stoichio_dict[element_dict[composition_key]] += substances[key].composition[composition_key] * mol
    
    if normalize:
        sum_stoichio = sum(stoichio_dict.values())
        for key in stoichio_dict.keys():
            stoichio_dict[key] /= sum_stoichio
    
    return stoichio_dict

def fitfunc_mixture_weights(wt_frac, substances, atomic_fractions):
    # wt_frac = weights/np.sum(weights)
    
    mixture_dict = {}
    for i in range(len(substances)):
        mixture_dict[substances[i]] = wt_frac[i]
    
    stoichio = mixture_stoichiometry(mixture_dict, normalize=True)
    
    goal_keys = list(atomic_fractions.keys())  # list of elements of interest (e.g. ['H', 'O'])
    
    residuals = np.zeros(len(goal_keys))
    for i in range(len(goal_keys)):
        key = goal_keys[i]
        residuals[i] = stoichio[key] - atomic_fractions[key]
    
    return residuals
    

def mixture_weights(atomic_fractions_dict, substances_list):
    x0 = np.ones(len(substances_list))/len(substances_list)
    # x0 = np.zeros(len(substances_list))
    fit = least_squares(fitfunc_mixture_weights, x0, bounds=(0,1), 
                        args=(substances_list, atomic_fractions_dict), 
                        gtol=1e-12, diff_step=1e-16)
    
    mixture_dict = {}
    for i in range(len(substances_list)):
        key = substances_list[i]
        mixture_dict[key] = fit.x[i]/np.sum(fit.x)
    
    return mixture_dict

# # Set up a Python dictionairy where the entries have the names of chemical compounds, and the values are their weights.
# mixture_dict = {'MgO': 0.3, 'SiO2': 0.3, 'H2O': 0.4}

# # The output of the function is a dictionairy containing the atomic fractions of the elements in the mixture.
# stoichio_dict = mixture_stoichiometry(mixture_dict)
# print(stoichio_dict)

# # The entries can be accessed by giving the name of the element as the key
# print(stoichio_dict['H'])

# # The opposite direction is available with the mass_fractions function by ChemPy
# stoichio_water = {'H': 2, 'O': 1}
# mass_frac = mass_fractions(stoichio_water)
# print(mass_frac)



# atomic_fractions_dict = {'Ca': 1.5, 'F': 0.03, 'Cl': 0.015}
# substances_list = ['CaF2', 'K2SO4', 'KCl']

# mix_wt = mixture_weights(atomic_fractions_dict, substances_list)
# print(mix_wt)
# mix_at = mixture_stoichiometry(mix_wt)
# print({key: mix_at[key] for key in atomic_fractions_dict.keys()})

major_oxides_list = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'Fe2O3', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'SO3']
major_oxides_substance_dict = {s: Substance.from_formula(s) for s in major_oxides_list}

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

mix_bp1 = {
    'SiO2': 47.2,
    'TiO2': 2.3,
    'Al2O3': 16.7,
    'FeO': 6.2,
    'Fe2O3': 5.9,
    'MnO': 0.21,
    'MgO': 6.5,
    'CaO': 9.2,
    'Na2O': 3.5,
    'K2O': 1.1,
    'P2O5': 0.52,
    'SO3': 2.1,
    'H2O': 0.46,
    'CO2': 0.05
    }

mix_olivine = {
    'Mg2SiO4': 0,
    'MgFeSiO4': 1,
    'Fe2SiO4': 1
    }

subst_jezero = Substance(name='Exolith JEZ-1', composition=comp_from_mix(mix_jezero))
subst_olivine = Substance(name='Olivine', composition=comp_from_mix(mix_olivine))
substances = {s.name: s for s in [subst_jezero, subst_olivine]}

wt_olivine = np.linspace(0, 1, num=21)
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

