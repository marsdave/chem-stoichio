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
