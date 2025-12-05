# Load a csv file and fit the model to the data
import numpy as np
import dill
from lmfit import Model, Parameters
from scipy import constants as spc
import matplotlib.pyplot as plt

from model import *

k_B = spc.Boltzmann #/ spc.e  # Boltzmann constant in eV/K
e = spc.e  # Elementary charge in Coulombs
h = spc.h  # Planck constant in J.s
m_0 = spc.m_e  # Electron rest mass in kg
m_star_inas = 0.023 * m_0  # Effective mass for InAs in kg

A_star_inas = 4 * np.pi * e * m_star_inas * k_B**2 / h**3  # Richardson constant for InAs in A/m^2/K^2

A_contact = 200e-9 * 300e-9 * 3/2  # Effective contact area in m^2

A_star = A_star_inas * A_contact  # Effective Richardson constant in A/K^2


data = np.genfromtxt('example_data.csv',comments='#', delimiter=',', skip_header=1)
data = data[1:, :]
I = data[:, 0]  # Current in Amperes
V = data[:, 1]  # Voltage in Volts


model = Model(Back_to_back_schottky_diode_current_eff_barrier_for_fit)

model.set_param_hint('B_1', value=0.4, vary=True, min=0.000001, max=2)
model.set_param_hint('B_2', value=0.4, vary=True, min=0.000001, max=2)
model.set_param_hint('n_1', value=1.2, vary=True, min=1, max=2)
model.set_param_hint('n_2', value=1.2, vary=True, min=1, max=2)
model.set_param_hint('R', value=4000000, vary=True, min=1000, max=1e15)
model.set_param_hint('T', value=300, vary=False, min=0, max=1000)
model.set_param_hint('A_star', value=A_contact*A_star_inas, vary=False)

model.print_param_hints()
result = model.fit(I, V=V)
result.params.pretty_print()
print(result.fit_report())

with open('example.dill', 'wb') as f:
    dill.dump(result, f)

