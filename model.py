import numpy as np
from scipy import constants as spc
import scipy.optimize as opt

k_B = spc.Boltzmann #/ spc.e  # Boltzmann constant in eV/K
e = spc.e  # Elementary charge in Coulombs

def schottky_diode_current_ideal(V, B, T, n, A_star):
    """Calculate the current through a Schottky diode using the thermionic emission theory.

    Parameters:
    V (float): Applied voltage in Volts.
    B (float): Schottky barrier height in eV.
    T (float): Temperature in Kelvin.
    n (float): Ideality factor.
    A_star (float): Richardson constant in A/m^2/K^2.

    Returns:
    float: Current in A.
    """
    J0 = A_star * T**2 * np.exp(-e * B / (k_B * T))
    J = J0 * (np.exp(e * V / (n * k_B * T)) - 1)
    return J

def schottky_diode_current_eff_B(V, B, T, A_star, n):
    """Calculate the current through a Schottky diode using the thermionic emission theory with effective barrier lowering.

    Parameters:
    V (float): Applied voltage in Volts.
    B (float): Schottky barrier height in eV.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.
    n (float): Ideality factor.

    Returns:
    float: Current in A.
    """


    B = B - V * (1/n - 1)
    # if B < 0:
    #     B = 1e-50  # prevent negative barrier

    J0 = A_star * T**2 * np.exp(-e * B / (k_B * T))
    if J0 == 0.0:
        J0 = 1e-50  # prevent zero current
    J = J0 * (np.exp(e * V / (k_B * T)) - 1)
    return J


def schottky_diode_voltage_eff_B(I, B, T, A_star, n, min_V=-1e0, max_V=1e0, max_iter=10000, tol=1e-6):
    """Calculate the voltage across a Schottky diode for a given current using the thermionic emission theory with effective barrier lowering.

    Parameters:
    I (float): Current in A.
    B (float): Schottky barrier height in eV.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.
    n (float): Ideality factor.
    I_min, I_max (floats): Lower and upper bounds for current (must bracket the root)
    tol (float) Convergence tolerance
    max_iter (int): Maximum number of iterations


    Returns:
    float: Voltage in Volts.
    """

    def g(V):
        return schottky_diode_current_eff_B(V, B, T, A_star, n) - I
    
    try:
        return opt.bisect(g, min_V, max_V, xtol=tol, maxiter=max_iter)
    except ValueError:
        if g(min_V) > 0 and g(max_V) > 0:
            return min_V
        elif g(min_V) < 0 and g(max_V) < 0:
            return max_V
        
        # raise ValueError("Function does not change sign between V_min and V_max — adjust bounds! g(V_min)={}, g(V_max)={}".format(g(min_V), g(max_V)))
    

def schottky_diode_voltage_ideal(I, B, T, A_star):
    """Calculate the voltage across a Schottky diode for a given current using the thermionic emission theory.

    Parameters:
    I (float): Current in A.
    B (float): Schottky barrier height in eV.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.

    Returns:
    float: Voltage in Volts.
    """

    J0 = A_star * T**2 * np.exp(-e * B / (k_B * T))
    if I < -J0:
        I = -J0 + 1e-20  # prevent log of negative number

    return (k_B * T / e) * np.log(I / J0 + 1)


def omic_resistance(I, R_s):
    """Calculate the voltage drop across an ohmic resistor.

    Parameters:
    I (float): Current in A.
    R_s (float): Series resistance in Ohms.

    Returns:
    float: Voltage drop in Volts.
    """
    return I * R_s


def bisection_series(V_total, inv_funcs, I_min, I_max, tol=1e-18, max_iter=200000):
    """
    Solve for current I in a series circuit using the bisection method.
    
    Parameters:
    -----------
    V_total : float
        Total applied voltage
    inv_funcs : list of callables
        Each function f_inv(I) returns the voltage across a device for current I
    I_min, I_max : floats
        Lower and upper bounds for current (must bracket the root)
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum number of iterations
    
    Returns:
    --------
    I : float
        Solved current
    V_drops : list of floats
        Voltage across each device
    """

    def g(I):
        return sum(f(I) for f in inv_funcs) - V_total

    g_min, g_max = g(I_min), g(I_max)
    if g_min * g_max > 0:
        raise ValueError("Function does not change sign between I_min and I_max — adjust bounds!")

    for n in range(max_iter):
        I_mid = 0.5 * (I_min + I_max)
        g_mid = g(I_mid)

        if abs(g_mid) < tol or (I_max - I_min) < tol:
            V_drops = [f(I_mid) for f in inv_funcs]
            return I_mid, V_drops

        # Choose the subinterval that contains the root
        if g_min * g_mid < 0:
            I_max = I_mid
            g_max = g_mid
        else:
            I_min = I_mid
            g_min = g_mid

    raise RuntimeError("Bisection did not converge")

def Back_to_back_schottky_diode_current(V, B_1, B_2, T, A_star, R=0): 
    """Calculate the current through back-to-back Schottky diodes.

    Parameters:
    V (float): Applied voltage in Volts.
    B_1 (float): Schottky barrier height of first diode in eV.
    B_2 (float): Schottky barrier height of second diode in eV.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.
    R (float): Series resistance in Ohms.

    Returns:
    float: Current in A.
    """
    if R == 0:
        inv_funcs = [
            lambda I: schottky_diode_voltage_ideal(I, B_1, T, A_star),
            lambda I: -schottky_diode_voltage_ideal(-I, B_2, T, A_star)
        ]
    else:
        inv_funcs = [
            lambda I: schottky_diode_voltage_ideal(I, B_1, T, A_star),
            lambda I: omic_resistance(I, R_s=R),
            lambda I: -schottky_diode_voltage_ideal(-I, B_2, T, A_star)
        ]
    
    I, V_drops = bisection_series(V, inv_funcs, I_min=-1e-4, I_max=1e-4)
    return I, V_drops

def Back_to_back_schottky_diode_current_eff_barrier(V, B_1, B_2, n_1, n_2, R=0, T=300, A_star=120e4): 
    """Calculate the current through back-to-back Schottky diodes.

    Parameters:
    V (float): Applied voltage in Volts.
    B_1 (float): Schottky barrier height of first diode in eV.
    B_2 (float): Schottky barrier height of second diode in eV.
    n_1 (float): Ideality factor of first diode.
    n_2 (float): Ideality factor of second diode.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.

    Returns:
    float: Current in A.
    list: Voltage drops across every device.
    """
    if R == 0:
        inv_funcs = [
            lambda I: schottky_diode_voltage_eff_B(I, B_1, T, A_star, n_1),
            lambda I: -schottky_diode_voltage_eff_B(-I, B_2, T, A_star, n_2)
        ]
    else:
        inv_funcs = [
            lambda I: schottky_diode_voltage_eff_B(I, B_1, T, A_star, n_1),
            lambda I: omic_resistance(I, R_s=R), 
            lambda I: -schottky_diode_voltage_eff_B(-I, B_2, T, A_star, n_2)
        ]
    
    I, V_drops = bisection_series(V, inv_funcs, I_min=-1e-3, I_max=1e-3)


    return I , V_drops

def Back_to_back_schottky_diode_current_eff_barrier_for_fit(V, B_1, B_2, n_1, n_2, R=0, T=300, A_star=120e4):
    """Calculate the current through back-to-back Schottky diodes, but only return the current. (Needed for the fit)

    Parameters:
    V (float): Applied voltage in Volts.
    B_1 (float): Schottky barrier height of first diode in eV.
    B_2 (float): Schottky barrier height of second diode in eV.
    n_1 (float): Ideality factor of first diode.
    n_2 (float): Ideality factor of second diode.
    T (float): Temperature in Kelvin.
    A_star (float): Richardson constant in A/m^2/K^2.

    Returns:
    float: Current in A.
    """

    result = [Back_to_back_schottky_diode_current_eff_barrier(v, B_1, B_2, n_1, n_2, R, T, A_star) for v in V]

    I = [res[0] for res in result]
    V_drops = [res[1] for res in result]
    return np.array(I)

def plot_example_back_to_back(eff=False):
    """Plot a example of the model.
    
    Parameters:
    eff (bool): Wether to use barrier hight lowering
    """
    import matplotlib.pyplot as plt

    V = np.linspace(-0.5, 0.5, 100)
    if eff:
        results = [Back_to_back_schottky_diode_current_eff_barrier(v, B_1=0.2, B_2=0.2, T=300, A_star=2.5e-9, n_1=1.2, n_2=1.2, R=40000) for v in V]
    else:
        results = [Back_to_back_schottky_diode_current(v, B_1=0.2, B_2=0.2, T=300, A_star=2.5e-9, R=40000) for v in V]

    J = [res[0] for res in results]
    V_1, V_2, V_3 = zip(*[res[1] for res in results])


    plt.figure(figsize=(5,4))
    plt.plot(V, np.array(J)*1e6)
    plt.xlabel('Voltage (V)')
    plt.ylabel(r'Current ($\mu$A)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    plt.close()
    plt.figure(figsize=(5,4))
    plt.plot(V, V_1, label='V_1')
    plt.plot(V, V_2, label='V_2')
    plt.plot(V, V_3, label='V_3')
    plt.plot(V, list(map(sum, zip(V_1, V_2, V_3))), label='V_total', color='black')
    plt.xlabel('Total Voltage (V)')
    plt.ylabel('Voltage Drop (V)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_example_schottky():
    """Plot a example of the characteristic curve of a Schottky diode with barrier hight lowering.
    """
    import matplotlib.pyplot as plt

    V = np.linspace(-1, 0.1, 100)
    J = [schottky_diode_current_eff_B(v, B=0.4, T=300, A_star=2.5e-9, n=1.3) for v in V]

    plt.plot(V, np.array(J)*1e6)
    plt.xlabel('Voltage (V)')
    plt.ylabel(r'Current Density ($\mu$A)')
    plt.title('Schottky Diode I-V Characteristic')
    plt.grid(True)
    plt.show()
    # plt.close()


def plot_inverse_example_schottky():
    """Plot a example of the inverse characteristic curve of a Schottky diode with barrier hight lowering.
    """
    import matplotlib.pyplot as plt

    I = np.linspace(-1e-4, 1e-3, 100)
    V = [schottky_diode_voltage_eff_B(i, B=0.1, T=300, n=1.1, A_star=2.5e-9) for i in I]

    plt.plot(I, V)
    plt.xlabel('Current (A)')
    plt.ylabel('Voltage (V)')
    plt.title('Inverse Schottky Diode I-V Characteristic')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    # plot_example_schottky()
    # plot_inverse_example_schottky()
    plot_example_back_to_back(True)