import numpy as np

m_e = 9.11e-28
m_u = 1.66e-24
c = 3.e10
k = 1.38e-16
h = 6.63e-27
a = 5.67e-15
e = 4.8e-10 # esu

# crystallization
Gamma_C = 170

mu_e = 2
mu_I = 1.0

mu = 1.0 / (1.0/mu_e + 1.0/mu_I)

Z = 1.0

# Fermi momentum constant (from number density integral)
B_fermi = (8.0*np.pi/3.0)*m_u*(m_e*c/h)**3

def deg_ideal(rho, mu_e=2.0):
    """temperature on the line separating degeneracy and ideal gas as
    a function of rho"""

    # dimensionless Fermi momentum
    x_F = (rho / mu_e / B_fermi)**(1./3.)

    # set E_F = 3/2 k T
    T = ((2./3.)/k) * m_e * c**2 * (np.sqrt(1 + x_F**2) - 1)

    return T

def rad_ideal(rho, mu=0.5):
    """temperature on the line separating radiation and ideal gas as a
    function of rho

    """

    T = ((3 * k / (m_u * a)) * (rho / mu))**(1./3.)
    return T

def crystallization(rho):
    """temperature where crystallization of ions sets in"""

    T = ((rho / mu_I) * (Z**2 * e**2 / (Gamma_C * k))**3 * (4.0 * np.pi / 3.0 / m_u))**(1./3.)
    return T

