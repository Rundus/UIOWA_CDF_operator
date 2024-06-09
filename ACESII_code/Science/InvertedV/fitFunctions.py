from numpy import power, cos, radians, exp, sin,pi
from ACESII_code.class_var_func import m_e,q0

def diffNFlux_for_mappedMaxwellian(x, n, T, beta, V, alpha):
    Vpara_sqrd = (2 * x * power(cos(radians(alpha)), 2) / m_e) - 2 * V / m_e + (1 - 1 / beta) * (2 * x / m_e) * (power(sin(radians(alpha)), 2))
    Vperp_sqrd = ((2 * x) / (beta * m_e)) * power(sin(radians(alpha)), 2)

    return (2 * x) * ((q0 / m_e) ** 2) * (1E2 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp((-m_e / (2 * T)) * (Vpara_sqrd + Vperp_sqrd))


def dist_Maxwellian(Vperp,Vpara,n,T):
    Emag = (0.5*m_e*(Vperp**2 + Vpara**2))/q0
    return (1E2 * n) * power(m_e / (2 * pi * q0 * T), 3 / 2) * exp(-Emag/T)

def calc_diffNFlux(Vperp,Vpara,dist):
    Emag = (0.5 * m_e * (Vperp ** 2 + Vpara ** 2)) / q0
    return (2 * Emag) * ((q0 / m_e) ** 2) * dist