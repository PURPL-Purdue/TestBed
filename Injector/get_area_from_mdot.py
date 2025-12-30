import numpy as np

def get_area_from_mdot(state, mdot, p_feed, p_c, rho, Cd, gamma=None):

    if state == 'gas':

        p_cr = (2 / (gamma + 1)) ** (gamma / (gamma -1))

        if  p_c / p_feed < p_cr: # Choked condition
            area = mdot / (Cd * np.sqrt(gamma * rho * p_feed * (2 / (gamma + 1)) ** ((gamma + 1)/(gamma - 1))))

        else: # Unchoked condition
            area = mdot / (Cd * np.sqrt(2 * rho * p_feed * (gamma / (gamma-1)) * (((p_c/p_feed) ** (2/gamma)) - ((p_c/p_feed) ** ((gamma + 1)/gamma)))))

        return area

    elif state == 'liquid':

        return mdot / (Cd * np.sqrt(2 * rho * (p_feed - p_c))) # Total fuel injection area (m^2)
    
    else:
        raise ValueError("Invalid state: choose 'gas' or 'liquid'")