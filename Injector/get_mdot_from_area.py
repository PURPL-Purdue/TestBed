import numpy as np

def get_mdot_from_area(state, area, p_feed, p_c, rho, Cd, gamma=None):
    """
    Inverse of get_injection_area().
    Computes mass flow rate mdot given injector area.
    """
    
    if state == "gas":
        if gamma is None:
            raise ValueError("gamma required for gas calculations")

        p_ratio = p_c / p_feed
        p_cr = (2 / (gamma + 1)) ** (gamma / (gamma - 1))

        # --- Choked ---
        if p_ratio < p_cr:
            mdot = area * Cd * np.sqrt(
                gamma * rho * p_feed * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
            )

        # --- Unchoked ---
        else:
            term = (gamma / (gamma - 1)) * ((p_ratio ** (2 / gamma)) -
                                            (p_ratio ** ((gamma + 1) / gamma)))
            mdot = area * Cd * np.sqrt(2 * rho * p_feed * term)

        return mdot

    elif state == "liquid":
        # Bernoulli orifice equation
        mdot = area * Cd * np.sqrt(2 * rho * (p_feed - p_c))
        return mdot

    else:
        raise ValueError("Invalid state: choose 'gas' or 'liquid'")
