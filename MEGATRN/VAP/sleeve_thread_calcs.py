from pint import Quantity as Q_
from matplotlib import pyplot as plt
from dataclasses import dataclass
import math

@dataclass
class SleeveThreadDesignPoint:
    F: Q_[float] # Required force to withstand (Let FS = 2)
    sigma = Q_[float] # UTS
    tau = Q_[float] # Shear UTS
    K = float # Thread form factor for shear area

@dataclass
class SleeveThreadGeometry:
    d_mean: Q_[float] # Mean diameter
    L_eng: Q_[float] # Thread engagement length

def calcMeanDiameter(F: Q_[float], sigma: Q_[float], tau: Q_[float], K: float, L_eng: Q_[float]):
    d_mean = F / (math.pi * L_eng * K * tau)
    return d_mean

class SleeveThreads:
    def __init__(self, design_point: SleeveThreadDesignPoint, geometry: SleeveThreadGeometry):
        self.d_mean = calcMeanDiameter(design_point.F, design_point.sigma, design_point.tau, design_point.K, geometry.L_eng)
        pass
pass