from pint import Quantity as Q_
from matplotlib import pyplot as plt
from dataclasses import dataclass
import math

# Some orifice height variable = thread engagement length
orifice_height = Q_(10.0, "millimeters") # Placeholder

@dataclass
class SleeveThreadDesignPoint:
    F: Q_[float] # Required force to withstand (Let FS = 2)
    sigma_UTS: Q_[float] # UTS
    tau_UTS: Q_[float] # Shear UTS
    K: float # Thread form factor for shear area

@dataclass
class SleeveThreadGeometry:
    L_eng: Q_[float] # Thread engagement length

def calcMeanDiameter(F: Q_[float], sigma_UTS: Q_[float], tau_UTS: Q_[float], K: float, L_eng: Q_[float]) -> Q_[float]:
    d_mean = (F / (math.pi * L_eng * K * tau_UTS)).to("millimeters")
    return d_mean

def calcShearArea(d_mean: Q_[float], L_eng: Q_[float], K: float) -> Q_[float]:
    A_shear = (math.pi * d_mean * L_eng * K).to("millimeters**2")
    return A_shear

def calcShearStress(F: Q_[float], A_shear: Q_[float]) -> Q_[float]:
    tau = (F / A_shear).to("MPa")
    return tau

class SleeveThreads:
    def __init__(self, design_point: SleeveThreadDesignPoint, geometry: SleeveThreadGeometry):
        self.F = design_point.F
        self.sigma_UTS = design_point.sigma_UTS
        self.tau_UTS = design_point.tau_UTS
        self.K = design_point.K
        self.L_eng = geometry.L_eng

        self.d_mean = calcMeanDiameter(self.F, self.sigma_UTS, self.tau_UTS, self.K, self.L_eng)

        self.A_shear = calcShearArea(self.d_mean, self.L_eng, self.K)
        self.tau = calcShearStress(self.F, self.A_shear)

if __name__ == "__main__":
    sleeveDesign = SleeveThreadDesignPoint(F = Q_(100, "N"), sigma_UTS = Q_(480, "MPa"), tau_UTS = Q_(240, "MPa"), K = 0.5)
    sleeveGeometry = SleeveThreadGeometry(L_eng = orifice_height)

    threads = SleeveThreads(sleeveDesign, sleeveGeometry)