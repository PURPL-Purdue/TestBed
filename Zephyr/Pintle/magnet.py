from math import pi, sqrt

def calcB(s):
    return 0.5 * Br * ( (s + L) / sqrt(OR**2 + (s+L)**2) - s / sqrt(OR**2 + s**2) - (s + L) / sqrt(IR**2 + (s+L)**2) + s / sqrt(IR**2 + s**2) )

in_to_m = 0.0254
lb_to_kg = 0.4536
N_to_lbf = 0.2248
lbf_to_lbins = 386.09

IR = 0.875 / 2 * in_to_m
P = 245 / N_to_lbf / in_to_m / in_to_m
x = 0.25 * in_to_m
Br = 1.46 # N55
t = 0.127 * 10**-3 # 36AWG
L = 0.125 * in_to_m
resistivity = 1.724 * 10**-8

ohmperm = resistivity / pi / (t/2)**2

sleeve_ir = (0.843 + 0.100) / 2 * in_to_m

wraps = 25
n = wraps * L / t

print(f"P = {P / 1000.0:.1f}kPa")
print(f"Number of turns: {n}, {wraps} wraps")
print(f"Resistance: {ohmperm * 1000:.1f} mOhm/m")

max_or = 1 * in_to_m
steps = 18
for i in range(1, steps):
    OR = IR + i * (max_or - IR) / steps
    
    grad = ( calcB(x + 0.0001) - calcB(x) ) / 0.0001
    mur = 500 * (OR**2 - IR**2) / OR**2 + IR**2 / OR**2
    F = P * pi * (OR**2 - sleeve_ir**2) 
    I = F / mur / n / grad / (pi * OR**2) 

    print(f"OD: {2*OR / in_to_m:.3f}in, mur: {mur:.0f}, F: {F * N_to_lbf:.1f}lbf, I: {1000 * I:.0f}mA, {I**2 * ohmperm * 2*pi*OR*n:.0f}W")

