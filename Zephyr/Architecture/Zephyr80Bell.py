import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from bisect import bisect_left

#Modified based on design by ___

"""
 Bell Nozzle Generator with high-precision angle display (degrees)
"""

# sp.heat, area_ratio, throat_radius, length percentage
def bell_nozzle(k, aratio, Rt, l_percent):
    # upto the nozzle designer, usually -135
    entrant_angle   = -135
    ea_radian       = math.radians(entrant_angle)

    # nozzle length percentage
    if l_percent == 60:     Lnp = 0.6
    elif l_percent == 80:   Lnp = 0.8
    elif l_percent == 90:   Lnp = 0.9    
    else:                   Lnp = 0.8

    # find wall angles (theta_n, theta_e) for given aratio (ar)
    angles = find_wall_angles(aratio, Rt, l_percent)
    nozzle_length = angles[0]; theta_n = angles[1]; theta_e = angles[2]

    data_interval = 100

    # entrant functions
    ea_start = ea_radian
    ea_end   = -math.pi/2    
    angle_list = np.linspace(ea_start, ea_end, data_interval)
    xe = []; ye = []
    for i in angle_list:
        xe.append(1.5 * Rt * math.cos(i))
        ye.append(1.5 * Rt * math.sin(i) + 2.5 * Rt)

    # exit section
    ea_start = -math.pi/2
    ea_end   = theta_n - math.pi/2
    angle_list = np.linspace(ea_start, ea_end, data_interval)    
    xe2 = []
    ye2 = []
    for i in angle_list:
        xe2.append(0.382 * Rt * math.cos(i))
        ye2.append(0.382 * Rt * math.sin(i) + 1.382 * Rt)

    # bell section
    Nx = 0.382 * Rt * math.cos(theta_n - math.pi/2)
    Ny = 0.382 * Rt * math.sin(theta_n - math.pi/2) + 1.382 * Rt 
    Ex = Lnp * ((math.sqrt(aratio) - 1) * Rt) / math.tan(math.radians(15))
    Ey = math.sqrt(aratio) * Rt 
    m1 = math.tan(theta_n); m2 = math.tan(theta_e)
    C1 = Ny - m1*Nx; C2 = Ey - m2*Ex
    Qx = (C2 - C1)/(m1 - m2)
    Qy = (m1*C2 - m2*C1)/(m1 - m2)    

    int_list = np.linspace(0, 1, data_interval)
    xbell = []
    ybell = []
    for t in int_list:        
        xbell.append((1-t)**2 * Nx + 2*(1-t)*t*Qx + t**2*Ex)
        ybell.append((1-t)**2 * Ny + 2*(1-t)*t*Qy + t**2*Ey)

    # create negative values for the other half of nozzle
    nye  = [-y for y in ye]
    nye2 = [-y for y in ye2]
    nybell = [-y for y in ybell]

    return angles, (xe, ye, nye, xe2, ye2, nye2, xbell, ybell, nybell)

# find wall angles (theta_n, theta_e) in radians for given aratio
def find_wall_angles(ar, Rt, l_percent = 80):
    aratio      = [ 4,    5,    10,   20,   30,   40,   50,   100]
    theta_n_60  = [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 35.0, 40.0]    
    theta_n_80  = [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
    theta_n_90  = [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0]
    theta_e_60  = [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2]
    theta_e_80  = [14.0, 13.0, 11.0,  9.0,  8.5,  8.0,  7.5,  7.0]
    theta_e_90  = [11.5, 10.5,  8.0,  7.0,  6.5,  6.0,  6.0,  6.0]    

    f1 = ((math.sqrt(ar) - 1) * Rt) / math.tan(math.radians(15))

    if l_percent == 60:
        theta_n = theta_n_60; theta_e = theta_e_60; Ln = 0.6 * f1
    elif l_percent == 80:
        theta_n = theta_n_80; theta_e = theta_e_80; Ln = 0.8 * f1        
    elif l_percent == 90:
        theta_n = theta_n_90; theta_e = theta_e_90; Ln = 0.9 * f1    
    else:
        theta_n = theta_n_80; theta_e = theta_e_80; Ln = 0.8 * f1

    x_index, _ = find_nearest(aratio, ar)

    if round(aratio[x_index], 1) == round(ar, 1):
        return Ln, math.radians(theta_n[x_index]), math.radians(theta_e[x_index])

    if (x_index>2):
        ar_slice = aratio[x_index-2:x_index+2]        
        tn_slice = theta_n[x_index-2:x_index+2]
        te_slice = theta_e[x_index-2:x_index+2]
    elif ((len(aratio)-x_index) <= 1):
        ar_slice = aratio[x_index-2:len(aratio)]        
        tn_slice = theta_n[x_index-2:len(aratio)]
        te_slice = theta_e[x_index-2:len(aratio)]
    else:
        ar_slice = aratio[0:x_index+2]        
        tn_slice = theta_n[0:x_index+2]
        te_slice = theta_e[0:x_index+2]

    tn_val = interpolate(ar_slice, tn_slice, ar)
    te_val = interpolate(ar_slice, te_slice, ar)                        

    return Ln, math.radians(tn_val), math.radians(te_val)

# simple linear interpolation
def interpolate(x_list, y_list, x):
    if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
        raise ValueError("x_list must be in strictly ascending order!")
    intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
    slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals]

    if x <= x_list[0]:
        return y_list[0]
    elif x >= x_list[-1]:
        return y_list[-1]
    else:
        i = bisect_left(x_list, x) - 1
        return y_list[i] + slopes[i] * (x - x_list[i])

# find the nearest index in the list
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]  

# nozzle contour plot
def plot_nozzle(ax, title, Rt, angles, contour):
    nozzle_length, theta_n, theta_e = angles
    xe, ye, nye, xe2, ye2, nye2, xbell, ybell, nybell = contour

    ax.set_aspect('equal')
    ax.plot(xe, ye, linewidth=2.5, color='g')
    ax.plot(xe, nye, linewidth=2.5, color='g')
    ax.plot(xe2, ye2, linewidth=2.5, color='r')
    ax.plot(xe2, nye2, linewidth=2.5, color='r')
    ax.plot(xbell, ybell, linewidth=2.5, color='b')
    ax.plot(xbell, nybell, linewidth=2.5, color='b')

    # draw angles
    draw_angle_arc(ax, theta_n, [xe2[-1], nye2[-1]], r'$\theta_n$')
    draw_angle_arc(ax, theta_e, [xbell[-1], nybell[-1]], r'$\theta_e$')

    ax.axhline(color='black', lw=0.5, linestyle="dashed")
    ax.axvline(color='black', lw=0.5, linestyle="dashed")
    ax.grid()
    ax.minorticks_on()
    ax.grid(which='major', linestyle='-', linewidth='0.5')
    ax.grid(which='minor', linestyle=':', linewidth='0.5')
    plt.title(title, fontsize=9)

# theta in rad, origin, degree symbol
def draw_angle_arc(ax, theta_rad, origin, degree_symbol=r'$\theta$'):
    length = 50
    startx, starty = origin
    endx = startx + np.cos(-theta_rad) * length * 0.5
    endy = starty + np.sin(-theta_rad) * length * 0.5
    ax.plot([startx, endx], [starty, endy], linewidth=0.5, color='k')
    arc_obj = Arc([startx, starty], 1, 1, angle=0, theta1=0, theta2=math.degrees(theta_rad), color='k')
    ax.add_patch(arc_obj)
    angle_deg = math.degrees(theta_rad)
    ax.text(startx + 0.5, starty + 0.5, f"{degree_symbol} = {angle_deg:.4f}°")

# ----------------  MAIN ------------------
if __name__=="__main__":
    
    #USER CONFIG
    k = 1.21
    l_percent = 80      #60, 80, or 90
    aratio = 4.897
    throat_radius = 83.3/2
    angles, contour = bell_nozzle(k, aratio, throat_radius, l_percent)
    title = f'Bell Nozzle [Area Ratio = {aratio:.4f}, Throat Radius = {throat_radius:.4f}]'
    fig, ax = plt.subplots(figsize=(8,6))
    plot_nozzle(ax, title, throat_radius, angles, contour)
    plt.show()

