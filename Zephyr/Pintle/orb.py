fittings = []
plugs = []

def mixed_to_float(mixed):
    parts = mixed.split(' ')
    if len(parts) == 1:
        frac = parts[0].split('/')
        if len(frac) == 1:
            return float(frac[0])
        return int(frac[0]) / int(frac[1])
    else:
        num,den = parts[1].split('/')
        return float(parts[0]) + int(num) / int(den)

with open("orbs.csv", "r") as file:
    lines = file.readlines()
    for line in lines:
        parts = line.split(',')
        num = int(parts[0])
        tube = mixed_to_float(parts[1])
        thread = parts[2].split('-')
        size = mixed_to_float(thread[0])
        pitch = int(thread[1])
        specs = parts[3].split(' ')
        fittings.append({
            'num': num,
            'tube': tube,
            'thread_diam': size,
            'thread_pitch': 1 / pitch,
            'spotface_diam': float(specs[0]) / 25.4,
            'min_id': float(specs[1]) / 25.4,
            'gauging': float(specs[2]) / 25.4,
            'chamfer_od': float(specs[3]) / 25.4,
            'chamfer_length': float(specs[4]) / 25.4,
            'tap_length': float(specs[5]) / 25.4,
            'spotface_depth': float(specs[6]) / 25.4,
            'thread_length': float(specs[7]) / 25.4,
            'chamfer_angle': float(specs[8]),
            'secondary_chamfer_angle': 45
        })

with open("orb_plug.csv", "r") as file:
    lines = file.readlines()
    for line in lines:
        parts = line.split(',')
        num = int(parts[0])
        thread_diam = mixed_to_float(parts[1])
        specs = parts[2].split(' ')
        allen = mixed_to_float(specs[0])
        thread_length = float(specs[1])
        total_length = float(specs[2])
        total_diam = float(specs[3])
        for fitting in fittings:
            if fitting['num'] == num:
                chamfer_od = fitting['chamfer_od']
                chamfer_length = fitting['chamfer_length']
                chamfer_angle = fitting['chamfer_angle']
                secondary_chamfer_angle = fitting['secondary_chamfer_angle']
                thread_pitch = fitting['thread_pitch']
                break

        plugs.append({
            'num': num,
            'thread_diam': thread_diam,
            'allen': allen,
            'thread_length': thread_length,
            'thread_pitch': thread_pitch,
            'total_length': total_length,
            'total_diam': total_diam,
            'chamfer_od': chamfer_od,
            'chamfer_length': chamfer_length,
            'chamfer_angle': chamfer_angle,
            'secondary_chamfer_angle': secondary_chamfer_angle,
        })
        

def orb_min_diam(diam):
    for fitting in fittings:
        if fitting['tube'] >= diam:
            return fitting
    return None

def plug_min_diam(diam):
    for plug in plugs:
        if plug['thread_diam'] >= diam:
            return plug
    return None

def dash_num(num):
    for fitting in fittings:
        if fitting['num'] == num:
            return fitting
    return None

