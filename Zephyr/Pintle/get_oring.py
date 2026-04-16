orings = []
radial_glands = []

def load_orings():
    with open("orings.csv", "r") as file:
        lines = file.readlines()
        for line in lines:
            num,ID,OD = line.split()
            orings.append({
                'num': int(num),
                'ID': mixed_to_float(ID),
                'OD': mixed_to_float(OD)
            })

    with open("oring_glands.csv", "r") as file:
        lines = file.readlines()
        i = 0
        while i < len(lines):
            line1 = lines[i]
            i += 1
            line2 = lines[i]
            i += 1

            num, max_rod_od, max_bore_id, max_gland_id = line1.split(' ')
            line2.lstrip()
            min_rod_od, min_bore_id, min_gland_id = line2.split(' ')
            radial_glands.append({
                'num': int(num),
                'rod_OD_min': float(min_rod_od),
                'rod_OD_max': float(max_rod_od),
                'bore_ID_min': float(min_bore_id),
                'bore_ID_max': float(max_bore_id),
                'gland_ID_min': float(min_gland_id),
                'gland_ID_max': float(max_gland_id)
            })

def mixed_to_float(mixed):
    parts = mixed.split('-')
    if len(parts) == 1:
        frac = parts[0].split('/')
        if len(frac) == 1:
            return float(frac[0])
        return int(frac[0]) / int(frac[1])
    else:
        num,den = parts[1].split('/')
        return float(parts[0]) + int(num) / int(den)

def inner_diam(diam):
    for ring in orings:
        if ring['ID'] >= diam:
            return ring

def radial_min_rod_od(diam):
    for gland in radial_glands:
        if diam <= gland['rod_OD_max']:
            return {'num': gland['num'], 'rod_od': gland['rod_OD_max']}



def get_face_gland_size(num):
    if num // 100 == 0:
        depth = 0.052
        width = 0.104
    elif num // 100 == 1:
        depth = 0.077
        width = 0.139
    elif num // 100 == 2:
        depth = 0.104
        width = 0.182
    elif num // 100 == 3:
        depth = 0.157
        width = 0.280
    elif num // 100 == 4:
        depth = 0.206
        width = 0.352
        
    return {'depth': depth, 'width': width, 'spacing': width}

def get_radial_gland_size(num):
    for gland in radial_glands:
        if gland['num'] == num:
            if num == 1:
                width = 0.092
            elif num == 2:
                width = 0.097
            elif num == 3:
                width = 0.107
            elif num >= 4 and num <= 7:
                width = 0.117
            elif num >= 8 and num <= 28:
                width = 0.107
            elif num >= 104 and num <= 109:
                width = 0.155
            elif num >= 110 and num <= 149:
                width = 0.145
            elif num >= 210 and num <= 247:
                width = 0.190
            elif num >= 325 and num <= 349:
                width = 0.275
            elif num >= 425 and num <= 460:
                width = 0.350
            else:
                print(f"Unsupported O-ring #{num:03d}")
                exit(1)

            return {
                    'id': gland['bore_ID_min'],
                    'depth': gland['gland_ID_min'] - gland['bore_ID_min'],
                    'width': width,
                    'spacing': width
            }
