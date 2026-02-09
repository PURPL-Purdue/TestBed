orings = []

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

def mixed_to_float(mixed):
    parts = mixed.split('-')
    if len(parts) == 1:
        frac = parts[0].split('/')
        if len(frac) == 1:
            return frac
        return int(frac[0]) / int(frac[1])
    else:
        num,den = parts[1].split('/')
        return float(parts[0]) + int(num) / int(den)

def inner_diam(diam):
    for ring in orings.rev():
        if ring['ID'] <= diam:
            return ring

def min_inner(diam):
    for ring in orings:
        if ring['ID'] >= diam:
            return ring

def get_gland_size():
    return {'depth': 0.112, 'width': 0.190, 'spacing': 0.250}
