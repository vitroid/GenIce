
# Multiplicity

uniq = dict()

values = {"0":+1, "1":-1}

def isPeriodic(bits):
    s = sum([values[bit] for bit in bits])
    return s % 3 == 0


def search0(bits):
    L = len(bits)
    doubled = bits + bits
    for i in range(L):
        rotated = doubled[i:i+L]
        if rotated in uniq:
            return rotated
    return None

inv = {"0":"1", "1":"0"}

def search(bits):
    x = search0(bits)
    if x is None:
        b2 = ""
        for bit in bits:
            b2 += inv[bit]
        x = search0(b2)
        if x is None:
            b2 = ""
            for bit in bits:
                b2 = bit + b2
            x = search0(b2)
            if x is None:
                b2 = ""
                for bit in bits:
                    b2 = inv[bit] + b2
                x = search0(b2)
    return x



for i in range(2**12):
    bits = "{0:012b}".format(i)
    if isPeriodic(bits):
        x = search(bits)
        if x is None:
            uniq[bits] = 1
        else:
            uniq[x] += 1

# hexagonal: alternative stacking
# cubic: progressive stacking
# So the layer types can be determined by the sequence of two successive digits

def layertypes(bits):
    L = len(bits)
    bits = bits + bits[0:1]
    code = ""
    for i in range(L):
        digits = bits[i:i+2]
        if digits in ["00", "11"]:
            code += "c"
        else:
            code += "h"
    return code
        
for bits in uniq:
    print("{0} x {1}".format(layertypes(bits), uniq[bits]))

