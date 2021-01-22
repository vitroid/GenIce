
desc={"ref": {"IV": 'Avogadro'},
      "usage": "No options available.",
      "brief": "Ice IV."
      }



import genice2.lattices
from genice2.cell import cellvectors

class Lattice(genice2.lattices.Lattice):
    def __init__(self):
        self.cell = """
        15.070097094562836 0 0
        5.129552870714296 14.170233370911644 0
        5.129552870714296 3.598426774295546 13.705722838877596
        """
        self.waters = """
        0.29274999999999984 0.2927499999999998 0.29274999999999984
        0.29274999999999984 0.29274999999999984 0.79275
        0.2927499999999997 0.7927499999999997 0.29274999999999984
        0.2927499999999998 0.7927499999999997 0.79275
        0.7927499999999997 0.2927499999999998 0.29274999999999984
        0.79275 0.29274999999999984 0.79275
        0.7927499999999996 0.7927499999999997 0.29274999999999984
        0.7927499999999998 0.7927499999999997 0.79275
        0.45724999999999993 0.45725000000000016 0.45725000000000016
        0.45725000000000016 0.45725000000000016 0.9572500000000003
        0.45725000000000016 0.9572500000000002 0.45725000000000016
        0.45725000000000016 0.9572500000000002 0.9572500000000003
        0.9572500000000002 0.45725000000000016 0.45725000000000016
        0.9572500000000004 0.45725000000000016 0.9572500000000003
        0.9572499999999999 0.9572500000000002 0.45725000000000016
        0.9572499999999999 0.9572500000000002 0.9572500000000003
        0.20725000000000016 0.20725000000000016 0.20725000000000016
        0.2072500000000001 0.2072500000000001 0.7072500000000002
        0.2072500000000001 0.70725 0.20725000000000016
        0.20725000000000005 0.70725 0.7072500000000002
        0.70725 0.20725000000000016 0.20725000000000016
        0.7072500000000002 0.2072500000000001 0.7072500000000002
        0.7072500000000002 0.70725 0.20725000000000016
        0.70725 0.70725 0.7072500000000002
        0.04274999999999983 0.042749999999999844 0.042749999999999844
        0.042749999999999816 0.042749999999999844 0.54275
        0.04274999999999985 0.5427499999999998 0.042749999999999844
        0.04274999999999979 0.5427499999999998 0.54275
        0.5427499999999997 0.042749999999999844 0.042749999999999844
        0.5427499999999997 0.042749999999999844 0.54275
        0.5427499999999997 0.5427499999999998 0.042749999999999844
        0.5427499999999998 0.5427499999999998 0.54275
        0.13020000000000032 0.19454999999999956 0.4401999999999999
        0.13020000000000043 0.1945499999999995 0.9402
        0.13020000000000032 0.6945499999999996 0.4401999999999999
        0.13020000000000032 0.6945499999999996 0.9402
        0.6302000000000003 0.19454999999999956 0.4401999999999999
        0.6302000000000004 0.1945499999999995 0.9402
        0.6302000000000003 0.6945499999999996 0.4401999999999999
        0.6302000000000003 0.6945499999999996 0.9402
        0.0554500000000004 0.1197999999999997 0.3098000000000001
        0.05545000000000047 0.11979999999999971 0.8098000000000002
        0.0554500000000004 0.6197999999999997 0.3098000000000001
        0.05545000000000053 0.6197999999999997 0.8098000000000002
        0.5554500000000004 0.1197999999999997 0.3098000000000001
        0.5554500000000004 0.11979999999999971 0.8098000000000002
        0.5554500000000003 0.6197999999999997 0.3098000000000001
        0.5554500000000003 0.6197999999999997 0.8098000000000002
        0.30545000000000044 0.059800000000000075 0.3697999999999997
        0.3054500000000004 0.05980000000000013 0.8697999999999997
        0.30545000000000033 0.5598 0.3697999999999997
        0.3054500000000005 0.5598000000000001 0.8697999999999997
        0.8054500000000003 0.059800000000000075 0.3697999999999997
        0.8054500000000002 0.05980000000000013 0.8697999999999997
        0.8054500000000001 0.5598 0.3697999999999997
        0.8054500000000002 0.5598000000000001 0.8697999999999997
        0.3802000000000002 0.19019999999999987 0.44454999999999956
        0.3802000000000004 0.19019999999999995 0.9445499999999996
        0.3802000000000003 0.6901999999999999 0.44454999999999956
        0.3802000000000003 0.6901999999999999 0.9445499999999996
        0.8802000000000001 0.19019999999999987 0.44454999999999956
        0.8802000000000003 0.19019999999999995 0.9445499999999996
        0.8802000000000001 0.6901999999999999 0.44454999999999956
        0.8802000000000003 0.6901999999999999 0.9445499999999996
        0.44019999999999987 0.1302000000000003 0.19454999999999958
        0.44020000000000004 0.1302000000000003 0.6945499999999996
        0.4402 0.6302000000000001 0.19454999999999958
        0.4401999999999999 0.6302000000000002 0.6945499999999996
        0.9401999999999997 0.1302000000000003 0.19454999999999958
        0.9401999999999997 0.1302000000000003 0.6945499999999996
        0.9401999999999999 0.6302000000000001 0.19454999999999958
        0.9402000000000001 0.6302000000000002 0.6945499999999996
        0.3098000000000001 0.055450000000000436 0.1197999999999997
        0.3098000000000001 0.05545000000000047 0.6197999999999997
        0.30979999999999996 0.5554500000000003 0.1197999999999997
        0.3098000000000001 0.5554500000000004 0.6197999999999997
        0.8098000000000001 0.055450000000000436 0.1197999999999997
        0.8097999999999999 0.05545000000000047 0.6197999999999997
        0.8097999999999999 0.5554500000000003 0.1197999999999997
        0.8098000000000001 0.5554500000000004 0.6197999999999997
        0.36979999999999963 0.3054500000000004 0.059800000000000075
        0.3697999999999997 0.30545000000000044 0.5598000000000001
        0.36979999999999963 0.8054500000000003 0.059800000000000075
        0.3697999999999996 0.8054500000000004 0.5598000000000001
        0.8697999999999996 0.3054500000000004 0.059800000000000075
        0.8697999999999996 0.30545000000000044 0.5598000000000001
        0.8697999999999996 0.8054500000000003 0.059800000000000075
        0.8697999999999996 0.8054500000000004 0.5598000000000001
        0.4445499999999995 0.38020000000000026 0.19019999999999992
        0.44454999999999945 0.3802000000000002 0.6902
        0.4445499999999994 0.8802000000000003 0.19019999999999992
        0.44454999999999945 0.8802000000000002 0.6902
        0.9445499999999994 0.38020000000000026 0.19019999999999992
        0.9445499999999994 0.3802000000000002 0.6902
        0.9445499999999994 0.8802000000000003 0.19019999999999992
        0.9445499999999997 0.8802000000000002 0.6902
        0.19454999999999956 0.44019999999999987 0.13020000000000032
        0.19454999999999956 0.4401999999999998 0.6302000000000003
        0.19454999999999956 0.9401999999999999 0.13020000000000032
        0.1945499999999995 0.9401999999999999 0.6302000000000003
        0.6945499999999996 0.44019999999999987 0.13020000000000032
        0.6945499999999994 0.4401999999999998 0.6302000000000003
        0.6945499999999994 0.9401999999999999 0.13020000000000032
        0.6945499999999993 0.9401999999999999 0.6302000000000003
        0.11979999999999963 0.3098000000000001 0.05545000000000045
        0.11979999999999963 0.3098 0.5554500000000006
        0.11979999999999963 0.8098000000000001 0.05545000000000045
        0.11979999999999952 0.8098000000000001 0.5554500000000006
        0.6197999999999998 0.3098000000000001 0.05545000000000045
        0.6197999999999997 0.3098 0.5554500000000006
        0.6197999999999996 0.8098000000000001 0.05545000000000045
        0.6197999999999995 0.8098000000000001 0.5554500000000006
        0.05980000000000009 0.36979999999999963 0.30545000000000044
        0.05980000000000016 0.3697999999999996 0.8054500000000004
        0.05980000000000009 0.8697999999999997 0.30545000000000044
        0.05980000000000005 0.8697999999999996 0.8054500000000004
        0.5598000000000001 0.36979999999999963 0.30545000000000044
        0.5598000000000001 0.3697999999999996 0.8054500000000004
        0.5598 0.8697999999999997 0.30545000000000044
        0.5597999999999999 0.8697999999999996 0.8054500000000004
        0.19019999999999992 0.44454999999999956 0.3802000000000003
        0.19019999999999992 0.44454999999999956 0.8802000000000003
        0.19019999999999992 0.9445499999999996 0.3802000000000003
        0.1901999999999998 0.9445499999999996 0.8802000000000003
        0.6901999999999997 0.44454999999999956 0.3802000000000003
        0.6901999999999999 0.44454999999999956 0.8802000000000003
        0.6901999999999999 0.9445499999999996 0.3802000000000003
        0.6901999999999998 0.9445499999999996 0.8802000000000003
        """
        self.coord = "relative"
        self.bondlen = 3
        self.density = 1.3072141048893433

        self.cell = cellvectors(a=15.070097094562836,
                           b=15.070097094562836,
                           c=15.070097094562838,
                           A=70.1,
                           B=70.1,
                           C=70.1)
