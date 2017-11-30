
def plaintext(a):
    s = ""
    for x in a:
        s += "{0:.4f} ".format(x)
    return s



def Line(v0,v1):
    return "l " + plaintext(v0) + plaintext(v1) + "\n"


def Text(v,txt):
    return "t " + plaintext(v) + " " + txt + "\n"


def Circle(v0):
    return "c " + plaintext(v0) + "\n"


def Arrow(v0,v1):
    return "s " + plaintext(v0) + plaintext(v1) + "\n"



def Polygon(vertices):
    s = "p {0} ".format(len(vertices))
    for v in vertices:
        s += plaintext(v)
    return s + "\n"



def Color(x):
    return "@ {0}\n".format(int(x))


def SetPalette(x,R,G,B):
    return "@ {0} {1} {2} {3}\n".format(int(x),R,G,B)


def Size(x):
    return "r {0}\n".format(float(x))


def Layer(x):
    return "y {0}\n".format(int(x))


def ArrowType(x):
    return "a {0}\n".format(int(x))


def NewPage():
    return "\n"


