from __future__ import print_function
#convert numpy array to a plain string
def plaintext(a):
    s = ""
    for x in a:
        s += "{0} ".format(x)
    return s



def Line(v0,v1):
    return "l " + plaintext(v0) + plaintext(v1) + "\n"


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


def Size(x):
    return "r {0}\n".format(float(x))


def Layer(x):
    return "y {0}\n".format(int(x))


def ArrowType(x):
    return "a {0}\n".format(int(x))


def NewPage():
    return "\n"
