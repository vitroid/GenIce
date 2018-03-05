# parallelepiped cell


import numpy as np
import sys
import logging

def rel_wrap(relvec):
    return relvec - np.floor(relvec + 0.5)


def rel_wrapf(relvec):
    return relvec - np.floor(relvec)


class Cell():
    mat = np.zeros(9).reshape(3,3)
    inv = None
    def __init__(self, desc=None, celltype=None):
        if celltype is not None:
            self.parse(desc,celltype)
        else:
            # copy the Cell class
            self.mat = desc.mat.copy()
            self.inv = desc.inv.copy()

            
    def abs2rel(self, absvecs):
        return np.dot(absvecs, self.inv)

        
    def rel2abs(self, relvec):
        return np.dot(relvec, self.mat)

    
    def abs_wrap(self, absvec):
        return self.rel2abs(rel_wrap(self.abs2rel(absvec)))


    def abs_wrapf(self, absvec):
        return self.rel2abs(rel_wrapf(self.abs2rel(absvec)))

    
    def volume(self):
        return np.linalg.det(self.mat)
    

    def scale(self, x):
        self.mat *= x
        self.inv = np.linalg.inv(self.mat)
        

    def scale2(self, x):
        for d in range(3):
            self.mat[d, :] = self.mat[d, :] * x[d]
        self.inv = np.linalg.inv(self.mat)


    def parse(self, desc, celltype):
        logger = logging.getLogger()
        if celltype == "rect":
            if type(desc) is str:
                vec = np.fromstring(desc, sep=" ")
            elif type(desc) is list:
                vec = np.array(desc)
            logger.debug("parse_cell 1: {0}".format(self.mat))
            self.mat = np.diag(vec)
        elif celltype == "monoclinic":
            if type(desc) is str:
                desc = np.fromstring(desc, sep=" ")
            elif type(desc) is list:
                desc = np.array(desc)
            beta = desc[3] * pi / 180.
            M = np.array(((desc[0] * 1.0, desc[1] * 0.0, desc[2] * cos(beta)),
                          (desc[0] * 0.0, desc[1] * 1.0, desc[2] * 0.0),
                          (desc[0] * 0.0, desc[1] * 0.0, desc[2] * sin(beta))))
            # all the vector calculations are done in transposed manner.
            self.mat = M.transpose()
        elif celltype == "triclinic":
            """
            Put the vectors like following:
            cell = "ax 0 0 bx by 0 cx cy cz"
            when you define a unit cell in lattices/
            """
            if type(desc) is str:
                desc = np.fromstring(desc, sep=" ")
                logger.debug(desc)
            elif type(desc) is list:
                desc = np.array(desc)
            self.mat = np.reshape(desc, (3, 3))
            # assert cell[0, 1] == 0 and cell[0, 2] == 0 and cell[1, 2] == 0
        else:
            logger.error("unknown cell type: {0}".format(celltype))
            sys.exit(1)
        self.inv = np.linalg.inv(self.mat)

