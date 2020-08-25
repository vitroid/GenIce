# coding: utf-8

desc={"ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
      "brief": "Directed graph of HBs.",
      "usage": "No options available."
      }


from logging import getLogger
import genice.formats

class Format(genice.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {4:self.hook4}


    def hook4(self, lattice):
        logger = getLogger()
        logger.info("Hook4: Output the hydrogen bond network.")

        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(lattice.reppositions))
        for i,j,k in lattice.spacegraph.edges(data=True):
            s += "{0} {1}\n".format(i,j)
        s += "-1 -1\n"
        s = "\n".join(lattice.doc) + "\n" + s
        print(s,end="")
        logger.info("Hook4: end.")
