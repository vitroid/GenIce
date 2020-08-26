# coding: utf-8

desc={"ref": {"NGPH": "https://vitroid.github.io/@NGPH"},
      "brief": "Directed graph of HBs.",
      "usage": "No options available."
      }


from logging import getLogger
import genice2.formats

class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def hooks(self):
        return {4:self.hook4}


    def hook4(self, ice):
        logger = getLogger()
        logger.info("Hook4: Output the hydrogen bond network.")

        s = ""
        s += "@NGPH\n"
        s += "{0}\n".format(len(ice.reppositions))
        for i,j,k in ice.spacegraph.edges(data=True):
            s += "{0} {1}\n".format(i,j)
        s += "-1 -1\n"
        s = "\n".join(ice.doc) + "\n" + s
        self.output = s
        logger.info("Hook4: end.")
