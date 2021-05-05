# convert a module to a class
# コメント行はそのまま。
# 先頭の文字列もそのまま
# descもそのまま
# それ以外の部分は全部classに入れる。
# クラス変数にはself.を付ける。
# インポート行は先頭にまとめる。

import sys


def detect_header(lines):
    i = 0
    inblock = False
    for line in lines:
        if len(line) > 5 and line[:4] in ("logg",
                                          "atom",
                                          "bond",
                                          "dens",
                                          "wate",
                                          "cell",
                                          "coor",
                                          "pair",
                                          "doub",
                                          "cage"):
            return i
        i += 1
    return -1


with open(sys.argv[1], "r") as source:
    lines = source.readlines()

x = detect_header(lines)
header = lines[:x]
body = lines[x:]

newbody = ["\n",
           "class Lattice(genice.lattices.Lattice):\n",
           "    def __init__(self):\n"]
imports = ["import genice.lattices\n"]
for line in body:
    if len(line) > 5 and line[:4] in ("from", "impo"):
        imports.append(line)
    else:
        for word in (
            "waters",
            "bondlen",
            "density",
            "cell",
            "coord",
            "cages",
            "pairs",
                "fixed"):
            line = line.replace(word, "self." + word)
        # recovery
        for word in ("waters_",):
            line = line.replace("self." + word, word)
        for word, alt in {"parse_self.": "parse_",
                          "waters_and_self.": "waters_and_",
                          "estimate_self.": "estimate_",
                          "self.cellve": "cellve",
                          }.items():
            line = line.replace(word, alt)

        newbody.append(" " * 8 + line)

header += imports
for line in header + newbody:
    print(line.rstrip())
