import json
import requests
import sys

# citations = [
#     ["Baez 1998", "10.1063/1.469938", ""],
#     ["Rosso 2016", "10.1038/ncomms13394", ""],
#     ["Sikiric 2010", "10.1107/S0108767310022932", ""],
#     ["Falenty 2014", "10.1038/nature14014", ""],
#     ["Fan 2010", "", "Xiaofeng Fan, Dan Bing, Jingyun Zhang, Zexiang Shen, Jer-Lai Kuo,  Computational Materials Science 49 (2010) S170–S175."],
#     ["Fennell 2005", "10.1021/ct050005s", ""],
#     ["Hirata 2017", "10.1021/acs.langmuir.7b01764", ""],
#     ["Huang 2016", "10.1126/sciadv.1501010", ""],
#     ["Huang 2017", "10.1016/j.cplett.2017.01.035", ""],
#     ["Jeffrey 1984", "", "G. A. Jeffrey, In Inclusion Compounds; J. L. Atwood, J. E. D. Davies, D. D. MacNicol, Eds.; Academic Press: London, 1984, Vol. 1, Chap. 5.",],
#     ["Karttunen 2011", "10.1021/ic102178d", ""],
#     ["Koga 2001", "10.1038/35090532", ""],
#     ["Kosyakov 1999", "10.1007/BF02903652", ""],
#     ["Kuhs 1998", "10.1063/1.448109", ""],
#     ["Liu 2019", "10.1073/pnas.1900739116", ""],
#     ["Lobban 1998", "10.1038/34622", ""],
#     ["Londono 1988", "10.1038/332141a0", ""],
#     ["Londono 1993", "10.1063/1.464942", ""],
#     ["Matsui 2017", "/10.1063/1.4994757", ""],
#     ["Matsui 2019", "10.1063/1.5083021", ""],
#     ["Matsumoto 2019", "10.1063/1.5096556", ""],
#     ["Mochizuki 2014", "10.1039/c4cp01616e", ""],
#     ["Mousseau 2001", "10.1016/S1359-0286(02)00005-0", ""],
#     ["Nakamura 2015", "10.1021/acs.jpcb.5b09544", ""],
#     ["Russo 2014", "10.1038/NMAT3977", ""],
#     ["Salzmann 2006", "10.1126/science.1123896", ""],
#     ["Smirnov 2013", "10.1021/jz401669d", ""],
#     ["Vos 1993", "10.1103/PhysRevLett.71.3150", ""],
#     ["Yagasaki 2018", "10.1021/acs.jpcb.8b04441", ""],
# ]


citations = json.load(sys.stdin)

updated = []
for key, doi, desc in citations:
    if len(doi) > 0 and len(desc) == 0:
        style = "iso690-author-date-en"
        headers = {"Accept": f"text/x-bibliography; style={style}"}
        r = requests.get(f"https://doi.org/{doi}", headers=headers)
        if r.status_code == 200:
            desc = r.content.decode("UTF-8")
        else:
            desc = ""
    updated.append([key, doi, desc])

print(json.dumps(updated, indent=4))
