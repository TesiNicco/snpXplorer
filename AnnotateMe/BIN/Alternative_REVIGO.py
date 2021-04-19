## TRY A DIFFERENT APPROACH FOR THE CLUSTERING OF GO TERMS
## MAJOR DRAWBACKS OF CURRENT METHODS ARE:

from pygosemsim import download
from pygosemsim import graph
from pygosemsim import similarity
import numpy as np
import csv
import networkx as nx
import sys

# get go database
#download.download("go-basic", "http://current.geneontology.org/ontology/go.obo")

# load go graph from downloaded file
G = graph.from_resource("/usr/local/lib/python3.9/site-packages/pygosemsim/_resources/go-basic")

# take precalculated lower bounds
similarity.precalc_lower_bounds(G)

# calculate semantic similarity -- this is an example
# similarity.lin(G, "GO:0004340", "GO:0019158")

# read file with all go terms and corresponding pvalues -- take only the significant ones
go_list = []
go_p = {}
with open(sys.argv[1]) as finp:
    for line in finp:
        line = line.rstrip().split()
        go, p = line[0], float(line[1])
        go_list.append(go)
        go_p[go] = p

# then loop
dist_mt = []                        # this will be a list of lists representative of a triangular matrix with all distances inside
for i in range(len(go_list)):
    tmp_list = []                   # temporary list to store each line's result
    term1 = go_list[i]              # this is term1
    for j in range(len(go_list)):
        term2 = go_list[j]          # this is term2
        dist = similarity.lin(G, term1, term2)
        tmp_list.append(dist)
    dist_mt.append(tmp_list)

# finally write output matrix
with open(sys.argv[2], "w") as outf:
    header = ",".join(go_list)
    outf.write(header)
    outf.write("\n")
    wr = csv.writer(outf)
    wr.writerows(dist_mt)

# then do the same using resnik
#dist_mt = []                        # this will be a list of lists representative of a triangular matrix with all distances inside
#for i in range(len(go_list)):
#  tmp_list = []                   # temporary list to store each line's result
#  term1 = go_list[i]              # this is term1
#  for j in range(len(go_list)):
#     term2 = go_list[j]          # this is term2
#     dist = similarity.resnik(G, term1, term2)
#     tmp_list.append(dist)
#  dist_mt.append(tmp_list)

# finally write output matrix
#with open(sys.argv[3], "w") as outf:
#   header = ",".join(go_list)
#   outf.write(header)
#   outf.write("\n")
#   wr = csv.writer(outf)
#   wr.writerows(dist_mt)

