import csv
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import scipy.cluster.hierarchy as hcluster
import re
from sklearn.cluster import KMeans
from collections import Counter
import pandas as pd
style.use('ggplot')


with open("blastn_dmoj_converted.txt") as file:
    lines = []
    for line in file:
        # The rstrip method gets rid of the "\n" at the end of each line
        lines.append(line.rstrip().split(","))

for line in lines:
    line[0] = re.sub(r'\t', ' ', line[0])

m = []
for i in range(1, len(lines)):
    m.append(lines[i][0].split())

start = np.zeros([1, len(lines)-1])
end = np.zeros([1, len(lines)-1])
gene = []
bit_score = np.zeros([1, len(lines)-1])
ind = 0
scaffold = []
for l in m:
    start[0, ind] = int(l[8])
    end[0, ind] = int(l[9])
    bit_score[0, ind] = float(l[-4])
    gene.append(l[-3])
    scaffold.append(l[1])

    ind += 1

start = np.array(start[0])
end = np.array(end[0])
bit_score = np.array(bit_score[0])
scaffold = np.array(scaffold)
gene_unique = list(set(gene))

ma = np.column_stack((gene, scaffold, start, end, bit_score))

gene = np.array(gene)
d = dict()
for g in gene_unique:
    # l = []
    # for x in ma:
    #     if x[0] == g:
    #         l.append(x[1:])
    s = list(np.where(gene==g)[0])
    d[g] = [list(ma[s,1:5])]

a_s = np.array([x[1] for x in d[gene_unique[0]][0]])
a_e = np.array([x[2] for x in d[gene_unique[0]][0]])
a = np.column_stack((a_s, a_e))

thresh = 100000
thresh_bit = 100


def cluster(d):
    Table = dict()
    for d1 in d:
        Table[d1] = []
        s = np.array([x[1] for x in d[d1][0]])
        e = np.array([x[2] for x in d[d1][0]])
        scaffold = np.array([x[0] for x in d[d1][0]])
        bit = np.array([x[3] for x in d[d1][0]])
        a = np.column_stack((s, e))
        if len(d[d1][0]) == 1:
            Table[d1] = []
        else:
            clusters = hcluster.fclusterdata(a, thresh, criterion="distance")
            k = len(set(clusters))

            clf = KMeans(n_clusters=k)
            clf.fit(a)
            centroids = clf.cluster_centers_
            labels = clf.labels_

            for i in range(k):
                p = list(np.where(labels == i)[0])
                sca = scaffold[p]
                c = Counter(sca)

                for j in c:
                    count = c[j]
                    pp = list(np.where(scaffold == j)[0])
                    bit_scr = np.sum(np.array([float(x) for x in bit[pp]]))
                    if count >= 1 and bit_scr > thresh_bit:
                        sta = min(s[pp])
                        end = max(e[pp])
                        interval = np.array([sta, end])
                        Table[d1].append([j, interval, bit_scr])

    return Table

table = cluster(d)

table_l = []
for t in table:
    if table[t] != []:
        for i in range(len(table[t])):
            li = [t+'.'+str(i+1), table[t][i][0], float(table[t][i][1][0]), float(table[t][i][1][1]),
                  float(table[t][i][2])]
            table_l.append(li)
    else:
        table_l.append([t, np.nan, np.nan, np.nan, np.nan])

df = pd.DataFrame(list(table_l))
df = pd.DataFrame(list(table_l), columns=['Gene', 'Scaffold_ID', 'Begin', 'End', 'Cummulative_Bit_Score'])

df.to_csv('clustered_data.txt', sep='\t', index=False,
                 header="Gene\t\Scaffold_ID\tBegin\tEnd\tCummulative_Bit_Score")
