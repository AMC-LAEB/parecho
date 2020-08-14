#!/usr/bin/env python

# Phylogenetic placement of new sequences to PhyCLIP clustered reference tree
# Authors: Edyth Parker and Alvin X. Han

from scipy.stats import levene
import re
import os
import subprocess
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
import json
import multiprocessing as mp
import itertools
import ete3

# parse phyclip tree output (NEXUS format)
def parse_phyclip_output(filename):

    id_to_taxon = {}
    cluster_to_taxon = {}

    fhandle = open(filename, "r").readlines()
    for line in fhandle:

        # parse for id, taxon and cluster
        try:
            id, taxon = re.search("^\s+(\d+)\s+([^,;]+)[,;]*$", line).group(1, 2)
            taxon = re.sub("(^'|'$|\*)", "", taxon)

            try:
                cluster = re.search("_cluster([\d\.a-z]+)", taxon).group(1)
                taxon = re.sub("_cluster[\d\.a-z]+", "", taxon)
            except:
                cluster = "unclustered"

            taxon = taxon.strip()
            id_to_taxon[id] = taxon
            try:
                cluster_to_taxon[cluster].append(taxon)
            except:
                cluster_to_taxon[cluster] = [taxon]
        except:
            pass

        # parse tree
        try:
            tree = re.search("tree[^(]+(\([^;]+;)", line).group(1)
            tree = re.sub("\[[^\]]+\]", "", tree)
        except:
            pass

    # invalid file
    if len(id_to_taxon) == 0:
        sys.exit(1)

    # replace id with taxon name
    new_tree = []
    prev_end = 0
    tree = re.sub("'", "", tree)
    for expr in re.finditer("[(,](\d+):", tree):
        new_tree.append(tree[prev_end:expr.start()+1])

        id = expr.group(1)
        new_tree.append(id_to_taxon[id])

        prev_end = expr.end()-1
    new_tree.append(tree[prev_end:])

    return "".join(new_tree), cluster_to_taxon

# searches clade tree strings of clusters
def cluster_to_clade_tree(ref_tree, c_to_t):

    eteTree = ete3.Tree(ref_tree)
    # resolve polytomy
    eteTree.resolve_polytomy()
    eteTree.ladderize()

    data = {"CLUSTER":["REF"], "TRUNK":[True],
            "TSTRING":[eteTree.write(format=5)], "REFSEQ":[eteTree.get_leaf_names()]}

    for n, node in enumerate(eteTree.traverse(strategy='levelorder')):
        n = str(n)
        if n in c_to_t.keys():
            ref_taxa = c_to_t[n]
            leaves = node.get_leaf_names()

            unclustered_taxa_subtended = list(set(c_to_t["unclustered"])&set(leaves))
            if len(unclustered_taxa_subtended) > 0:
                ref_taxa = list(set(ref_taxa)|set(unclustered_taxa_subtended))

            clade_tree_string = node.write(format=5)

            # terminal clade tree
            if set(leaves) != set(ref_taxa):
                trunk_clade = False

            # trunk clusters
            else:
                cTree = ete3.Tree(clade_tree_string)
                cTree.prune(ref_taxa)
                clade_tree_string = cTree.write(format=5)
                trunk_clade = True

            data["CLUSTER"].append(n)
            data["TRUNK"].append(trunk_clade)
            data["TSTRING"].append(clade_tree_string)
            data["REFSEQ"].append(ref_taxa)

    return pd.DataFrame(data), eteTree.get_leaf_names()

def parse_aln(filename):

    data = {}
    fhandle = open(filename, "r").readlines()

    for key, group in itertools.groupby(fhandle, key=lambda _: re.search("^>", _)):
        if key:
            header = re.sub("^>", "", next(group)).strip()
            try:
                data["HEADER"].append(header)
            except:
                data["HEADER"] = [header]
        else:
            sequence = "".join([line.strip() for line in list(group)])
            try:
                data["SEQUENCE"].append(sequence)
            except:
                data["SEQUENCE"] = [sequence]

    return pd.DataFrame(data)

def weighted_high_median(a, wts):
    N = len(a)
    wtotal = 0
    wdiscardedlow = 0

    for i in range(N):
        wtotal += wts[i]

    nn = N
    while True:
        assert (nn > 0 and len(a) == nn)

        trial = sorted(a)[int(nn/2)]
        # Count up the weight to the left of and at the trial point.
        # Weight to the right of it isn't needed
        wleft = wtrial = 0
        for i in range(nn):
            if a[i] < trial:
                wleft += wts[i]
            elif a[i] == trial:
                wtrial += wts[i]

        if 2*(wdiscardedlow + wleft) > wtotal:
            # Trial value is too high
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] < trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates

        elif 2*(wdiscardedlow + wleft + wtrial) > wtotal:
            # Trial value is just right
            return trial

        else:
            # Trial value is too low
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] > trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates
            wdiscardedlow += wleft+wtrial

        a=a[:nn]
        wts=wts[:nn]

def qn(data):
    # sort data
    data = np.sort(data)

    n = len(data)
    h = int(n/2) + 1
    k = int(h*(h-1)/2)

    left = np.arange(n+1,1,-1)
    right = np.full(n,n, dtype= np.int64)

    work = np.zeros(n) # dtype = np.float64
    weight = np.zeros(n, np.int64)
    P = np.zeros(n, np.int64)
    Q = np.zeros(n, np.int64)

    jhelp = int((n*(n+1))/2)
    knew = k+jhelp
    nL = jhelp
    nR = n*n
    found = False
    Qn = 0*data[0]

    while (nR-nL) > n:
        j = 1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                weight[j-1] = right[i] - left[i] + 1
                jhelp = left[i] + int(weight[j-1]/2)
                work[j-1] =data[i] - data[n+1-jhelp-1]
                j += 1

        trial = weighted_high_median(work[:j-1], weight[:j-1])

        j=0
        for i in range(n-1, -1, -1):
            while (j < n) and (data[i]-data[n-j-1] < trial):
                j += 1
            P[i] = j

        j = n+1
        for i in range(n):
            while data[i]-data[n-j+2-1] > trial: # 55
                j -= 1
            Q[i] = j

        sumP = sum(P)
        sumQ = sum(Q)-n # 60

        if knew <= sumP:
            right[:] = P[:]
            nR = sumP
        elif knew > sumQ:
            left[:] = Q[:]
            nL = sumQ
        else:
            Qn = trial
            found = True
            break

    if found == False:
        j=1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                for jj in range(left[i], right[i]+1, 1):
                    work[j-1] = data[i]-data[n-jj]
                    j += 1

        Qn = sorted(work[:j])[knew-nL-1]

    if n<10:
        nscale = [0, .399, .994, .512, .844, .611, .857, .669, .872][n-1]
    elif n%2 == 1:
        nscale = n/(n+1.4)
    else:
        nscale = n/(n+3.8)

    Qn = Qn*2.2219*nscale

    return Qn

class PhylogeneticPlacement():

    def __init__(self, fdat, lwr, ncpu, ezd, qval_thres, query_fdat, cluster_to_cladeTree, cluster_to_taxa, raxml_bin, submodel):
        self.fdat = fdat
        self.lwr = lwr
        self.ncpu = ncpu
        self.ezd = ezd
        self.qval_thres = qval_thres
        self.exponent = np.int(np.abs(np.floor(np.log10(ezd))))
        self.query_sequences = query_fdat
        self.cluster_to_cladeTree = cluster_to_cladeTree
        self.cluster_to_taxa = cluster_to_taxa
        self.taxon_to_cluster = {taxon:cluster for cluster, taxa in cluster_to_taxa.items() for taxon in taxa}
        self.raxml_bin = raxml_bin
        self.model = submodel

    def epa(self, cluster_to_analyse, xpercentile_edges):

        clade_info = self.cluster_to_cladeTree[self.cluster_to_cladeTree.CLUSTER == cluster_to_analyse]

        # write sequence alignment file
        reference_taxa = clade_info.iloc[0]["REFSEQ"]
        reference_sequences = self.fdat[self.fdat.HEADER.isin(reference_taxa)]
        with open("temp.fdat.fa", "w") as output:
            for i, row in reference_sequences.iterrows():
                output.write(">{}\n{}\n".format(row["HEADER"], row["SEQUENCE"]))

            for i, row in self.query_sequences.iterrows():
                output.write(">{}\n{}\n".format(row["HEADER"], row["SEQUENCE"]))

        # write tree file
        cladeTree = clade_info.iloc[0]["TSTRING"]
        with open("temp.ref_tree.tre", "w") as output:
            output.write(cladeTree)

        # run epa, slow-version across the top xpercentile_edges
        #print (xpercentile_edges)
        cmd = [self.raxml_bin, "-f", "v", "-G", str(xpercentile_edges), "-s", "temp.fdat.fa", "-t", "temp.ref_tree.tre", "-m", self.model, "-n", "results"]
        if re.search("pthread", self.raxml_bin, re.I): # threaded versions of raxml used
            cmd += ["-T", str(self.ncpu)]
        if os.path.isfile("./phyloplace_epa.log"):
            write_state = "a"
        else:
            write_state = "w"
        with open("phyloplace_epa.log", write_state) as logfile:
            subprocess.call(cmd, stdout=logfile)

        # remove all temporary files
        subprocess.call("rm temp.*", shell=True)

        return

    def analyse_jplace(self):

        with open("RAxML_portableTree.results.jplace") as fhandle:
            j_dat = json.load(fhandle)

        # read tree
        jtree = j_dat["tree"]

        # get nearest neighbouring leaf of edge number
        leaf_to_edge = {}
        new_jtree = []
        prev_end = 0
        max_edge_number = -1

        for expr in re.finditer("([^\)(,]+|\))(:[^:]+){(\d+)}", jtree):

            new_jtree.append(jtree[prev_end:expr.start()])

            edge_number = int(expr.group(3))
            if edge_number > max_edge_number:
                max_edge_number = int(edge_number)

            if expr.group(1) == ")":
                # node edge
                new_jtree.append("){}".format(edge_number))
            else:
                # leaf edge
                leaf_name = expr.group(1)
                new_jtree.append("{}".format(leaf_name))
                leaf_to_edge[leaf_name] = edge_number

            new_jtree.append(expr.group(2))
            prev_end = expr.end()

        new_jtree.append(jtree[prev_end:])

        # read clade tree as ete3 tree object
        eteCladeTree = ete3.Tree("".join(new_jtree), format=1)

        # assign edge numbers to reference clusters
        edge_to_ref_cluster = {}
        edge_to_n = {}
        self.leaf_to_n = {}
        self.n_to_node = {}
        self.node_to_n = {}
        self.n_to_edge_length = {}
        self.n_to_parent = {}
        cluster_to_n_list = {}
        desc_n_to_cluster = {}

        for n, node in enumerate(eteCladeTree.traverse(strategy="levelorder")):
            self.n_to_node[n] = node  # save ete node information
            self.node_to_n[node] = n

            if n != 0:
                parent = node.up
                self.n_to_parent[n] = self.node_to_n[parent]
                d = np.float64(np.around(node.get_distance(parent), decimals=self.exponent))
                if d <= self.ezd:
                    self.n_to_edge_length[n] = np.float64(0.)
                else:
                    self.n_to_edge_length[n] = d
            else:
                self.n_to_parent[n] = -1

            if node.is_leaf(): # leaf node
                leafname = node.name
                edge_number = leaf_to_edge[leafname]
                cluster = self.taxon_to_cluster[leafname]
                self.leaf_to_n[leafname] = n

                edge_to_n[int(edge_number)] = n
                if cluster == "unclustered":
                    edge_to_ref_cluster[int(edge_number)] = leafname # on edge leading to unclustered leaf
                else:
                    edge_to_ref_cluster[int(edge_number)] = cluster
            else:

                edge_number = node.name

                if edge_number == "": # root
                    continue
                edge_to_n[int(edge_number)] = n

                for cluster, taxa in self.cluster_to_taxa.items():
                    if cluster == "unclustered":
                        if (set(taxa)&set(node.get_leaf_names())) == set(node.get_leaf_names()): # node wholly subtends an unclustered set of taxa
                            edge_to_ref_cluster[edge_number] = "unclustered"
                        continue

                    if set(taxa) <= set(node. get_leaf_names()): # node subtends all taxa in cluster
                        try:
                            cluster_to_n_list[cluster].append(n)
                        except:
                            cluster_to_n_list[cluster] = [n]

                    else: # leaves of node subtend more than just taxa in the cluster

                        # leaves subtended by node minus unclustered sequences
                        leaves = set(node.get_leaf_names())-set(self.cluster_to_taxa["unclustered"])
                        # remaining leaves must all be in cluster taxa
                        if leaves < set(taxa):
                            try:
                                desc_n_to_cluster[n] = cluster
                            except:
                                desc_n_to_cluster[n] = cluster

        # get most desc node that subtends cluster (ie cluster root)
        n_to_cluster = {max(n_list):cluster for cluster, n_list in cluster_to_n_list.items()}

        self.cluster_to_taxon_pairdist = {}  # cluster to leaf pair distance distribution
        self.max_wcl = -1
        for cluster, taxa in self.cluster_to_taxa.items():
            if cluster == "unclustered":
                continue

            taxon_pair_list = list(itertools.combinations(taxa, 2))
            self.cluster_to_taxon_pairdist[cluster] = np.zeros(len(taxon_pair_list), dtype=np.float64)

            for _, (i, j) in enumerate(taxon_pair_list):
                self.cluster_to_taxon_pairdist[cluster][_] = self._get_distance(self.leaf_to_n[i], self.leaf_to_n[j])

            mean_x = np.mean(self.cluster_to_taxon_pairdist[cluster])
            if mean_x > self.max_wcl:
                self.max_wcl = mean_x

        # sort edge numbers of internal branchs
        for edge, n in edge_to_n.items():
            if edge in edge_to_ref_cluster: # ref cluster already found
                continue
            # n is root of cluster
            try:
                ref_cluster = n_to_cluster[n]
                edge_to_ref_cluster[edge] = ref_cluster
                continue
            except:
                pass

            # n is a descendant of a cluster
            try:
                ref_cluster = desc_n_to_cluster[n]
                edge_to_ref_cluster[edge] = "{}".format(ref_cluster)
                continue
            except:
                pass

            edge_to_ref_cluster[edge] = "trunk"

        # fields
        fields = j_dat["fields"]

        # get placements as dataframe
        placementsdf = pd.DataFrame(j_dat["placements"])

        # consolidate most likely placements
        mldat = {}
        for r, row in placementsdf.iterrows():
            query = row.n[0]
            try:
                mldat['query'].append(query)
            except:
                mldat['query'] = [query]

            placements = row.p
            p_dat = {}
            for dat in placements:
                for i, item in enumerate(dat):
                    field = fields[i]
                    try:
                        p_dat[field].append(item)
                    except:
                        p_dat[field] = [item]

            p_datdf = pd.DataFrame(p_dat)

            # get entry with largest value of lwr
            maxlwr_entry = p_datdf.iloc[p_datdf["like_weight_ratio"].idxmax()]
            try:
                mldat['lwr'].append(maxlwr_entry.like_weight_ratio)
            except:
                mldat['lwr'] = [maxlwr_entry.like_weight_ratio]

            edge = int(maxlwr_entry.edge_num)
            try:
                mldat['edge_n'].append(edge_to_n[edge])
            except:
                mldat['edge_n'] = [edge_to_n[edge]]

            try:
                mldat['putative_cluster'].append(edge_to_ref_cluster[edge])
            except:
                mldat['putative_cluster'] = [edge_to_ref_cluster[edge]]

            # distal length
            dl = np.around(maxlwr_entry.distal_length, decimals=self.exponent)
            if dl <= self.ezd:
                dl = 0.

            # pendant length
            pl = np.around(maxlwr_entry.pendant_length, decimals=self.exponent)
            if pl <= self.ezd:
                pl = 0.

            try:
                mldat['add_edge_len'].append(pl - dl)
            except:
                mldat['add_edge_len'] = [pl - dl]

        mldat = pd.DataFrame(mldat)
        self.query_to_assigned_cluster = {'query':[], 'assigned_cluster':[], 'lwr':[]}
        # get all unique putative clusters to be appended
        for pc in mldat.putative_cluster.unique():
            fil_mldat = mldat[mldat["putative_cluster"] == pc] # mldat filtered for putative cluster

            if pc == "trunk": # trunk
                # check if there are queries placed on the same branch
                for edge_n in fil_mldat.edge_n.unique():
                    # get nearest cluster to edge
                    clustern_to_d = {n:self._get_distance(n, edge_n) for n in n_to_cluster.keys()}
                    closest_putative_cluster = n_to_cluster[min(clustern_to_d, key=clustern_to_d.get)]

                    curr_dat = fil_mldat[fil_mldat.edge_n == edge_n]

                    #if len(curr_dat) == 1:
                    self.compare_against_putative_cluster(curr_dat, closest_putative_cluster, 1)

            elif pc in self.taxon_to_cluster: # outliers
                 # check if there are queries placed on the same branch
                for edge_n in fil_mldat.edge_n.unique():
                    # get nearest cluster to outlier edge
                    clustern_to_d = {n:self._get_distance(n, self.leaf_to_n[pc]) for n in n_to_cluster.keys()}
                    closest_putative_cluster = n_to_cluster[min(clustern_to_d, key=clustern_to_d.get)]

                    curr_dat = fil_mldat[fil_mldat.edge_n == edge_n]

                    #if len(curr_dat) == 1:
                    self.compare_against_putative_cluster(curr_dat, closest_putative_cluster, pc)

            else: # reference clusters
                self.compare_against_putative_cluster(fil_mldat, pc)

        #query_to_qval = self.multiple_testing_correction(query_to_pval)
        #potential_outliers_from_ref_clusters = [query for query, qval in query_to_qval.items() if qval < self.qval_thres]

        return

    def compare_against_putative_cluster(self, df, put_clus, dispersion_or_outlier_edge = 0):
        taxon_pair_dist_array = self.cluster_to_taxon_pairdist[put_clus]
        np.mean(taxon_pair_dist_array)

        for _, row in df.iterrows():

            self.query_to_assigned_cluster['query'].append(row.query)
            self.query_to_assigned_cluster['lwr'].append(row.lwr)

            add_edge_dist_array = np.zeros(len(self.cluster_to_taxa[put_clus]), dtype=np.float64)
            for _t, taxon in enumerate(self.cluster_to_taxa[put_clus]):
                add_edge_dist_array[_t] = self._get_distance(row.edge_n, self.leaf_to_n[taxon]) + row.add_edge_len

            if isinstance(dispersion_or_outlier_edge, int):
                if np.mean(np.concatenate((add_edge_dist_array, taxon_pair_dist_array))) <= self.max_wcl:
                    if dispersion_or_outlier_edge == 1: # placed on trunk/outlier edge
                        if levene(taxon_pair_dist_array, add_edge_dist_array).pvalue > self.qval_thres:
                            assigned_cluster = [put_clus] # assigned to nearest cluster
                        else:
                            assigned_cluster = [put_clus, 'outlier'] # outlier of nearest cluster
                    else:
                        assigned_cluster = [put_clus] # assigned to cluster and <= max_wcl
                else:
                    assigned_cluster = [put_clus, 'outlier'] # assigned to cluster and > max_wcl

            else: # placed on an outlier edge
                # check if we can place both outlier and query in nearest putative cluster
                outlier_edge_dist_array = np.zeros(len(self.cluster_to_taxa[put_clus]), dtype=np.float64)
                for _t, taxon in enumerate(self.cluster_to_taxa[put_clus]):
                    outlier_edge_dist_array[_t] = self._get_distance(self.leaf_to_n[dispersion_or_outlier_edge], self.leaf_to_n[taxon])

                if np.mean(np.concatenate((add_edge_dist_array, outlier_edge_dist_array, taxon_pair_dist_array))) <= self.max_wcl:
                    assigned_cluster = [put_clus, dispersion_or_outlier_edge]  # assigned to cluster and <= max_wcl
                else:
                    assigned_cluster = [put_clus, "outlier"]

            self.query_to_assigned_cluster["assigned_cluster"].append("_".join(assigned_cluster))

        return

    def multiple_testing_correction(self, pval_dict):
        dictkeys = pval_dict.keys()
        qval_list = sm.stats.multipletests([pval_dict[key] for key in dictkeys], method='fdr_bh')[1].tolist()
        return {q:qval_list[_] for _, q in enumerate(dictkeys)}

    def _get_distance(self, n_i, n_j):

        dist = 0
        mrca = self._mrca(n_i, n_j)

        while n_i != mrca:
            dist += self.n_to_edge_length[n_i]
            n_i = self.n_to_parent[n_i]

        while n_j != mrca:
            dist += self.n_to_edge_length[n_j]
            n_j = self.n_to_parent[n_j]

        return dist

    def _mrca(self, _i, _j):

        visited = []
        while _i > -1:
            visited.append(_i)
            _i = self.n_to_parent[_i]

        mrca = _j
        while mrca not in visited:
            mrca = self.n_to_parent[mrca]

        return mrca

    def analyse(self):

        self.epa("REF", 0.2) # analyse reference tree first
        self.analyse_jplace()

        return pd.DataFrame(self.query_to_assigned_cluster)

def main():
    # check number of cpus
    ncpu = mp.cpu_count()

    # parse arguments
    parser = argparse.ArgumentParser(description="Phylogenetic placement of new sequences to PhyCLIP-clustered reference tree.")
    parser.add_argument("-t", "--tree", required=True, type=str, help="PhyCLIP's output NEXUS-format tree file of reference sequences.")
    parser.add_argument("-a", '--aln', required=True, type=str, help="FASTA sequence alignment of both reference and query sequences.")
    parser.add_argument("-c", "--cpu", default=ncpu, type=str, help="Number of threads for RAxML (default = %(default)s)")
    parser.add_argument("-l", "--lwr", default=0.7, type=float, help="Minimum likelihood weight ratio threshold (default = %(default)s)")
    parser.add_argument("-z", "--ezd", default=1e-6, type=float, help="Equivalent zero distance (default = %(default)s)")
    parser.add_argument("-q", "--qval", default=0.05, type=float, help="Multiple-testing corrected p-value threshold (default = %(default)s)")
    parser.add_argument("-b", "--raxml_bin", default="raxmlHPC-PTHREADS-AVX2", type=str, help="Path to raxmlHPC binary (default = %(default)s)")
    params = parser.parse_args()

    # read PhyCLIP output tree file of reference sequences
    print ("Reading PhyCLIP output tree file...")
    try:
        reference_tree, cluster_to_taxa = parse_phyclip_output(params.tree)
    except:
        raise Exception("Invalid tree input.")

    # get tree information
    cluster_to_cladeTree, reference_taxa = cluster_to_clade_tree(reference_tree, cluster_to_taxa)

    # parse alignment fasta
    print ("\nParsing alignment...")
    fdat = parse_aln(params.aln)

    submodel =  "GTRGAMMA" if len(set(fdat["SEQUENCE"].iloc[0].upper())-set(list("NATGCRYSWKMBDHV-"))) == 0 else "PROTGAMMAGTR"
    print ("...%s sequences detected."%("nucleotide" if submodel == "GTRGAMMA" else "amino acid"))

    # filter for query sequences
    query_fdat = fdat[~fdat["HEADER"].isin(reference_taxa)]

    # phylogenetic placement
    print ("\nPerforming phylogenetic placement using RAxML-EPA (%s substitution model)..."%(submodel))
    po = PhylogeneticPlacement(fdat, params.lwr, params.cpu, params.ezd, params.qval, query_fdat, cluster_to_cladeTree, cluster_to_taxa, params.raxml_bin, submodel)
    query_to_assigned_cluster = po.analyse()
    query_to_assigned_cluster.to_csv('phyloplace_results.csv', header=True, index=None)

    print ("...done.\n")

if __name__ == "__main__":
    main()
