import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.util import *
import math
from scipy.stats.stats import pearsonr, spearmanr
from scipy.special import comb
import multiprocessing
import networkx
import obonet
import statsmodels.stats.multitest as multitest
import dcor
from minepy import MINE

def is_constant(vec):
    first = vec[0]
    for i in range(len(vec)):
        if vec[i] != first:
            return False
    return True

def leave_one_out(pid, coding_rows, regulators, show, methods, return_dict):
    coding_noncoding_pairs = []
    mine = MINE()
    if show:
        print("Regulators list:\n\t"+str(regulators.head())+"\nLen="+str(len(regulators)))
        print(str(coding_rows.head())+"\nLen="+str(len(coding_rows)))
    for i in get_iterator(range(len(regulators)), show=show):
        row2 = regulators.iloc[i]
        reads2 = np.array(row2.values[1:],dtype=np.float32)
        rows = coding_rows[coding_rows[coding_rows.columns[0]] != row2[coding_rows.columns[0]]]
        for i in range(len(rows)):
            row1 = coding_rows.iloc[i]
            reads1 = np.array(row1.values[1:],dtype=np.float32)
            
            if "MIC" in methods:
                mine.compute_score(reads1, reads2)
                mic_corr = mine.mic()
                if(mic_corr > 0.5) or (mic_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], mic_corr, "MIC"))
            if "PRS" in methods:
                pearson_corr, p = pearsonr(reads1, reads2)
                if(pearson_corr > 0.5) or (pearson_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], pearson_corr, "PRS"))
            if "SPR" in methods:
                spearman_corr, p = spearmanr(reads1, reads2)
                if(spearman_corr > 0.5) or (spearman_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], spearman_corr, "SPR"))
            if "DC" in methods:
                dc_corr = dcor.distance_correlation(reads1, reads2)
                if(dc_corr > 0.5) or (dc_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], dc_corr, "DC"))
            
    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def try_find_coexpression_process(pid, coding_rows, nc_rows, show, methods, return_dict):
    coding_noncoding_pairs = []
    mine = MINE()
    #print("Methods="+str(methods))
    for i in get_iterator(range(len(coding_rows)), show=show):
        row1 = coding_rows.iloc[i]
        reads1 = np.array(row1.values[1:],dtype=np.float32)
        for name, row2 in nc_rows.iterrows():
            reads2 = np.array(row2.values[1:],dtype=np.float32)
            if "MIC" in methods:
                mine.compute_score(reads1, reads2)
                mic_corr = mine.mic()
                if(mic_corr > 0.5) or (mic_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], mic_corr, "MIC"))
            if "PRS" in methods:
                pearson_corr, p = pearsonr(reads1, reads2)
                if(pearson_corr > 0.5) or (pearson_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], pearson_corr, "PRS"))
            if "DC" in methods:
                dc_corr = dcor.distance_correlation(reads1, reads2)
                if(dc_corr > 0.5) or (dc_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], dc_corr, "DC"))
            if "SPR" in methods:
                spearman_corr, p = spearmanr(reads1, reads2)
                if(spearman_corr > 0.5) or (spearman_corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], spearman_corr, "SPR"))
    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def get_ancestors(graph, parent_id):
    ancestors = networkx.ancestors(graph, parent_id)
    return ancestors

def pvalue(m, N, M, n):
        comb_sum = 0
        for i in range(m, min(n,M)+1):
            comb_sum += comb(M, i, exact=True)*comb(N-M, n-i, exact=True)
        p = comb_sum / comb(N, n, exact=True)
        return p

def pvalue_process(N, possible_gene_term, genes_coexpressed_with_ncRNA,
                    genes_annotated_with_term, print_progress, pid, return_dict):
    gene_term_pvalue = []
    n_lens = [len(genes_coexpressed_with_ncRNA[possible_gene_term[i][0]]) for i in range(len(possible_gene_term))]
    M_lens = [len(genes_annotated_with_term[possible_gene_term[i][1]]) for i in range(len(possible_gene_term))]
    m_lens = [len(genes_annotated_with_term[possible_gene_term[i][1]]
                .intersection(genes_coexpressed_with_ncRNA[possible_gene_term[i][0]]))
                for i in range(len(possible_gene_term))]
    if print_progress:
        print("Printing progress for thread " + str(pid))
    for i in get_iterator(possible_gene_term,show=print_progress):
        n = n_lens[i] # number of genes in test set
        #if n >= 5:
        M = M_lens[i] # number of genes with annotation
        #    if M >= 5:
        m = m_lens[i] # number of genes in test set with annotation
        #if n_lens[i] >= 5 and M_lens[i] >= 5 and m_lens[i] >= 1:
        pval = pvalue(m, N, M, n)
        gene_term_pvalue.append((possible_gene_term[i][0],
                            possible_gene_term[i][1],pval))
    return_dict[pid] = gene_term_pvalue

def get_valid_associations(genes_coexpressed_with_ncRNA, genes_annotated_with_term, possible_gene_term,
                            min_n, min_M, min_m):
    n_lens = [len(genes_coexpressed_with_ncRNA[gene]) for gene, term in possible_gene_term]
    M_lens = [len(genes_annotated_with_term[term]) for gene, term in possible_gene_term]
    m_lens = [len(genes_annotated_with_term[term]
                .intersection(genes_coexpressed_with_ncRNA[gene]))
                for gene, term in possible_gene_term]
    valid = [n_lens[i] >= min_n and M_lens[i] >= min_M and m_lens[i] >= min_m for i in range(len(possible_gene_term))]
    return valid


def parallel_pvalues(N, possible_gene_term, valid_gene_term, genes_coexpressed_with_ncRNA,
                    genes_annotated_with_term, threads):
    print("Dividing chunks for parallel p-value calculation")
    gene_terms = []
    for i in range(len(possible_gene_term)):
        if valid_gene_term[i]:
            gene_terms.append(possible_gene_term[i])
    print("Calculating pvalues for " + str(len(gene_terms)) + " associations.")
    #gene_terms = list(filter(lambda x: x == True))
    calcs_per_thread = int(len(gene_terms) / threads)+1
    calc_chunks = list(chunks(gene_terms, calcs_per_thread))
    print("Calculating p-values")
    gene_term_pvalue = []

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    processes = []
    showing = True

    i = 0
    for chunk in calc_chunks:
        p = multiprocessing.Process(target=pvalue_process,
                args=(N, chunk, genes_coexpressed_with_ncRNA,
                    genes_annotated_with_term, showing, i, return_dict, ))
        processes.append(p)
        print("Spawned process")
        p.start()
        if showing:
            showing = False
        i += 1

    for p in processes:
        p.join()
        print("Finished process " + str(p.pid))

    print("Merging results")
    for value in return_dict.values():
        gene_term_pvalue += value
    print("Saved results from processes")
    manager.shutdown()
    return gene_term_pvalue
