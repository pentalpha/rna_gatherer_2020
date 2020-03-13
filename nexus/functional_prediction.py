import os
import threading
import numpy as np
import pandas as pd
from tqdm import tqdm
from nexus.util import *
import math
from scipy.stats.stats import pearsonr, spearmanr
from scipy.special import comb
from scipy.stats import hypergeom
import multiprocessing
import networkx
import obonet
import statsmodels.stats.multitest as multitest
import dcor
from minepy import MINE
from sys import getsizeof


def is_constant(vec):
    first = vec[0]
    for i in range(len(vec)):
        if vec[i] != first:
            return False
    return True

def calc_sobolev(reads1,reads2):
    z1 = ((np.square(reads1) / np.sum(np.square(reads1))) -
        (np.square(reads2) / np.sum(np.square(reads2))))
    fourier_transform = np.fft.fft(z1)
    w = (2*math.pi
        *(np.array([x for x in range(0,len(fourier_transform))], dtype=np.float32))
        /(len(fourier_transform)))
    s = np.square(np.sum((1+w)*np.square(np.abs(fourier_transform))),0.5)
    return s

def calc_fisher_information(reads1,reads2):
    t = ((np.square(reads1) / np.sum(np.square(reads1))) *
        (np.square(reads2) / np.sum(np.square(reads2))))
    return np.arccos(np.sum(np.sqrt(t)))

def calc_mic(reads1,reads2,mine):
    mine.compute_score(reads1, reads2)
    return mine.mic()

def prs(reads1,reads2):
    pearson_corr, p = pearsonr(reads1, reads2)
    return pearson_corr

def spr(reads1,reads2):
    spearman_corr, p = spearmanr(reads1, reads2)
    return spearman_corr

def dc(reads1,reads2):
    return dcor.distance_correlation(reads1, reads2)

def leave_one_out(pid, coding_rows, regulators, show, method_ids, return_dict):
    coding_noncoding_pairs = []
    mine = MINE()
    mic = lambda reads1,reads2: calc_mic(reads1,reads2,mine)

    method_names = {"MIC": mic, "PRS": prs, "SPR": spr, "DC": dc, 
                    "FSH": calc_fisher_information, "SOB": calc_sobolev}
    methods = []
    names = []
    for method_name in method_ids:
        methods.append(method_names[method_name])
    '''if show:
        print("Regulators list:\n\t"+str(regulators.head())+"\nLen="+str(len(regulators)))
        print(str(coding_rows.head())+"\nLen="+str(len(coding_rows)))'''
    for i in get_iterator(range(len(regulators)), show=show):
        row2 = regulators.iloc[i]
        reads2 = np.array(row2.values[1:],dtype=np.float32)
        rows = coding_rows[coding_rows[coding_rows.columns[0]] != row2[coding_rows.columns[0]]]
        for i in range(len(rows)):
            row1 = coding_rows.iloc[i]
            reads1 = np.array(row1.values[1:],dtype=np.float32)
            for i in range(len(method_ids)):
                corr = methods[i](reads1,reads2)
                if(corr > 0.5) or (corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], corr, method_ids[i]))
            '''if "MIC" in methods:
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
            if "FSH" in methods:
                fisher_information = calc_fisher_information(reads1,reads2)
                if(fisher_information > 0.5) or (fisher_information < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], fisher_information, "FSH"))
            if "SOB" in methods:
                sobolev_metric = calc_sobolev(reads1,reads2)
                if(sobolev_metric > 0.5) or (sobolev_metric < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], sobolev_metric, "SOB"))'''
            
    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def try_find_coexpression_process(pid, coding_rows, nc_rows, show, methods, return_dict):
    coding_noncoding_pairs = []
    mine = MINE()
    mic = lambda reads1,reads2: calc_mic(reads1,reads2,mine)

    method_names = {"MIC": mic, "PRS": prs, "SPR": spr, "DC": dc, 
                    "FSH": calc_fisher_information, "SOB": calc_sobolev}
    methods = []
    names = []
    for method_name in method_ids:
        methods.append(method_names[method_name])
    #print("Methods="+str(methods))
    for i in get_iterator(range(len(coding_rows)), show=show):
        row1 = coding_rows.iloc[i]
        reads1 = np.array(row1.values[1:],dtype=np.float32)
        for name, row2 in nc_rows.iterrows():
            reads2 = np.array(row2.values[1:],dtype=np.float32)
            for i in range(len(rows)):
            row1 = coding_rows.iloc[i]
            reads1 = np.array(row1.values[1:],dtype=np.float32)

            for i in range(len(method_ids)):
                corr = methods[i](reads1,reads2)
                if(corr > 0.5) or (corr < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], corr, method_ids[i]))
            '''if "MIC" in methods:
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
            if "FSH" in methods:
                fisher_information = calc_fisher_information(reads1,reads2)
                if(fisher_information > 0.5) or (fisher_information < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], fisher_information, "FSH"))
            if "SOB" in methods:
                sobolev_metric = calc_sobolev(reads1,reads2)
                if(sobolev_metric > 0.5) or (sobolev_metric < -0.5):
                    coding_noncoding_pairs.append((row1[0], row2[0], sobolev_metric, "SOB"))'''
    print(str(len(coding_noncoding_pairs)) + " correlation pairs found.")
    return_dict[pid] = coding_noncoding_pairs

def get_ancestors(graph, parent_id):
    ancestors = networkx.ancestors(graph, parent_id)
    return ancestors

def get_descendants(graph, parent_id):
    descendants = networkx.descendants(graph, parent_id)
    return descendants

def pvalue(m, N, M, n):
    '''
    M = number of successes in the population
    m = number of drawn successes
    N = sample size
    n = n
    '''
    '''comb_sum = 0
    for i in range(m, min(n,M)+1):
        comb_sum += comb(M, i, exact=True)*comb(N-M, n-i, exact=True)
    p = comb_sum / comb(N, n, exact=True)
    return p'''
    #return hypergeom.sf(m, N, M, n)
    return hypergeom.sf(m, N, M, n)

def pvalue_process(N, params, pid, return_dict):
    '''gene_term_pvalue = []
    n_lens = [len(genes_coexpressed_with_ncRNA[possible_gene_term[i][0]]) for i in range(len(possible_gene_term))]
    M_lens = [len(genes_annotated_with_term[possible_gene_term[i][1]]) for i in range(len(possible_gene_term))]
    m_lens = [len(genes_annotated_with_term[possible_gene_term[i][1]]
                .intersection(genes_coexpressed_with_ncRNA[possible_gene_term[i][0]]))
                for i in range(len(possible_gene_term))]'''
    '''if print_progress:
        print("Printing progress for thread " + str(pid))
    for i in get_iterator(possible_gene_term,show=False):
        #gene_term, n, M, m = possible_gene_term[i]
        #n = n_lens[i] # number of genes in test set
        n = possible_gene_term[i][1]
        #if n >= 5:
        #M = M_lens[i] # number of genes with annotation
        M = possible_gene_term[i][2]
        #    if M >= 5:
        #m = m_lens[i] # number of genes in test set with annotation
        m = possible_gene_term[i][3]
        #if n_lens[i] >= 5 and M_lens[i] >= 5 and m_lens[i] >= 1:
        pval = pvalue(m, N, M, n)
        gene_term_pvalue.append((possible_gene_term[i][0][0],
                            possible_gene_term[i][0][1],pval))'''
    return_dict[pid] = [(i,pvalue(m, N, M, n)) for i, n, M, m in params]

def get_valid_associations(genes_coexpressed_with_ncRNA, genes_annotated_with_term, possible_gene_term,
                            min_n, min_M, min_m):
    n_lens = [len(genes_coexpressed_with_ncRNA[gene]) for gene, term in possible_gene_term]
    M_lens = [len(genes_annotated_with_term[term]) for gene, term in possible_gene_term]
    m_lens = [len(genes_annotated_with_term[term]
                .intersection(genes_coexpressed_with_ncRNA[gene]))
                for gene, term in possible_gene_term]
    valid = [n_lens[i] >= min_n and M_lens[i] >= min_M and m_lens[i] >= min_m for i in range(len(possible_gene_term))]
    return valid, n_lens, M_lens, m_lens


def parallel_pvalues(N, possible_gene_term, valid_gene_term, 
                    n_lens, M_lens, m_lens, 
                    threads, available_memory):
    #basic_process_mem_usage = (getsizeof(N)*2 + getsizeof(genes_coexpressed_with_ncRNA)
    #    + getsizeof(genes_annotated_with_term))
    #available_memory -= (basic_process_mem_usage * threads)
    min_for_chunks = threads * 24 * 1024
    if available_memory < min_for_chunks:
        print("Setting available memory to minimum")
        available_memory = min_for_chunks
    print("Memory available for pvalue calculation: " + str(available_memory/1024))
    print("Dividing chunks for parallel p-value calculation")
    #gene_terms = []
    print("Possible associations: " + str(len(possible_gene_term)))
    params = []
    for i in range(len(possible_gene_term)):
        if valid_gene_term[i]:
            #gene_terms.append(possible_gene_term[i])
            params.append((i,n_lens[i],M_lens[i],m_lens[i]))
    print("Calculating pvalues for " + str(len(params)) + " associations.")
    one_size = getsizeof(params[0])
    print("Size of one element: " + str(one_size/1024))
    params_per_run = int(available_memory / one_size)
    run_chunks = list(chunks(params, params_per_run))

    index_pvalue = []
    print("Pvalues per chunk: " + str(len(run_chunks[0])))
    for run_chunk in tqdm(run_chunks):
        #
        #run_chunk = list(filter(lambda x: x == True))
        calcs_per_thread = int(len(run_chunk) / threads)+1
        calc_chunks = list(chunks(run_chunk, calcs_per_thread))
        #print("Calculating p-values")

        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        processes = []
        showing = False
        i = 0
        for chunk in calc_chunks:
            p = multiprocessing.Process(target=pvalue_process,
                    args=(N, chunk, i, return_dict, ))
            processes.append(p)
            #print("Spawned process")
            p.start()
            if showing:
                showing = False
            i += 1

        for p in processes:
            p.join()
            #print("Finished process " + str(p.pid))

        #print("Merging results")
        for value in return_dict.values():
            index_pvalue += value
        #print("Saved results from processes")
        manager.shutdown()
    print("Aggregating pvalue information to (gene,term) pairs")
    gene_term_pvalue = []
    for i in tqdm(range(len(index_pvalue))):
        index,pvalue = index_pvalue[i]
        if pvalue != np.nan:
            gene,term = possible_gene_term[index]
            gene_term_pvalue.append((gene,term,float(pvalue)))
    return gene_term_pvalue
