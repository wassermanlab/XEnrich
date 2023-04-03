#! /usr/bin/env python

import os
import pandas as pd
import numpy as np
import re

def __extract_results(results_dir, query_size, human_TFs):
    gresults = []
    enriched_TFs = []
    depleted_TFs = []

    for file in os.listdir(results_dir):
        df = pd.read_table(os.path.join(results_dir, file), index_col=False)
        df = __readable(df, human_TFs)
        gresults.append(df)
    
        bg_size = int(re.findall(r'\d+', file.split("_")[-1])[0])
        enriched_TFs.append(__get_enriched(df, query_size, bg_size)[['TF', 'overlap_ratio']])
        depleted_TFs.append(__get_depleted(df, query_size, bg_size)[['TF', 'overlap_ratio']])

    df_enrichedTFs = __dfsummary(enriched_TFs, human_TFs)
    df_enrichedTFs['obs/exp'] = df_enrichedTFs['overlap_ratio'] / (1/6)
    df_depletedTFs = __dfsummary(depleted_TFs, human_TFs)
    df_depletedTFs['obs/exp'] = df_depletedTFs['overlap_ratio'] / (1/6)

    return([gresults, df_enrichedTFs, df_depletedTFs])


#Convert results table file name to more readable
def __readable(df, human_TFs):
    bed = [ b[b.rindex('/')+1:] for b in df['#file'] ]
    df['TF'] = [ a[0:a.rindex('.bed')] for a in bed ]
    df['Type'] = np.where(df['TF'].isin(human_TFs), 'TF', 'other')
    
    c = list(df.columns)
    c = c[-2:] + c[:-2]
    c.remove('#file')
    df = df[c]
    df['combo_right'] = -np.log10(df['fishers_right_tail']) * np.log2(df['odds_ratio'])
    df['combo_left'] = -np.log10(df['fishers_left_tail']) * np.log2(df['odds_ratio'])
    df['overlap_ratio'] = df['overlaps'] / df['file_size']
    
    return(df)


def __get_enriched(df, querysize, bgsize):
    enriched = df[(df['odds_ratio'] > 1) & 
                  (df['combo_right'] >= df['combo_right'].quantile(0.85)) & 
                  (df['fishers_left_tail'] >= df['fishers_left_tail'].quantile(0.85)) &
                  (df['overlap_ratio'] >= querysize/bgsize)]
    return(enriched)


def __get_depleted(df, querysize, bgsize):
    depleted = df[(df['odds_ratio'] < 1) & 
                  (df['combo_left'] <= df['combo_left'].quantile(0.15)) & 
                  (df['fishers_right_tail'] >= df['fishers_right_tail'].quantile(0.85)) &
                  (df['overlap_ratio'] <= querysize/bgsize)]
    return(depleted)


def __dfsummary(TF_list, human_TFs):
    df = pd.concat(TF_list, ignore_index=True)
    df = pd.DataFrame(
        df.value_counts('TF'), columns=['counts']).join(
            pd.DataFrame(df.groupby('TF').mean('overlap_ratio'))).sort_values('TF').reset_index()
    df['TypeTF'] = df['TF'].isin(human_TFs)

    return(df)


def __extract_conEnDep(df_enrich, df_deplete, iter_cutoff):
    conEnrich = df_enrich.query('counts >= @iter_cutoff')
    conDeplete = df_deplete.query('counts >= @iter_cutoff')

    return(conEnrich, conDeplete)
