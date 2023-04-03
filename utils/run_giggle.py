import os
import pandas as pd
import numpy as np

from Bio import bgzf
from pybedtools import BedTool
from pyGiggle.giggle import Giggle

from utils.dbdownload import __get_gccontent

def main():

    #set up variables
    outdir = os.getcwd() #"./"
    __makedir(outdir)

    query_dir = os.path.join(outdir, "query")
    __makedir(query_dir)

    genome = "hg38"
    biotype = "GM12878"

    TSS_bedcols = ['chr', 'TSS500up', 'TSS500down', 'Gene_Name']
    TSS_cols = TSS_bedcols + ['TSS_sequence', 'TSS_GCcontent']

    #read in data
    meta = pd.read_table('./db/meta_genes.tsv')
    select_TSS = pd.read_table('./db/hg38/hg38_ncbiRefSeqSelect.txt.gz', 
                                usecols=[2,3,4,5,12],
                                names=['chr', 'strand', 'txStart', 'txEnd', 'Gene_Name']).query('chr == "chrX"')
    select_TSS['TSS'] = np.where(select_TSS['strand'] == '-', select_TSS['txEnd'] , select_TSS['txStart'])

    select_TSS['TSS500up'] = select_TSS['TSS'] - 500
    select_TSS['TSS500down'] = select_TSS['TSS'] + 500
    select_TSS['TSS_sequence'], select_TSS['TSS_GCcontent'] = __get_gccontent(select_TSS[TSS_bedcols], './db/hg38/hg38.fa')
    
    #separate meta into different categories
    meta_bed = out_query_bgzf(meta, "allX", genome, query_dir, TSS_bedcols)

    par_genes = meta.query('Balaton_consensus_calls == "PAR"')['Gene_Name'].tolist()

    chrX_noPAR = meta.query('Balaton_consensus_calls != "PAR"')
    chrX_noPAR_bed = out_query_bgzf(chrX_noPAR, "noPAR", genome, query_dir, TSS_bedcols)

    escape = meta.query('Balaton_consensus_calls == ["E", "Mostly E"]')
    escape_genes = escape['Gene_Name'].tolist()
    escape_bed = out_query_bgzf(escape, "escape", genome, query_dir, TSS_bedcols)
    
    escape_upshell_bed = out_query_bgzf(escape, "Upshell", genome, query_dir, ['chr', 'shell10000up', 'TSS500up', 'Gene_Name'])
    escape_downshell_bed = out_query_bgzf(escape, "Downshell", genome, query_dir, ['chr', 'TSS500down', 'shell10000down', 'Gene_Name']) 

    subject = meta.query('Balaton_consensus_calls == ["S", "Mostly S"]')
    subject_genes = subject['Gene_Name'].tolist()
    subject_bed = out_query_bgzf(subject, "subject", genome, query_dir, TSS_bedcols)

    #Running Giggle for the different bedfiles & make results directory
    # escape vs chrX
    gcmatched_escapeVSchrX = __bg_gcmatched(escape, chrX_noPAR.query('Gene_Name != @escape_genes'), 20, 5, 3, TSS_cols)
    escapeVSchrX_resdir = __giggle_results(escape_bed, gcmatched_escapeVSchrX, "escape", "chrX", "TSS", biotype, 1000, TSS_cols, outdir)

    #subject vs chrX
    gcmatched_subjectVSchrX = __bg_gcmatched(subject, chrX_noPAR.query('Gene_Name != @subject_genes'), 20, 5, 3, TSS_cols)
    subjectVSchrX_resdir = __giggle_results(subject_bed, gcmatched_subjectVSchrX, "subject", "chrX", "TSS", biotype, 1000, TSS_cols, outdir)

    #upshell escape vs chrX
    upshell_cols = ['chr', 'shell10000up', 'TSS500up', 'Gene_Name', 'UP_sequence', 'UP_GCcontent']
    gcmatched_upshell_escapeVSchrX = __bg_gcmatched(escape, chrX_noPAR.query('Gene_Name != @escape_genes'), 20, 5, 3, upshell_cols)
    upshell_escapeVSchrX_resdir = __giggle_results(escape_upshell_bed, gcmatched_upshell_escapeVSchrX, "escape", "chrX", "UPSHELL", biotype, 10000, upshell_cols, outdir)

    #downshell escape vs chrX
    downshell_cols = ['chr', 'TSS500down', 'shell10000down', 'Gene_Name', 'DOWN_sequence', 'DOWN_GCcontent']
    gcmatched_downshell_escapeVSchrX = __bg_gcmatched(escape, chrX_noPAR.query('Gene_Name != @escape_genes'), 20, 5, 3, downshell_cols)
    downshell_escapeVSchrX_resdir = __giggle_results(escape_downshell_bed, gcmatched_downshell_escapeVSchrX, "escape", "chrX", "DOWNSHELL", biotype, 10000, downshell_cols, outdir)


def __makedir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)


def out_query_bgzf(df, prefix, genome, outdir, bedcols):
    bed_ext = ".bed.gz"

    filename = "_".join([genome, prefix, str(df.shape[0]) + "TSS" + bed_ext])
    filepath = os.path.join(outdir, filename)

    if not os.path.exists(filepath):
        bgzf_file = bgzf.BgzfWriter(filepath)
        bgzf_file.write(df[bedcols].to_csv(sep="\t", header=False, index=False))
        bgzf_file.close()

    return filepath


def __bg_gcmatched(query, database, iterations, numgenes, precision, colnames):
    ### For each gene in the query, find a number of genes with the same GC content 
    ### within a certain precision from the background (database)
    ### Export a list of dfs for the number of iterations specified

    gc_col = [c for c in colnames if 'GCcontent' in c][0]
    list_df = []
    
    for i in range(1, iterations + 1):
        db = database
        out_df = pd.DataFrame(columns=colnames)
        for row in query[gc_col].value_counts().reset_index().iterrows():
            gc_range = list(range(row[1][0]-precision, row[1][0]+precision+1))
            count = int(row[1][1]) * numgenes
            same_gc = db[db[gc_col].isin(gc_range)][colnames]
            if count < same_gc.shape[0]:
                sample = same_gc.sample(count, random_state=i)
                db = db.drop(sample.index) #non replaced (same gene cannot be added more than once)
                out_df = out_df.append(sample[colnames])
            else:
                db = db.drop(same_gc.index)
                out_df = out_df.append(same_gc[colnames])
        out_df = out_df.append(query[colnames]).sort_values(colnames[0:2]).reset_index(drop=True)
        list_df.append(out_df)
    
    return list_df


def __giggle_results(queryfile, bglist, querytype, bgtype, regiontype, biotype, regionlength, colnames, outdir):
    ### Export GIGGLE results to a directory in an interpretable manner

    prefix =  regiontype + "_" + querytype + "VSgcmatched" + bgtype + "_NRuniq"

    remap_giggle_dir = os.path.join(outdir, 'ReMap2022_Giggle')
    if not os.path.isdir(remap_giggle_dir):
        os.makedirs(remap_giggle_dir)

    prefix_dir = os.path.join(remap_giggle_dir, prefix)
    if not os.path.isdir(prefix_dir):
        os.makedirs(prefix_dir)

    results_dir = os.path.join(prefix_dir, "results") 
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    elif os.listdir(results_dir) != 0:
        return(results_dir)

    db_dir = os.path.join(outdir, "db")
    biotype_dir = os.path.join(db_dir, biotype)
    remap_dir = os.path.join(biotype_dir, "ReMap2022_NonRedun") 
        
    query_size = pd.read_table(queryfile, names=colnames[:4]).shape[0]
    
    #iterating for each of the background for the GIGGLE operation
    for i in range(len(bglist)):
        itername = prefix + str(i + 1)
        iter_dir = os.path.join(prefix_dir, itername)
        if not os.path.isdir(iter_dir):
            os.makedirs(iter_dir)
    
        #exporting background to file
        df = bglist[i][colnames]
        df_file = os.path.join(prefix_dir, itername + "_bg.tsv")
        df.to_csv(df_file, sep="\t", header=False, index=False)
    
        bg_size = df.shape[0]
        bg_length = bg_size * regionlength
    
        df_BT = BedTool.from_dataframe(df) #background file
    
        for tfbs_file in os.listdir(remap_dir):

            TF = tfbs_file.split('.')[0]
            file_path = os.path.join(remap_dir, tfbs_file)

            if os.path.isdir(file_path):
                continue

            #For shell regions, we want to make sure we are not counting hits that were already counted in the TSS regions
            #Remove parts of TF binding regions that already overlap with the TSS regions
            if regiontype.find("SHELL") > -1:
                TF_BT = BedTool(file_path)
                tss_file = queryfile.replace(queryfile.split("/")[-1], "hg38_noPAR_1122TSS.bed.gz")
                TF_noTSS = TF_BT.intersect(tss_file, v=True) #get anything TF region that's doesn't overlap TSS region
                TF_noTSS_BT = BedTool(TF_noTSS)
                inter = df_BT.intersect(TF_noTSS_BT, wa=True)
            else:
                TF_BT = BedTool(file_path)
                inter = df_BT.intersect(TF_BT, wa=True, u=True)
            
            inter_df = pd.read_table(inter.fn, usecols=[0,1,2], 
                                 names=['chr', 'start', 'end']).sort_values(['chr', 'start'])
        
            if inter_df.empty:
                continue
                
            inter_bgzf = bgzf.BgzfWriter(os.path.join(iter_dir, tfbs_file))
            inter_bgzf.write(inter_df.to_csv(sep="\t", header=False, index=False))
            inter_bgzf.close()
    
        index_dir = os.path.join(prefix_dir, itername + "_index")
        Giggle.create(index_dir, os.path.join(iter_dir, "*.bed.gz"))
    
        result_file = '_'.join(["ReMap", str(query_size), itername, biotype, str(bg_size) + regiontype + ".result"])
        result_cmd = "~/Desktop/Wasserman/giggle/bin/giggle search -i " + index_dir + " -q " + queryfile + " -s -g " + str(bg_length) + " > " + results_dir + "/" + result_file  
        os.system(result_cmd)
        
    return(results_dir)

if __name__ == "__main__":
    main()