from ast import Or
import os
import click
import gzip
import json

from urllib import request
import requests

import pandas as pd
import numpy as np

from pyliftover import LiftOver
from pybedtools import BedTool
from Bio import SeqIO, SeqUtils
import pybiomart
from pyjaspar import jaspardb

CONTEXT_SETTINGS = {
    "help_option_names": ["-h", "--help"],
}

@click.command(no_args_is_help=True, context_settings=CONTEXT_SETTINGS)
@click.option(
    "-g", "--genome",
    help="genome (e.g. hg38).",
    default="hg38",
    show_default=True
)
@click.option(
    "--biotype",
    help="biotype (e.g. GM12878).",
    required=True
)
@click.option(
    "--taxid",
    help="Taxonomic identifier (e.g. 9606).",
    type=int,
    default=9606,
    show_default=True
)

def cli(**params):

    # SET UP DIRECTORIES
    output_dir = os.getcwd() #this script should be called from the top

    db_dir = os.path.join(output_dir, "db")
    if not os.path.isdir(db_dir):
        os.makedirs(db_dir)

    genome_dir = os.path.join(db_dir, params['genome'])
    if not os.path.isdir(genome_dir):
        os.makedirs(genome_dir)

    biotype_dir = os.path.join(db_dir, params['biotype'])
    if not os.path.isdir(biotype_dir):
        os.makedirs(biotype_dir)

    # Get major files for the project
    fa_file, cpg_file, select_file, status_file, db_dir, remap_dir = __get_db(params['genome'], params['biotype'], params['taxid'], output_dir)
    TSS_bedcols = ['chr', 'TSS500up', 'TSS500down', 'Gene_Name']

    # METADATA FOR X GENES
    meta = pd.read_excel(status_file, 
                        sheet_name='List', 
                        usecols=['Gene Name', 'Balaton consensus calls',  'Transcript Start', 'Transcript Stop', 'Y homology'], 
                        skiprows=[461]).astype({'Transcript Start':'int', 'Transcript Stop':'int'})
    meta.columns=meta.columns.str.replace(" ", "_")
    meta['chr'] = 'chrX'

    ### Replace names {PPP2R3B-AS1 -> LINC00685} and {SPANXB2 -> SPANXB1}
    meta.loc[2,'Gene_Name'] = 'LINC00685'
    meta.loc[959,'Gene_Name'] = 'SPANXB1'

    meta['group'] = np.where(meta['Balaton_consensus_calls'].isin(['E', 'Mostly E']),
                            'escape',
                            np.where(meta['Balaton_consensus_calls'].isin(['S', 'Mostly S']), 
                            'subject',
                            np.where(meta['Balaton_consensus_calls'].isin(['VE', 'Mostly VE']), 'variable', None)))

    ## DETERMINE TSS REGION
    ### Using select file might remove a couple rows b/c they don't have the same exact number of rows from some of our genes (perhaps the U6)
    select_TSS = pd.read_table(select_file, 
                                usecols=[2,3,4,5,12],
                                names=['chr', 'strand', 'txStart', 'txEnd', 'Gene_Name']).query('chr == "chrX"')
    select_TSS['TSS'] = np.where(select_TSS['strand'] == '-', select_TSS['txEnd'] , select_TSS['txStart'])
    select_TSS = select_TSS.drop(columns=['txStart', 'txEnd'])

    meta = meta.merge(select_TSS, how='left')

    ## DEAL WITH MISSING VALUES
    meta['hg19_TSS'] = np.where(meta['strand'] == '-', meta['Transcript_Stop'] , meta['Transcript_Start'])
    lo = LiftOver('hg19', params['genome'])

    missing = meta.query('TSS.isnull()').reset_index()
    missing['TSS'] = [ lo.convert_coordinate('chrX', coord)[0][1] 
                            if lo.convert_coordinate('chrX', coord) != [] 
                            else None 
                            for coord in missing['hg19_TSS'] ]
    missing.loc[106, 'TSS'] = 45851017 #one specific change in coord identified manually (for hg38)
    missing = missing.set_index('index') #matching index back to meta table

    meta['TSS'] = meta['TSS'].fillna(missing['TSS']).astype(int)
    meta['TSS500up'] = meta['TSS'] - 500
    meta['TSS500down'] = meta['TSS'] + 500
    meta['TSS_sequence'], meta['TSS_GCcontent'] = __get_gccontent(meta[TSS_bedcols], fa_file)
    
    ## CPG ISLAND OVERLAP TSS REGION
    cpg = pd.read_table(cpg_file, 
                        usecols=[1,2,3,4,5],
                        names=['chr', 'start', 'end', 'name', 'length']).query('chr == "chrX"')
    cpg.drop(columns=['length']).to_csv(os.path.join(genome_dir, "chrX_CGI.bed"), header=False, index=False, sep="\t")
    # cpg = pd.read_table(cpg_file, names=['bin', 'chr', 'start', 'end', 'name', 'length', 'cpgNum', 'gcNum', 'perCpg', 'perGc', 'obsExp'])
    
    cpg_BT = BedTool.from_dataframe(cpg[['chr', 'start', 'end']])
    tss_BT = BedTool.from_dataframe(meta[TSS_bedcols])

    tss_cpg = pd.read_table(tss_BT.intersect(cpg_BT, wao=True).fn, #wao to give both A & B in the output - to be able to match back to meta
                            usecols=[0,1,2,3,7], 
                            names=TSS_bedcols + ['cpgoverlap_bp']).drop_duplicates().groupby(TSS_bedcols).agg(sum).reset_index()
    tss_cpg = tss_cpg.assign(CpG_island=tss_cpg['cpgoverlap_bp'] != 0) 
    meta = meta.merge(tss_cpg, how='left') #includes both bop overlap & T/F for CpG Island

    ### Dropping duplicates at this point, should reduce the number of genes down by about 10 -> 1133
    meta = meta.drop(columns=['Transcript_Stop', 'Transcript_Start', 'hg19_TSS']).drop_duplicates()

    ## SHELL REGIONS
    meta['shell10000up'] = meta['TSS500up'] - 10000
    meta['shell10000down'] = meta['TSS500down'] + 10000
    
    UP_bedcols = ['chr', 'shell10000up', 'TSS500up', 'Gene_Name']
    meta['UP_sequence'], meta['UP_GCcontent'] = __get_gccontent(meta[UP_bedcols], fa_file)

    DOWN_bedcols = ['chr', 'TSS500down', 'shell10000down', 'Gene_Name']
    meta['DOWN_sequence'], meta['DOWN_GCcontent'] = __get_gccontent(meta[DOWN_bedcols], fa_file)

    # TF METADATA
    remap_tfs_files = sorted(os.listdir(remap_dir))
    remap_tfs = [ s.replace('.bed.gz', '') for s in remap_tfs_files]

    ## RELEVENT INFO FROM ENSEMBL & JASPAR
    df_ensembl, df_panther = __get_tfclass(remap_tfs, params['taxid'])
    df_jaspar = __get_jaspar(params['taxid'])

    tf_meta = df_ensembl.drop(columns=['UniProtKB']).drop_duplicates().merge(
        df_ensembl.merge(df_panther, how='left').dropna(), how='left'
    ).sort_values('TF_Name').reset_index(drop=True)
    
    tf_meta = tf_meta.merge(df_jaspar[['TF_Name', 'JASPAR2022_Collection']].drop_duplicates(), how='left', on='TF_Name')

    meta.to_csv(os.path.join(db_dir, "meta_genes.tsv"), index=False, sep="\t")
    tf_meta.to_csv(os.path.join(db_dir, "meta_tfs.tsv"), index=False, sep="\t")


def __get_db(genome, biotype, taxid, output_dir):

    db_dir = os.path.join(output_dir, "db")
    if not os.path.isdir(db_dir):
        os.makedirs(db_dir)

    genome_dir = os.path.join(db_dir, genome)
    if not os.path.isdir(genome_dir):
        os.makedirs(genome_dir)

    biotype_dir = os.path.join(db_dir, biotype)
    if not os.path.isdir(biotype_dir):
        os.makedirs(biotype_dir)

    status_file = __download_xcistatus(db_dir) #Supplementary from Balaton
#    screen_file = __download_screen(biotype, biotype_dir)
    fa_file, cpg_file, select_file = __download_ucsc(genome, genome_dir)

    remap_dir = __extract_remap(biotype, taxid, biotype_dir)
#    __extract_segway(biotype, biotype_dir, genome)
    
    return(fa_file, cpg_file, select_file, status_file, db_dir, remap_dir) #took out the screen_file for now (may not be needed)

def __download_xcistatus(db_dir):
    status_url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4696107/bin/13293_2015_53_MOESM1_ESM.xlsx"
    status_file = os.path.join(db_dir, "Concensus_XCI_Status.xlsx")
    download_from_url(status_url, status_file)

    return(status_file)

def __download_screen(biotype, biotype_dir):
    if biotype == 'GM12878':
        screen_url = "https://api.wenglab.org/screen_v13/fdownloads/Seven-Group/ENCFF915DFR_ENCFF583MDQ_ENCFF340JIF_ENCFF852CRG.7group.bed"
        screen_file = os.path.join(biotype_dir, biotype + "_screen.bed.gz")
        download_from_url_togzip(screen_url, screen_file)

    return screen_file

def __download_ucsc(genome, genome_dir):
    goldenpath = "http://hgdownload.soe.ucsc.edu/goldenPath/" + genome + "/"
    genome_prefix = genome + "_"

    ### GENOME FASTA FILE
    fa_url = os.path.join(goldenpath, "bigZips/" + genome + ".fa.gz")
    fa_file = os.path.join(genome_dir, fa_url.split('/')[-1]).replace(".gz", "")
    if not os.path.exists(fa_file):
        with request.urlopen(fa_url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                with open(fa_file, 'wb') as outfile:
                    outfile.write(uncompressed.read())

    ### SELECT FILES WITH TRANSCRIPTS INFO 
    select_url = os.path.join(goldenpath, "database/ncbiRefSeqSelect.txt.gz")
    select_file = os.path.join(genome_dir, genome_prefix + select_url.split('/')[-1])
    download_from_url(select_url, select_file)

    ### CPG ISLANDS
    cpg_url = os.path.join(goldenpath, "database/cpgIslandExt.txt.gz")
    cpg_file = os.path.join(genome_dir, genome_prefix + cpg_url.split('/')[-1])
    download_from_url(cpg_url, cpg_file)

    return(fa_file, cpg_file, select_file)

def __extract_remap(biotype, taxid, bio_dir):
    # nonredun_dir = os.path.join(bio_dir, 'ReMap2020_NonRedun')
    nonredun_dir = os.path.join(bio_dir, 'ReMap2022_NonRedun')
    if not os.path.isdir(nonredun_dir):
        os.makedirs(nonredun_dir)

    chromosomes = [ "chr" + str(num) for num in list(range(1, 23)) ] + ["chrX"]

    if len(os.listdir(nonredun_dir)) == 0:

        ### DOWNLOAD REMAP FILE
        remap_api = "https://remap.univ-amu.fr/api/v1/datasets/findByBiotype/biotype=" + biotype + "&taxid=" + str(taxid)

        # nr_file = os.path.join(bio_dir, 'remap2020_GM12878_nr_macs2_hg38_v1_0.bed.gz') #REMAP2020
        nr_url = json.load(request.urlopen(remap_api)).get(biotype).get('bed_url').replace('_all_', '_nr_').replace('http', 'https') #REMAP2022
        nr_file = os.path.join(bio_dir, nr_url.split('/')[-1]) #REMAP2022
        download_from_url(nr_url, nr_file)
    
        remap_nrraw = pd.read_table(nr_file, usecols=[0,1,2,3], 
                                    names=['chr', 'start', 'end', 'TF'])
        remap_nrraw['TF'] = [ s.split(":")[0] for s in remap_nrraw['TF'] ]
        remap_nrraw = remap_nrraw.query('chr == @chromosomes')
        #__groupTFs(remap_nrraw, nonredun_dir)

        chrx_dir = os.path.join(bio_dir, nonredun_dir + "_chrX") # Make another directory for just chrX?
        if not os.path.isdir(chrx_dir):
            os.makedirs(chrx_dir)

        ### SPLIT FILE INTO TFS
        grouped_TF = remap_nrraw.groupby('TF')
        for group in grouped_TF:
            group[1].to_csv(os.path.join(nonredun_dir, group[0] + ".bed.gz"),
                            sep="\t", index=False, header=False, compression='gzip')
            group[1].query('chr == "chrX"').to_csv(os.path.join(chrx_dir, group[0] + ".bed.gz"),
                                                    sep="\t", index=False, header=False, compression='gzip')

    return(nonredun_dir)

def __get_tfclass(remap_tfs, taxid):
    chrom_symbol = [ str(s) for s in range(1, 23)]  + ["X"]

    en_server = pybiomart.Server(host='http://www.ensembl.org')
    en_mart = en_server['ENSEMBL_MART_ENSEMBL']['hsapiens_gene_ensembl']
    ensembl = en_mart.query(attributes=['external_gene_name', 
                                        'chromosome_name', 
                                        'entrezgene_id', 
                                        'ensembl_gene_id', 
                                        'hgnc_id',  
                                        'uniprot_gn_id', 
                                        'description'], use_attr_names=True)

    ### ZNF687 seem to be missing from BioMart (from query - can be found individually but not within the query)
    df_ensembl = ensembl.query(
        "(external_gene_name == @remap_tfs) & (chromosome_name == @chrom_symbol)"
    ).append(ensembl.query('external_gene_name == "ZNF687"')).reset_index(drop=True)
    df_ensembl.columns = ['TF_Name', 'chr', 'Entrez_ID', 'Ensembl_ID', 'HGNC_ID', 'UniProtKB', 'description']
    df_ensembl['HGNC_ID'] = [ s.replace('HGNC:', '') for s in df_ensembl['HGNC_ID'] ]
    df_ensembl['Entrez_ID'] = df_ensembl['Entrez_ID'].astype(int)

    ### Repetitive entries for MEF2B & BACH1 that needs to be taken care of
    df_ensembl = df_ensembl.query('Entrez_ID != [4207, 100379661]')

    ## PANTHER PROTEIN CLASS
    panther_response = requests.get("http://pantherdb.org/services/oai/pantherdb/geneinfo?", params={'geneInputList':", ".join(remap_tfs), 'organism':taxid})

    df_panther = pd.DataFrame(columns=['UniProtKB', 'protein_class'])
    panther_acc = []
    i=0
    for r in panther_response.json()['search']['mapped_genes']['gene']:  
        acc = r['accession'].split("=")[-1]
        panther_acc.append(acc)
        try:
            r['annotation_type_list']['annotation_data_type']
        except:
            continue
        else:
            for d in r['annotation_type_list']['annotation_data_type']:
                if d['content'] == "ANNOT_TYPE_ID_PANTHER_PC":
                    df_panther.loc[i, 'UniProtKB'] = acc
                    df_panther.loc[i, 'protein_class'] = d['annotation_list']['annotation']['name']
                    i += 1

    return(df_ensembl, df_panther)

def __get_jaspar(taxid):
    jdb2022 = jaspardb(release='JASPAR2022')

    core_motifs = jdb2022.fetch_motifs(tax_group = ['vertebrates'], collection="CORE", species=taxid)
    unvalid_motifs = jdb2022.fetch_motifs(tax_group = ['vertebrates'], collection="unvalidated", species=taxid)

    df2022 = pd.DataFrame(columns=['TF_Name', 'MatrixID', 'JASPAR2022_Collection', 'Class', 'Family', 'Accession'])

    df2022 = extract_motifs(core_motifs, df2022)
    df2022 = extract_motifs(unvalid_motifs, df2022)  
    df2022_uniq = df2022.query('TF_Name.str.contains("::") == False')

    return(df2022_uniq)

def extract_motifs(motif_list, motif_df):
    for motif in motif_list:
        df = {'TF_Name': motif.name, 'MatrixID': motif.matrix_id, 'JASPAR2022_Collection': motif.collection, 'Class': motif.tf_class, 'Family': motif.tf_family, 'Accession': motif.acc}
        motif_df = motif_df.append(df, ignore_index=True)

    return (motif_df)

def __makedirs(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

def download_from_url(url, filename):
    if not os.path.exists(filename):
        with request.urlopen(url) as response:
            with open(filename, 'wb') as outfile:
                outfile.write(response.read())

def download_from_url_togzip(url, filename):
    if not os.path.exists(filename):
        with request.urlopen(url) as response:
            with gzip.open(filename, 'wb') as outfile:
                outfile.write(response.read())

def __get_gccontent(df, fa_file): #df should have at least 3 columns (chr, start, end)
    df_BT = BedTool.from_dataframe(df)
    df_seq = df_BT.sequence(fi=fa_file)
    
    sequence = []
    gc = []
    fasta_sequences = SeqIO.parse(open(df_seq.seqfn),'fasta')
    for fasta in fasta_sequences:
        seq = str(fasta.seq).upper()
        sequence.append(seq)
        gc.append(round(SeqUtils.GC(seq)))
    return(sequence, gc)


if __name__ == "__main__":
    cli()


## NEED TO ADD VCF FOR GM12878
#ftp://ftp.sra.ebi.ac.uk/vol1/ERZ094/ERZ094050/NA12878.EVA.vcf.gz

#OR - Complete Genomics 2018
#ftp://ftp.sra.ebi.ac.uk/vol1/ERZ022/ERZ022005/vcfBeta-NA12878-200-37-ASM.vcf.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/ERZ021/ERZ021996/vcfBeta-NA12891-200-37-ASM.vcf.gz

