import sys
import re
import pandas as pd

"""
从krakren2输出的分类文件里找到丰度最高的属的taxID以及以该taxID为根的所有种、亚种等的所有taxID，
用来从kraken2的reads分类里筛选对应的reads名称，继而从原始fastq里提取这些reads用于denovo
"""


def get_mos_g():
    reads_classified = sys.argv[1]
    classified = sys.argv[2]
    taxon_reads = sys.argv[3]
    classified_df = pd.read_csv(classified, sep="\t", header=None)
    classified_df.columns = ['abundance', 'reads_num_to_clade', 'reads_num_to_taxon', 'classified', 'taxon', 'name']
    classified_subdf = classified_df.query("classified == 'G'").head(2)
    taxon_index = classified_subdf.index.to_list()
    print(taxon_index)
    taxon_of_most_g = classified_df['taxon'].to_list()[taxon_index[0]: taxon_index[1]]
    print(taxon_of_most_g)
    reads_classified_df = pd.read_csv(reads_classified, sep="\t", header=None)
    reads_classified_df.columns = ['is_classified', 'read_name', 'taxon', 'length', 'lca']
    reads_classified_subdf = reads_classified_df.query("taxon in @taxon_of_most_g")
    reads_classified_subdf[['read_name', 'taxon']].to_csv(taxon_reads, sep="\t", index=False, header=None)


if __name__ == '__main__':
    get_mos_g()
