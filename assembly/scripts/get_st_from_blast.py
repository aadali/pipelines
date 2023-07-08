import sys
import json
from os import path

import pandas as pd

"""
the below lines is the headers of blast out and what does it mean
sseqid slen length nident pident qseqid qstart qend qseq sstrand

 qseqid: Query Seq-id
 sseqid: Subject Seq-id
   slen: Subject sequence length
 qstart: Start of alignment in query
   qend: End of alignment in query
 sstart: Start of alignment in subject
   send: Send of alignment in subject
 length: Alignment length
 pident: Percentage of identical matches
sstrand: Subject Strand
"""

blast_out = "/home/a/big/ycq/analysis_out/test001/barcode22/03.annotate/blast.out"


def print_log(msg, log="error"):
    # 31 red; 32 green; 33 yellow; 34 blue;
    if log == "warning":
        print(f"\033[33m{msg}\033[0m")
    elif log == "error":
        print(f"\033[31m{msg}\033[0m")
        sys.exit(1)


class AlignAllele(object):
    """
    Struct to store a locus align info
    """

    def __init__(self, locus, slen, allele, qseqid, qstart, qend, pident, strand):
        self.locus = locus
        self.slen = slen
        self.allele = allele
        self.qseqid = qseqid
        self.qstart = qstart
        self.qend = qend
        self.pident = pident
        self.strand = strand

    def __repr__(self):
        info = dict(
            locus=self.locus,
            slen=self.slen,
            allele=self.allele,
            qseqid=self.qseqid,
            qstart=self.qstart,
            qend=self.qend,
            pident=self.pident,
            strand=self.strand
        )
        return f"{info}"


class ParseMLST(object):
    """
    data structure of st_alleles_js is like this:
    ########################################################
       {
            "Achromobacter_spp.": {
                "scheme": "eno;gltB;lepA;nrdA;nuoL;nusA;rpoB",
                "eno@2|gltB@4|lepA@59|nrdA@5|nuoL@2|nusA@1|rpoB@26": "2",
                "eno@2|gltB@2|lepA@62|nrdA@1|nuoL@8|nusA@1|rpoB@26": "3",
                "eno@2|gltB@2|lepA@68|nrdA@1|nuoL@8|nusA@1|rpoB@26": "4",
                "eno@7|gltB@2|lepA@59|nrdA@2|nuoL@3|nusA@1|rpoB@26": "5",
            },
            "Brachyspira_hyodysenteriae": {
                "scheme": "Bhy_adh;Bhy_alp;Bhy_est;Bhy_gdh;Bhy_glpK;Bhy_pgm;Bhy_thi",
                "Bhy_adh@2|Bhy_alp@6|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@17|Bhy_pgm@4|Bhy_thi@3": "2",
                "Bhy_adh@2|Bhy_alp@6|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@17|Bhy_pgm@10|Bhy_thi@20": "3",
                "Bhy_adh@2|Bhy_alp@3|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@4|Bhy_pgm@3|Bhy_thi@16": "4",
                "Bhy_adh@2|Bhy_alp@10|Bhy_est@3|Bhy_gdh@10|Bhy_glpK@8|Bhy_pgm@1|Bhy_thi@3": "5",
                "Bhy_adh@2|Bhy_alp@7|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@12|Bhy_pgm@1|Bhy_thi@3": "6",
            }
        }
    }
    ########################################################
    construct the string key from blast_out like above to find which st this query is.
    each specie contain a key named "scheme" whose value separated by semi meaning what locus is used to determine the ST of this specie
    the part before @ in the string key is locus name;
    the part after @ in the string key is the allele of this locus; the value is ST
    all locus is sorted by the name of gene
    :param blast_out:
    :param st_alleles_js:
    :param mlst: write the mlst result into this file
    :return:
    """

    def __init__(self, blast_out: str, st_alleles_js: str, mlst: str):
        self.blast_df = pd.read_csv(blast_out, sep="\t", header=None,
                                    usecols=range(10),
                                    names=['qseqid', 'sseqid', 'slen', 'qstart', 'qend',
                                           'sstart', 'send', 'length', 'pident', 'strand'],
                                    dtype={
                                        "qseqid": str,
                                        "sseqid": str,
                                        "slen": int,
                                        "qstart": int,
                                        "qend": int,
                                        "sstart": int,
                                        "send": int,
                                        "length": int,
                                        "pident": float,
                                        "strand": str
                                    })

        self.st_alleles = json.load(open(st_alleles_js, 'r'))
        self.mlst = mlst

    def _get_specie(self):
        scov = round((abs(self.blast_df['sstart'] - self.blast_df['send']) + 1) / self.blast_df['slen'], 3)
        criteria1 = scov == 1  # the alignment must cover the full length of the allele
        criteria2 = self.blast_df['pident'] > 95
        sub_blast_df = self.blast_df[criteria1 & criteria2].reset_index(drop=True)
        return sub_blast_df

    def write_mlst(self):
        sub_blast_df = self._get_specie()
        if sub_blast_df.empty:
            print_log(
                "No subseq can cover the full length of the all locus alleles or all the pident is less than 95%. "
                "Drop all blast record", log="warning")
            with open(self.mlst, 'w', encoding='utf-8') as outfile:
                outfile.write("#NOT FOUND MLST FOR THIS ASSEMBLY")
            return

        # sseqid format: {specie}@{locus}_{allele}
        locus_allele = sub_blast_df['sseqid'].str.rsplit(pat="_", n=1, expand=True)
        locus_allele.columns = ['specie_locus', 'allele']
        specie_locus = locus_allele['specie_locus'].str.rsplit(pat='@', n=1, expand=True)
        specie_locus.columns = ['specie', 'locus']

        sub_blast_df.insert(loc=0, column="specie", value=specie_locus['specie'])
        sub_blast_df.insert(loc=1, column="locus", value=specie_locus['locus'])
        sub_blast_df.insert(loc=2, column="allele", value=locus_allele['allele'])
        grouped = sub_blast_df.groupby("specie")

        # extract_match_alleles: store the extract match scheme allele, the scov and pident of each record is 100%
        # not_extract_match_alleles: store the match scheme allele, one or some scov and pident of records is(are) less than 100
        # if the former is not empty, the later will be ignored
        # else the lager will be used to instead of the former
        # if the two list is both empty, this means NO MLST WAS FOUND
        extract_match_alleles = []
        not_extract_match_alleles = []
        for specie in grouped.groups:
            this_specie_blast = grouped.get_group(specie)
            this_specie_blast = this_specie_blast.sort_values(by="pident", ascending=False)
            specie_locus = grouped.get_group(specie)['locus'].unique().tolist()
            scheme = ";".join(sorted(specie_locus))  # scheme is separated by ";" and must be sorted

            # scheme of this specie NOT FOUND
            if scheme != self.st_alleles[specie]['scheme']:
                print_log(f"[INFO] {specie} dropped. Expect scheme: {self.st_alleles[specie]['scheme']}, "
                          f"But Got scheme: {scheme}", log="warning")
                continue

            # scheme of this specie FOUND.
            # There two situation:
            #   1) all locus allele extract match
            #   2) the pident of one or more locus allele is less than 100%
            else:
                align_alleles = []
                for locus in scheme.split(";"):
                    this_specie_locus = this_specie_blast[(this_specie_blast["locus"] == locus) &
                                                          (this_specie_blast["specie"] == specie)]
                    first_line = this_specie_locus.head(1).reset_index(drop=True)
                    align_allele = AlignAllele(locus=locus,
                                               allele=first_line['allele'][0],
                                               slen=first_line['slen'][0],
                                               qseqid=first_line['qseqid'][0],
                                               qstart=first_line['qstart'][0],
                                               qend=first_line['qend'][0],
                                               pident=first_line['pident'][0],
                                               strand=first_line['strand'][0])
                    align_alleles.append(align_allele)
                st_key = "|".join(
                    [f"{it.locus}@{it.allele}" for it in align_alleles]
                )  # construct the string key in st_alleles_js

                if st_key in self.st_alleles[specie]:  # situation 1: extract match
                    st = self.st_alleles[specie][st_key]
                    for align_allele in align_alleles:
                        extract_match_alleles.append("\t".join([
                            align_allele.qseqid,
                            str(align_allele.qstart),
                            str(align_allele.qend),
                            specie,
                            align_allele.locus,
                            str(align_allele.slen),
                            str(align_allele.allele),
                            align_allele.strand,
                            str(align_allele.pident),
                            "1",
                            f"{st}"
                        ]))
                else:  # situation 2: not all extract match
                    for align_allele in align_alleles:
                        not_extract_match_alleles.append("\t".join([
                            align_allele.qseqid,
                            str(align_allele.qstart),
                            str(align_allele.qend),
                            specie,
                            align_allele.locus,
                            str(align_allele.slen),
                            str(align_allele.allele) if align_allele.pident == 1 else f"~{align_allele.allele}",
                            align_allele.strand,
                            str(align_allele.pident),
                            "100",
                            "~"
                        ]))
        if extract_match_alleles:
            with open(self.mlst, 'w', encoding='utf-8') as outfile:
                outfile.write("\t".join(
                    ["QuerySeqID", "QueryStart", "QueryEnd", "Specie", "Locus", "SubjectLen",
                     "Allele", "Strand", "PercentageIdent", "SubjectCov", "St"]) + "\n")
                outfile.write("\n".join(extract_match_alleles) + "\n")
            return

        if not_extract_match_alleles:
            with open(self.mlst, 'w', encoding='utf-8') as outfile:
                outfile.write("\t".join(
                    ["QuerySeqID", "QueryStart", "QueryEnd", "Specie", "Locus", "SubjectLen",
                     "Allele", "Strand", "PercentageIdent", "SubjectCov", "St"]) + "\n")
                outfile.write("\n".join(not_extract_match_alleles) + "\n")
            return

        with open(self.mlst, 'w', encoding='utf-8') as outfile:
            outfile.write("#All species were dropped because of no scheme was found\n")
            return


if __name__ == '__main__':
    # parser = ParseMLST(blast_out, "/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/st_alleles.json",
    #                    "/home/a/Desktop/a.tsv")
    # parser.write_mlst()
    useage = f"python3 {path.basename(__file__)} <blast_out> <st_alleles.json> <mlst_result.out>"
    if len(sys.argv) != 4:
        print_log(useage, log="error")
    else:
        parser = ParseMLST(sys.argv[1], sys.argv[2], sys.argv[3])
        parser.write_mlst()
