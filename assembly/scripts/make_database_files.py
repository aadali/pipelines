import os
import re
import json
import sys

from os import path
from collections import defaultdict

from Bio import SeqIO

"""
THIS SCRIPT WILL NOT BE USED IN ASSEMBLY PIPELIE

There some functions were used to make some pubmlst and vfdb database files. 
Such as scheme.txt, locus_vs_scheme.json, st_alleles.json of pubmlst
Such as vfdb.json of vfdb
"""


def get_vfdb_info(setA_nt, setA_pro):
    """
    get vfdb.json to store all the info of virulence genes or protein
    :param setA_nt: the core virulence genes from VFDB
    :param setA_pro: the core virulence protein from VFDB
    :return:
    """
    nuc_records = SeqIO.parse(open(setA_nt, "r"), "fasta")
    pro_records = SeqIO.parse(open(setA_pro, "r"), "fasta")
    seq_dict = {}
    for nuc_record in nuc_records:
        seq_dict[nuc_record.id] = [str(nuc_record.seq)]

    for pro_record in pro_records:
        seq_dict[pro_record.id].append(str(pro_record.seq))

    pat = re.compile("(VFG\d+\(gb.+?\)) \((.+)\) .+ \[(.+) \((VF\d+)\) - (.+) \((VFC\d+)\)\] \[(.+)\]", re.IGNORECASE)
    records = SeqIO.parse(open(setA_pro, 'r', encoding='utf-8'), "fasta")
    vfg = {}
    for record in records:
        mat = pat.search(record.description)
        try:
            gene_id = mat.group(1)
            gene_name = mat.group(2)
            vf_name = mat.group(3)
            vf_id = mat.group(4)
            vfc_name = mat.group(5)
            vfc_id = mat.group(6)
            bacteria = mat.group(7)
            vfg[gene_id] = {
                "gene_id": gene_id,
                "gene_name": gene_name,
                "vf_name": vf_name,
                "vf_id": vf_id,
                "vfc_name": vfc_name,
                "vfc_id": vfc_id,
                "bacteria": bacteria,
                "nt_seq": seq_dict[gene_id][0],
                'pro_seq': seq_dict[gene_id][1]
            }
        except (IndexError, AttributeError) as e:
            print(e, file=sys.stderr)
            print(f"{record.name} error when construct VFG dict", file=sys.stderr)
            continue
    json.dump(vfg, open("vfdb.json", 'w'))


def download(xmlfile):
    """
    parse the dbases.xml from https://pubmlst.org/data , generate the download.sh to download the profile.tsv and
    each alleles fasta file
    After the download.sh generated, do `bash download.sh`
    :param xmlfile: https://pubmlst.org/static/data/dbases.xml
    """
    # xmlfile = "../database/pubmlst/mlst_dbases.xml"
    f = open(xmlfile, "r")
    xml = f.read()
    mat = re.findall("<species>(.+?)</mlst>", xml, flags=re.DOTALL)
    # shell = open("download.sh", "w", encoding='utf-8')
    for record in mat:
        specie = re.findall("(.+)\n<mlst>", record)
        specie = specie[0].replace(" ", "_").replace("(", "_").replace(")", "_").replace("/", "#")
        profies_url = re.findall("<profiles>\n<url>(.+?)</url>", record)
        locus = re.findall("<locus>\n(.+?)\n<url>", record)
        locus_url = re.findall("<locus>\n.+?\n<url>(.+?)</url>", record)
        if len(locus) != len(locus_url):
            raise Exception("For species: count of locus is not eq count of locus url")
        print(record)
        print(locus)
        print(locus_url)
        # break
        with open(f"{specie}.download.sh", "w", encoding="utf-8") as shell:
            shell.write(f"mkdir {specie}\n")
            shell.write(f"wget -O ./{specie}/profile.tsv {profies_url[0]}\n")
            for each_locus, each_locus_url in zip(locus, locus_url):
                shell.write(f"wget -O ./{specie}/{each_locus}.fasta {each_locus_url}\n")
                shell.write("\n\n")
        assert len(locus) == len(locus_url)


def merge_alleles(db_dir, pubmlst_alleles):
    """
    merge all alleles of all species to one fasta that will be used to make blast database
    :param db_dir: directory that contain all species profile and alleles from download function
    :param pubmlst_alleles: the file name of merged fasta
    :return:
    """
    lines = []
    species = os.listdir(db_dir)
    for specie in species:
        files = os.listdir(path.join(db_dir, specie))
        fasta_files = filter(lambda x: x.endswith(".fasta"), files)
        for fasta_file in fasta_files:
            with open(path.join(path.join(db_dir, specie, fasta_file)), 'r') as infile:
                for line in infile:
                    if line.startswith(">"):
                        line = ">" + specie + "@" + line[1:]
                    lines.append(line)
    with open(pubmlst_alleles, 'w') as outfile:
        outfile.write("".join(lines))


def get_scheme(xmlfile):
    """
    get tsv file: column1 is the scheme name, column2 is the locuses belong to the scheme
    :param xmlfile: dbases.xml from https://pubmlst.org/data
    :return:
    """
    # xml_file = "/home/a/big/ycq/projects/assembly/database/pubmlst/mlst_dbases.xml"

    file = open(xmlfile, 'r', encoding='utf-8')
    xml = file.read()

    mat = re.findall("<species>(.+?)</mlst>", xml, flags=re.DOTALL)
    schemes = {}
    for record in mat:
        specie = re.findall("(.+)\n<mlst>", record)
        specie = specie[0].replace(" ", "_").replace("(", "_").replace(")", "_").replace("/", "#")
        profies_url = re.findall("<profiles>\n<url>(.+?)</url>", record)
        locus = re.findall("<locus>\n(.+?)\n<url>", record)
        locus_url = re.findall("<locus>\n.+?\n<url>(.+?)</url>", record)
        schemes[specie] = locus

    with open("schemes.txt", 'w', encoding='utf-8') as outfile:
        for key, value in schemes.items():
            outfile.write(key + "\t")
            outfile.write(",".join(value) + '\n')


def get_st_alleles(database_dir, scheme_file):
    """
    get the json file. All the info about alleles vs st of each specie will be collected.
    the  unique  combination of alleles is the key, just like this:
        {
            "Achromobacter_spp.": {
                "eno@2|gltB@4|lepA@59|nrdA@5|nuoL@2|nusA@1|rpoB@26": "2",
                "eno@2|gltB@2|lepA@62|nrdA@1|nuoL@8|nusA@1|rpoB@26": "3",
                "eno@2|gltB@2|lepA@68|nrdA@1|nuoL@8|nusA@1|rpoB@26": "4",
                "eno@7|gltB@2|lepA@59|nrdA@2|nuoL@3|nusA@1|rpoB@26": "5",
            },
            "Brachyspira_hyodysenteriae": {
            "Bhy_adh@2|Bhy_alp@6|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@17|Bhy_pgm@4|Bhy_thi@3": "2",
            "Bhy_adh@2|Bhy_alp@6|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@17|Bhy_pgm@10|Bhy_thi@20": "3",
            "Bhy_adh@2|Bhy_alp@3|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@4|Bhy_pgm@3|Bhy_thi@16": "4",
            "Bhy_adh@2|Bhy_alp@10|Bhy_est@3|Bhy_gdh@10|Bhy_glpK@8|Bhy_pgm@1|Bhy_thi@3": "5",
            "Bhy_adh@2|Bhy_alp@7|Bhy_est@3|Bhy_gdh@1|Bhy_glpK@12|Bhy_pgm@1|Bhy_thi@3": "6",
            }
        }
    :param database_dir: the path where download.sh worked. Only the scheme directory allowed existing in this directory
    :param scheme_file: the outfile of get_scheme(xmlfile)
    :return:
    """
    dirs = os.listdir(database_dir)
    species = list(filter(lambda x: path.isdir(path.join(database_dir, x)), dirs))
    species_number = len(species)
    # specie= [specie for specie in species]
    d = {}
    with open(scheme_file, 'r') as infile:
        for line in infile:
            new_line = line.strip()
            key, values = new_line.split("\t")
            d[key] = values.split(",")

    if len(d) != species_number:
        raise Exception("the number specie in scheme_file is not equal species in database_dir")
    d2 = {}
    """
    NOTE: THE profile.tsv of Mycoplasma_hominis MAY BE eST, BUT WHAT WE NEED IS ST, NOT eST, 
    YOU SHOULD DOWNLOAD THIS PROFILE FROM pubmlst.com MANUALLY
    """
    count = 1
    for specie in species:
        locus_number = len(d[specie])
        with open(path.join(f"{database_dir}/{specie}", "profile.tsv"), 'r') as infile:
            lines = infile.readlines()
            total_st = len(lines) - 1
            header = lines[0].strip("\n").split("\t")
            locus = header[1:1 + locus_number]
            this_specie_scheme = ";".join(sorted(locus))
            for line in lines[1:]:
                new_line = line.strip("\n")
                st_alleles = new_line.split("\t")
                st = st_alleles[0]
                alleles = st_alleles[1:1 + locus_number]
                locus_allele = [f"{loci}@{allele}" for loci, allele in zip(locus, alleles)]
                st_key = "|".join(list(sorted(locus_allele)))
                if specie in d2:
                    if st_key in d2[specie]:
                        print(f"duplicate st key for {specie}: {st_key}: {st}\t\t\t{count}")
                        count += 1
                        raise Exception("duplicate st key")
                    d2[specie][st_key] = st
                else:
                    d2[specie] = {}
                    d2[specie]['scheme'] = this_specie_scheme
                    d2[specie][st_key] = st
        if len(d2[specie]) != total_st+1:
            print(f"ERROR: for {specie}, expect {total_st}, got {len(d2[specie])}")
    json.dump(d2, open("../database/pubmlst/st_alleles.json", 'w'), indent=4)


def get_locus_vs_scheme(scheme_file):
    """
    get the json file: each locus belongs to which schemes
    :param scheme_file: the outfile of get_scheme(xmlfile)
    :return:
    """
    d = defaultdict(list)
    with open(scheme_file, 'r') as infile:
        for line in infile:
            scheme, locuses = line.strip().split("\t")
            locuses = locuses.split(",")
            for locus in locuses:
                if scheme not in d[locus]:
                    d[locus].append(scheme)

    json.dump(d, open("locus_vs_scheme.json", 'w'), indent=4)


if __name__ == '__main__':
    pass
    # get_scheme(xmlfile="/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/mlst_dbase.xml")
    get_st_alleles("/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/db", "/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/schemes.txt")
    # merge_alleles("/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/db", "pubmlst_alleles.fasta")
    # get_locus_vs_scheme("/home/a/big/ycq/projects/pipelines/assembly/database/pubmlst/schemes.txt")
    # get_st_alleles("/home/a/big/ycq/db/assembly_database/pubmlst/db", "/home/a/big/ycq/db/assembly_database/pubmlst/schemes.txt")
