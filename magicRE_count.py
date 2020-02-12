"""
Perform an alignment of reads to the pseudogenome,
Count unique, fractional, and total reads mapped.
Write output to tsv files.
"""
import argparse
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path


def align_to_pseudogenome(out_dir, pseudogenome_db, sample_name, read_type, read_handle, read_handle2, threads, paired_end, max_words):
    """Align reads to the pseudogenome using magicBLAST"""
    outfile = Path(out_dir, f"__{sample_name}_magicBLAST_results.tab")

    if read_type == "sra":
        align_cmd = f"magicblast -sra {read_handle} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -splice F -max_db_word_count {max_words} > {outfile}"
    elif read_type == "fasta":
        if paired_end:
            align_cmd = f"magicblast -query {read_handle} -query_mate {read_handle2} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -no_discordant -splice F -max_db_word_count {max_words} > {outfile}"
        else:
            align_cmd = f"magicblast -query {read_handle} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -splice F -max_db_word_count {max_words} > {outfile}"
    elif read_type == "fastq":
        if paired_end:
            align_cmd = f"magicblast -query {read_handle} -query_mate {read_handle2} -db {pseudogenome_db} -infmt fastq -num_threads {threads} -outfmt tabular -no_unaligned -no_discordant -splice F -max_db_word_count {max_words} > {outfile}"
        else:
            align_cmd = f"magicblast -query {read_handle} -db {pseudogenome_db} -infmt fastq -num_threads {threads} -outfmt tabular -no_unaligned -splice F -max_db_word_count {max_words} > {outfile}"
    else:
        raise TypeError("Unknown read_type: must be one of [sra, fasta, fastq]")
    
    # run the alignment
    print("Aligning sample to the pseudogenome. Handing off to magicblast...")
    print(f"Effective magicBLAST command:\n{align_cmd}")
    subprocess.run(align_cmd, shell=True)

    return outfile


def combine_counts(counts1, counts2):
    """Generic function for combining two dictionaries of count data"""
    combined_counts = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

    for class_ in counts1.keys():
        for family in counts1[class_].keys():
            for elem in counts1[class_][family].keys():
                combined_counts[class_][family][elem] += counts1[class_][family][elem]
        
    for class_ in counts2.keys():
        for family in counts2[class_].keys():
            for elem in counts2[class_][family].keys():
                combined_counts[class_][family][elem] += counts2[class_][family][elem]
            
    return combined_counts


def count_magicBLAST_result(magicBLAST_outfile, id_threshold, debug):
    """Count the unique, total and fractional counts from the magicBLAST results"""
    print("Generating counts from alignment to the pseudogenome...")
    uniq_counts = defaultdict(float)
    total_counts = defaultdict(float)
    frac_counts = defaultdict(float)

    with open(magicBLAST_outfile, 'r') as infile:
        for line in infile:
            if line.startswith("#"):  # skip header lines
                continue
            else:
                l = line.strip().split()
                read_id =l[0]
                sub_family = l[1]
                percent_identity = float(l[2])
                N_alignments = int(l[17])

                # calculate total counts
                if percent_identity >= id_threshold:
                    total_counts[sub_family] += 1
                
                # calculate unique counts
                if percent_identity >= id_threshold and N_alignments == 1:
                    uniq_counts[sub_family] += 1

                # calculate fractional counts -- still must be combined with uniq counts
                counted = set()
                if percent_identity >= id_threshold and N_alignments > 1:
                    if read_id not in counted:
                        frac_counts[sub_family] += (1 / N_alignments)
                        counted.add(read_id)
    
    if not debug:
        print("Removing large intermediate files...")
        subprocess.run(f"rm {magicBLAST_outfile}", shell=True)

    return uniq_counts, frac_counts, total_counts


def create_element_mapping(repnames_bedfile):
    """Create a mapping of the element names to their classes and families"""
    elem_key = defaultdict(lambda : defaultdict(str))
    with open(repnames_bedfile, "r") as bed:
        for line in bed:
            l = line.strip().split("\t")
            name = l[3]
            class_ = l[4]
            family = l[5]
            elem_key[name]["class"] = class_
            elem_key[name]["family"] = family
    
    return elem_key


def map_counts(counts, element_mapper):
    """Create a dictionary of reads using the class > family > element hierarchy.
    Since the output from the alignment step function returns element names and their counts, 
    these element names must be re-associated with the classes and families they are part of. 
    This association is performed by the element_mapper.
    """
    hierarchy = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for elem in counts.keys():
        class_ = element_mapper[elem]['class']
        family = element_mapper[elem]['family']
        hierarchy[class_][family][elem] += counts[elem]

    return hierarchy


def summarize(count_data):
    """Collapse the counts to the class and family level."""
    class_data = defaultdict(float)
    family_data = defaultdict(lambda : defaultdict(float))
    
    # class level counts
    for class_ in count_data:
        for family in count_data[class_]:
            for elem in count_data[class_][family]:
                count = float(count_data[class_][family][elem])
                class_data[class_] += count
    
    # family level counts
    for class_ in count_data:
        for family in count_data[class_]:
            for elem in count_data[class_][family]:
                count = float(count_data[class_][family][elem])
                family_data[class_][family] += count
    
    return class_data, family_data


def write_class_summary(class_counts, outfile):
    """Write class-level counts to outfile. Output truncates floats to ints"""
    with open(outfile, "w") as f:
        for class_ in class_counts.keys():
            count = int(class_counts[class_])
            f.write(f"{class_}\t{count}\n")


def write_family_summary(family_counts, outfile):
    """Write family-level counts to outfile. Output truncates floats to ints."""
    with open(outfile, "w") as f:
        for class_ in family_counts.keys():
            for family in family_counts[class_]:
                count = int(family_counts[class_][family])
                f.write(f"{class_}\t{family}\t{count}\n")


def write_output(count_data, outfile):
    """Write count data to the outfile. Output truncates floats to ints."""
    print(f"Writing results to {outfile}...")
    with open(outfile, "w") as f:
        for class_ in count_data:
            for family in count_data[class_]:
                for elem in count_data[class_][family]:
                    count = int(count_data[class_][family][elem])
                    f.write(f"{class_}\t{family}\t{elem}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(description="Align reads to pseudogenome using magicBLAST and return counts of repeat elements")
    parser.add_argument('sampleName', help='The name of the sample to be processed.')
    parser.add_argument('setupDir', help="Path to the magicRE_Setup directory.")
    parser.add_argument('queryType', default='fastq', choices=['sra', 'fasta', 'fastq'], help='Defines the input format to be aligned. One of [SRA, fasta, fastq]. Defines the input format to be aligned.')
    parser.add_argument('query', help="Either the SRA accession number or the path to the fasta or fastq file to be aligned to the pseudogenome.")
    parser.add_argument('--query2', default=None, help="If pairedEnd is set, designate this to be the path to second read pair.")
    parser.add_argument('--percentIdentity', default=90, type=float, metavar='90.0', help='Threshold at which to include mapped reads in the count data. Mapped reads must be >= percentIdentity to be included.')
    parser.add_argument('--threads', default=1, type=int, help='Number of threads to use for alignment.')
    parser.add_argument('--pairedEnd', dest='pairedEnd', action='store_true', help='Designate this option for paired-end sequencing.')
    parser.add_argument('--summarize', dest='summarize', action='store_true', help='In addition to the repeat name level output, produce collapsed counts at the class and family levels for each of the count types.')
    parser.add_argument('--outDirectory', default=".", help="Specify where to save the count results. Defaults to the current directory/sampleName.")
    parser.add_argument('--debug', dest='debug', action='store_true', help='Select this option to prevent the removal of temporary files; useful for debugging.')
    parser.add_argument('--maxWords', default=9999999, type=int, metavar='9999999', help='16-base words that appear in the genome more than this number of times will be filtered. This is a work-around for setting -limit_lookup F (allowing repeats) since setting -limit_lookup F in the alignment command results in the process being killed.')
    parser.set_defaults(pairedEnd=False, summarize=False, debug=False)
    args = parser.parse_args()

    sample_name = args.sampleName
    setup_dir = args.setupDir
    query_type = args.queryType
    query = args.query
    query2 = args.query2
    perc_id = args.percentIdentity
    threads = args.threads
    paired_end = args.pairedEnd
    summarize_data = args.summarize
    out_path = args.outDirectory
    max_words = args.maxWords
    debug = args.debug

    # setup file i/o ----------------------------------------------------------
    out_dir = Path(out_path, sample_name + "_magicRE_results")
    pseudogenome_db = Path(setup_dir, "pseudogenome")
    repnames_bedfile = Path(setup_dir, "repnames.bed")
    total_counts_outfile = Path(out_dir, sample_name + "_total_counts.tsv")
    uniq_counts_outfile = Path(out_dir, sample_name + "_unique_counts.tsv")
    fractional_counts_outfile = Path(out_dir, sample_name + "_fractional_counts.tsv")
    class_total_counts_outfile = Path(out_dir, sample_name + "_class_total_counts.tsv")
    class_uniq_counts_outfile = Path(out_dir, sample_name + "_class_unique_counts.tsv")
    class_fractional_counts_outfile = Path(out_dir, sample_name + "_class_fractional_counts.tsv")
    family_total_counts_outfile = Path(out_dir, sample_name + "_family_total_counts.tsv")
    family_uniq_counts_outfile = Path(out_dir, sample_name + "_family_unique_counts.tsv")
    family_fractional_counts_outfile = Path(out_dir, sample_name + "_family_fractional_counts.tsv")

    if not out_dir.exists():
        out_dir.mkdir(parents=True)

    # main routine ------------------------------------------------------------
    res_path = align_to_pseudogenome(out_dir, pseudogenome_db, sample_name, query_type, query, query2, threads, paired_end, max_words)
    raw_uniq, raw_frac, raw_total = count_magicBLAST_result(res_path, perc_id, debug)
    elem_mapper = create_element_mapping(repnames_bedfile)
    uniq_counts = map_counts(raw_uniq, elem_mapper)
    frac_counts = map_counts(raw_frac, elem_mapper)
    total_counts = map_counts(raw_total, elem_mapper)
    fractional_counts = combine_counts(uniq_counts, frac_counts)  # fractional counts must be combined with unique counts before writing out
    write_output(fractional_counts, fractional_counts_outfile)
    write_output(total_counts, total_counts_outfile)
    write_output(uniq_counts, uniq_counts_outfile)

    if summarize_data:
        print("Writing class-level and family-level summarized counts...")
        class_total, family_total = summarize(total_counts)
        class_uniq, family_uniq = summarize(uniq_counts)
        class_fractional, family_fractional = summarize(fractional_counts)
        write_class_summary(class_total, class_total_counts_outfile)
        write_class_summary(class_uniq, class_uniq_counts_outfile)
        write_class_summary(class_fractional, class_fractional_counts_outfile)
        write_family_summary(family_total, family_total_counts_outfile)
        write_family_summary(family_uniq, family_uniq_counts_outfile)
        write_family_summary(family_fractional, family_fractional_counts_outfile)


if __name__ == '__main__':
    main()