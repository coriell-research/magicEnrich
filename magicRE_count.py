"""
Perform an alignment of reads to the pseudogenome,
Count unique, fractional, and total reads mapped.
Write output to tsv files.
"""
import argparse
import shlex
from collections import defaultdict
from pathlib import Path


def align_to_pseudogenome(pseudogenome_db, sample_name, read_type, read_handle, read_handle2, threads, paired_end, debug):
    """Align reads to the pseudogenome using magicBLAST"""
    tmp_dir = Path(Path(pseudogenome_db).parent, "__{sample_name}_tmp"))
    outfile = Path(tmp_dir, f"__{sample_name}_magicBLAST_results.tab")

    if read_type == "sra":
        align_cmd = f"magicblast -sra {read_handle} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -splice F -limit_lookup F > {outfile}"
    elif read_type == "fasta":
        if paired_end:
            align_cmd = f"magicblast -query {read_handle} -query_mate {read_handle2} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -splice F -limit_lookup F > {outfile}"
        else:
            align_cmd = f"magicblast -query {read_handle} -db {pseudogenome_db} -num_threads {threads} -outfmt tabular -no_unaligned -splice F -limit_lookup F > {outfile}"
    elif read_type == "fastq":
        if paired_end:
            align_cmd = f"magicblast -query {read_handle} -query_mate {read_handle2} -db {pseudogenome_db} -infmt fastq -num_threads {threads} -outfmt tabular -no_unaligned -splice F -limit_lookup F > {outfile}"
        else:
            align_cmd = f"magicblast -query {read_handle} -db {pseudogenome_db} -infmt fastq -num_threads {threads} -outfmt tabular -no_unaligned -splice F -limit_lookup F > {outfile}"
    else:
        raise TypeError("Unknown read_type: must be one of [sra, fasta, fastq]")
    
    # run the alignment
    subprocess.run(shlex.split(align_cmd), shell = True)

    if not debug:
        subprocess.run(f"rm -r {tmp_dir}")


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


def count_magicBLAST_result(magicBLAST_outfile, id_threshold):
    """Count the unique, total and fractional counts from the magicBLAST results"""
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


def map_counts(multi_counts, element_mapper):
    """Create a dictionary of reads using the class > family > element hierarchy.
    Since the output from the alignment step function returns element names and their counts, 
    these element names must be re-associated with the classes and families they are part of. 
    This association is performed by the element_mapper.
    """
    hierarchy = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for elem in multi_counts.keys():
        class_ = element_mapper[elem]['class']
        family = element_mapper[elem]['family']
        hierarchy[class_][family][elem] += multi_counts[elem]

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


def main():
    parser = argparse.ArgumentParser(description="Align reads to pseudogenome using magicBLAST and return counts of repeat elements")
    parser.add_argument('sampleName', help='The name of the sample to be processed.')
    parser.add_argument('setupDir', help="Path to the magicRE_Setup directory.")
    parser.add_argument('readType', help='Defines the input format to be aligned. One of [SRA, fasta, fastq]. Defines the input format to be aligned.', choices=['SRA', 'fasta', 'fastq'])
    parser.add_argument('--percentIdentity', default=90, type=float, help='Threshold at which to include mapped reads in the count data. Mapped reads must be >= percentIdentity to be included.')
    parser.add_argument('--threads', default=1, type=int, help='Number of threads to use for alignment.')
    parser.add_argument('--pairedEnd', dest='pairedEnd', action='store_true', help='Designate this option for paired-end sequencing.')
    parser.add_argument('--summarize', dest='summarize', action='store_true', help='In addition to the repeat name level output, produce collapsed counts at the class and family levels for each of the count types.')
    parser.add_argument('--debug', dest='debug', action='store_true', help='Select this option to prevent the removal of temporary files; useful for debugging')
    parser.set_defaults(pairedEnd=False, summarize=False, debug=False)
    parser.parse_args()


if __name__ == '__main__':
    main()