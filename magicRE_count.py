"""
Perform an alignment of reads to the pseudogenome
"""
from collections import defaultdict
from pathlib import Path
import argparse


def count(pseudogenome_sam):
    """Return Unique, fractional and total counts from parsing the pseudogenome"""
    print("Counting reads overlapped to pseudogenome...")
    uniq_counts = defaultdict(float)
    total_counts = defaultdict(float)
    frac_numerator = defaultdict(lambda: defaultdict(float))
    frac_counts = defaultdict(float)

    # calculate uniq counts and tally multi counts
    with open(pseudogenome_sam, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                continue  # skip headers
            else:
                l = line.strip().split("\t")
                read_id = l[0]
                rep_element = l[2]
                mapq = int(l[4])

                # tally total counts for all rep elemnets in the alignment
                total_counts[rep_element] += 1

                if mapq == 255:  # calculate unique counts
                    uniq_counts[rep_element] += 1
                else:            # tally the sub families that a read multimaps to
                    frac_numerator[read_id][rep_element] = 1  # every element gets numerator 1 (even if observed multiple times)

    # calculating fractional counts
    for read_id in frac_numerator.keys():
        N_subfamilies = len(frac_numerator[read_id].keys())
        for sub_family in frac_numerator[read_id].keys():
            frac_counts[sub_family] += (1 / N_subfamilies)  # frac count == 1 / N_subfamilies read aligns to

    return uniq_counts, frac_counts, total_counts


def create_element_mapping(repnames_bedfile):
    """Create a mapping of the element names to their classes and families"""
    print("Creating an element mapping...")
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
    """Create a dictionary of multi-mapped reads using the class > family > element hierarchy.
    Since the output from the count_multi function returns element names and their counts, 
    these element names must be re-associated with the classes and families they are part of. This 
    association is performed by the element_mapper.
    """
    print("Mapping sub-family counts back to their classes and families...")
    hierarchy = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for elem in multi_counts.keys():
        class_ = element_mapper[elem]['class']
        family = element_mapper[elem]['family']
        hierarchy[class_][family][elem] += multi_counts[elem]

    return hierarchy


def write_output(count_data, outfile):
    """Write count data to the outfile. Output truncates floats to ints."""
    print(f"Writing results to {outfile}...")
    with open(outfile, "w") as f:
        for class_ in count_data:
            for family in count_data[class_]:
                for elem in count_data[class_][family]:
                    count = int(count_data[class_][family][elem])
                    f.write(f"{class_}\t{family}\t{elem}\t{count}\n")


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
    parser = argparse.ArgumentParser(prog="count_from_pseudoalignment.py",
                                    usage="count_from_pseudoalignment.py file_prefix in.sam repnames.bed")
    parser.add_argument('file_prefix', help="Sample name to prefix the results files with")
    parser.add_argument('in_sam', help="SAM file of reads aligned to the pseudogenome")
    parser.add_argument('annotation_file', help="bedfile of repeat annotation file")
    args = parser.parse_args()

    file_prefix = args.file_prefix
    samfile = args.in_sam
    repnames_bed = args.annotation_file

    # Set up outfile names-----------------------------------------------------
    sam_parent = Path(samfile).parent
    uniq_out = Path(sam_parent, f"{file_prefix}_unique_counts.tsv")
    frac_out = Path(sam_parent, f"{file_prefix}_fractional_counts.tsv")
    total_out = Path(sam_parent, f"{file_prefix}_total_counts.tsv")
    uniq_family_out = Path(sam_parent, f"{file_prefix}_family_unique_counts.tsv")
    frac_family_out = Path(sam_parent, f"{file_prefix}_family_fractional_counts.tsv")
    total_family_out = Path(sam_parent, f"{file_prefix}_family_total_counts.tsv")
    uniq_class_out = Path(sam_parent, f"{file_prefix}_class_unique_counts.tsv")
    frac_class_out = Path(sam_parent, f"{file_prefix}_class_fractional_counts.tsv")
    total_class_out = Path(sam_parent, f"{file_prefix}_class_total_counts.tsv")

    uniq, frac, tot = count(samfile)
    elem_mapper = create_element_mapping(repnames_bed)
    uniq_counts = map_counts(uniq, elem_mapper)
    frac_counts = map_counts(frac, elem_mapper)
    total_counts = map_counts(tot, elem_mapper)

    write_output(uniq_counts, uniq_out)
    write_output(frac_counts, frac_out)
    write_output(total_counts, total_out)

    class_total, family_total = summarize(total_counts)
    class_uniq, family_uniq = summarize(uniq_counts)
    class_fractional, family_fractional = summarize(frac_counts)

    write_class_summary(class_total, total_class_out)
    write_class_summary(class_uniq, uniq_class_out)
    write_class_summary(class_fractional, frac_class_out)
    write_family_summary(family_total, total_family_out)
    write_family_summary(family_uniq, uniq_family_out)
    write_family_summary(family_fractional, frac_family_out)


if __name__ == '__main__':
    main()