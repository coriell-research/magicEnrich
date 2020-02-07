from collections import defaultdict
from pathlib import Path
from multiprocessing import Pool


def parse_repeatmasker(annotation_file):
    """Extract information for repnames.bed file from repeatmasker input"""
    repeat_data = defaultdict(dd)
    name_idx, chrom_idx, start_idx, end_idx, class_fam_idx = 9, 4, 5, 6, 10
    rep_bed = Path(outdir, 'repnames.bed')

    with open(rep_bed, 'w') as outfile:
        with open(annotation_file, "r") as annot:
            for line in annot:
                if line.strip().startswith("SW"):
                    continue
                if line.strip().startswith("score"):
                    continue
                if line == "\n":
                    continue
            
                l = line.strip().split()
                name = l[name_idx]
                start = int(l[start_idx])
                end = int(l[end_idx])
                chrom = l[chrom_idx]
                class_family = l[class_fam_idx]
                try:
                    class_fam = class_family.split("/")
                    class_ = class_fam[0]
                    family = class_fam[1]
                except:
                    class_ = class_family    # make the class and family names the same
                    family = class_family    # for Simple_repeats, Low_complexity, etc.
                
                repeat_data[name]["start"].append(start)
                repeat_data[name]["end"].append(end)
                repeat_data[name]["count"].append(chrom)
                repeat_data[name]["chrom"].append(chrom)
    
    # convert the list of chromosomes to a count
    for k in repeat_data.keys():
        repeat_data[k]["count"] = len(repeat_data[k]["count"])
    
    return repeat_data


def parse_bedfile(annotation_file):
    """Extract information for repnames.bed file from bedfile input"""
    name_idx, chrom_idx, start_idx, end_idx, class_idx, family_idx = 3, 0, 1, 2, 4, 5
    repeat_data = defaultdict(dd)
    rep_bed = Path(outdir, 'repnames.bed')
    
    with open(rep_bed, 'w') as outfile:
        with open(annotation_file, "r") as annot:
            for line in annot:
                l = line.strip().split()
                name = l[name_idx]
                start = int(l[start_idx])
                end = int(l[end_idx])
                chrom = l[chrom_idx]
                class_ = l[class_idx]
                family = l[family_idx]
                
                repeat_data[name]["start"].append(start)
                repeat_data[name]["end"].append(end)
                repeat_data[name]["count"].append(chrom)
                repeat_data[name]["chrom"].append(chrom)
    
    # convert the list of chromosomes to a count
    for k in repeat_data.keys():
        repeat_data[k]["count"] = len(repeat_data[k]["count"])
    
    return repeat_data


def process_annotation(annotation_file, is_bed):
    """Extract the repeat name and information from annotation file"""
    print(f"Processing the annotation file : {annotation_file}...")

    if not is_bed:
        print(f"\t{annotation_file} is set to repeat masker format.")
        data = parse_repeatmasker(annotation_file)
    else:
        print(f"\t{annotation_file} is set to bedfile format.")
        data = parse_bedfile(annotation_file)
    
    print(f"Finished processing {annotation_file}.")
    
    return data


def generate_chromosome(genome, repeat_dict, repeat, flank_length, spacer_size):
    """Create a function that generates a 'chromosome' for the given repeat element by
    concatenating every genomic instance of the element along with its flanking region
    separated by a spacer
    """
    print(f"Creating pseudogenome for {repeat}...")
    pseudo_genome = []
    for i in range(repeat_dict[repeat]["count"]):
        chrom = repeat_dict[repeat]["chrom"][i]
        start_pos = max(repeat_dict[repeat]["start"][i] - flank_length, 0)
        end_pos = min(repeat_dict[repeat]["end"][i] + flank_length, len(genome[chrom].seq))

        spacer = "N" * spacer_size
        sequence = str(genome[chrom].seq[start_pos:end_pos])
        pseudo_genome.append(spacer)
        pseudo_genome.append(sequence)

    pseudo_genome_string = "".join(pseudo_genome)
    # create SeqRecord for writing to fasta - required to write fasta out
    record = SeqRecord(Seq(pseudo_genome_string, IUPAC.unambiguous_dna), 
                        id = repeat, 
                        name = "",
                        description = "")
    
    return record


def generate_pseudogenomes(genome_fasta, outfile, flank_length, spacer_size, repeat_dict, threads):
    """Generate a single pseudogenome file for all repeats in the annotation. 
    Each repeat constitutes one entry in the final fasta file.
    """
    print(f"Generating pseudogenome from {genome_fasta} using flank length = {flank_length} and spacer size = {spacer_size}...")
    print("This may take a long time...")
    g = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    # create lists of arguments for starmap
    repeats = [repeat for repeat in repeat_dict.keys()]
    genomes = [g for _ in range(len(repeats))]
    dicts = [repeat_dict for _ in range(len(repeats))]
    flanks = [flank_length for _ in range(len(repeats))]
    spacers = [spacer_size for _ in range(len(repeats))]

    # zip all arguments together into iterable
    args = list(zip(genomes, dicts, repeats, flanks, spacers))

    # Multiprocess each repeat instance
    with Pool(processes=threads) as pool:
        results = pool.starmap_async(generate_chromosome, args)
        records = results.get()
    
    print(f"Writing pseudogenome fasta file to : {outfile}...")
    SeqIO.write(records, outfile, "fasta")
    print("Done.")


def create_magicDB(pseudogenome_fasta):
    """Create a magic-BLAST db from the pseudogenome fasta file"""
    pass


def main():
    parser = argparse.ArgumentParser(description="Preparation of repetitive element psuedogenomes and repnames.bed file", usage="python magicRE_build.py path/to/repeat_masker.out path/to/genome.fasta --threads 8")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('annotationFile', help='The annotation file that contains repeat masker annotation for the genome of interest and may be downloaded at RepeatMasker.org. Example: hg38.fa.out. Custom bed files are also accepted.')
    parser.add_argument('genomeFasta', help='File name and path for genome of interest in fasta format. Example: /data/hg38.fa')
    parser.add_argument('--threads', default=1, type=int, metavar='1', help='Number of cores to be used for STAR index building.')
    parser.add_argument('--gapLength', default=200, type=int, metavar='200', help='Length of the spacer used to build repeat psuedogeneomes.')
    parser.add_argument('--flankLength', default=25, type=int, metavar='25', help='Length of the flanking region adjacent to the repeat element that is used to build repeat psuedogenomes. The flanking length should be set according to the length of your reads.')
    parser.add_argument('--isBed', dest='isBed', action='store_true', help="Is the annotation file a bed file? Bedfiles are also a compatible format. The file needs to be a tab seperated .bed with optional fields. Ex. format: chr\\tstart\\tend\\trepeat_name\\tclass\\tfamily.")
    parser.set_defaults(isBed=False)
    args = parser.parse_args()

    gap_length = args.gapLength
    flank_length = args.flankLength 
    annotation_file = args.annotationFile
    genome_fasta = args.genomeFasta
    threads = args.threads
    is_bed = args.isBed

    pseudo_genome_fasta = Path(setup_folder, "pseudogenome.fasta")
    create_setup_dir(setup_folder)
    repeat_data = process_annotation(annotation_file, setup_folder, is_bed)
    generate_pseudogenomes(genome_fasta, pseudo_genome_fasta, flank_length, gap_length, repeat_data, threads)
    create_magicDB(pseudo_genome_fasta)


if __name__ == '__main__':
    main()