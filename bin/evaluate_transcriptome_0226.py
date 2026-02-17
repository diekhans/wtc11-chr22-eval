import pipettor
import argparse

# BED_READS_FILE = 'WTC11.ENCFF370NFS.chr22.genomealigned.bed'
# BED_READS_INTERVALS_FILE = 'WTC11.ENCFF370NFS.chr22.genomealigned.intervals.bed'  # will be created
# ANNOT_FILE = 'gencode.v38.annotation.chr22.gtf'


def parse_args():
    parser = argparse.ArgumentParser(description='''for evaluating transcriptomes compared to raw reads and reference annotation''')
    parser.add_argument('--gtf', required=True, 
                        help="gtf annotation file")
    parser.add_argument('--bed_reads', required=True, 
                        help="bed file of raw reads used to run transcriptome. can be generated with bedtools bamtobed")
    parser.add_argument('--isoforms', required=True, 
                        help="bed file of isoforms that you want to evaluate")
    parser.add_argument('--se_end_window', default=100, 
                        help="window in bp to compare ends of single exon isoforms")
    parser.add_argument('-o', '--output', required=True, help='output prefix')
  
    args = parser.parse_args()
    return args


def get_chromtoint(file):
    totreads = 0
    chromtoint = {}
    for line in open(file): 
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        if chrom not in chromtoint: chromtoint[chrom] = []
        chromtoint[chrom].append((start, end))
        totreads += 1
    return chromtoint, totreads

def get_regions(chromtoint, outfilename):
    totregions= 0
    out = open(outfilename, 'w')
    for chrom in chromtoint:
        intlist = sorted(chromtoint[chrom])
        newints = []
        laststart, lastend = 0, 0
        for s, e in intlist:
            if s > lastend:
                if lastend != 0:
                    newints.append((laststart, lastend))
                laststart = s
            if e > lastend:
                lastend = e
        newints.append((laststart, lastend))
        for s, e in newints:     
            out.write(f'{chrom}\t{s}\t{e}\n')   
            totregions += 1 
    out.close()
    return totregions

def get_intersect_count(filea, fileb):
    c = 0
    dr = pipettor.DataReader()
    pipettor.run([('bedtools', 'intersect', '-f', '0.5', '-u', '-a', filea, '-b', fileb)], stdout=dr)
    for line in dr.data.rstrip('\n').split('\n'):
        c += 1
    return c

def extract_sj_info(file, se_end_window):
    read_sjc, read_se_ends = {}, {}
    for line in open(file): ##bed reads file
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        cs = chrom #(chrom, strand) #can't trust read strand
        if cs not in read_sjc:
            read_sjc[cs] = {}
            read_se_ends[cs] = {}
        if len(introns) > 0:
            if introns not in read_sjc[cs]: read_sjc[cs][introns] = 0
            read_sjc[cs][introns] += 1
        else:
            roundedends = (se_end_window * round(start/se_end_window), \
                            se_end_window * round(end/se_end_window))
            if roundedends not in read_se_ends[cs]: read_se_ends[cs][roundedends] = 0
            read_se_ends[cs][roundedends] += 1
    return read_sjc, read_se_ends



def calculate_sjc_metrics(found_sjc, read_sjc):    

    found_subsets = {}
    for cs in found_sjc:
        found_subsets[cs] = set()
        for sjc in found_sjc[cs]:
            for slen in range(len(sjc)-1, 0, -1):
                for i in range(0, len(sjc)-slen+1):
                    found_subsets[cs].add(sjc[i:i+slen])

    tot_sjc, sup_sjc = 0, 0
    subset_sjc = 0
    for cs in read_sjc:
        if cs in found_sjc:
            for sjc in read_sjc[cs]:
                if sjc in found_sjc[cs]:
                    sup_sjc += read_sjc[cs][sjc]
                elif sjc in found_subsets[cs]:
                    subset_sjc += read_sjc[cs][sjc]
                tot_sjc += read_sjc[cs][sjc]
        else:
            for sjc in read_sjc[cs]:
                tot_sjc += read_sjc[cs][sjc]
    
    return tot_sjc, sup_sjc, subset_sjc,

def calculate_se_metrics(found_se_ends, read_se_ends):
    tot_se, sup_se = 0, 0
    for cs in read_se_ends:
        if cs in found_se_ends:
            for se in read_se_ends[cs]:
                if se in found_se_ends[cs]:
                    sup_se += read_se_ends[cs][se]
                tot_se += read_se_ends[cs][se]
        else:
            for se in read_se_ends[cs]:
                tot_se += read_se_ends[cs][se]

    return tot_se, sup_se


def process_gtf(gtf):
    transcripttoexons = {}
    for line in open(gtf):  # gtf reference file
        if line[0] != '#':
            line = line.split('\t')
            if line[2] == 'exon':
                chrom, strand = line[0], line[6]
                start, end = int(line[3]), int(line[4])
                tname = line[-1].split('transcript_id "')[1].split('"')[0]
                if (chrom, strand) not in transcripttoexons: transcripttoexons[(chrom, strand)] = {}
                if tname not in transcripttoexons[(chrom, strand)]: transcripttoexons[(chrom, strand)][tname] = []
                transcripttoexons[(chrom, strand)][tname].append((start, end))
    return transcripttoexons

def categorize_reference_isoforms(transcripttoexons):
    refjuncs, refjuncchains, refseends = {}, {}, {}
    for chrom, strand in transcripttoexons:
        refjuncs[(chrom, strand)] = set()
        refjuncchains[(chrom, strand)] = set()
        refseends[(chrom, strand)] = set()
        for tname in transcripttoexons[(chrom, strand)]:
            exons = sorted(transcripttoexons[(chrom, strand)][tname])
            if len(exons) > 1:
                introns = [(exons[x][1], exons[x+1][0]-1) for x in range(len(exons)-1)]  # remove 1 from exon start coord to match bed
                refjuncchains[(chrom, strand)].add(tuple(introns))
                #print(tuple(introns))
                refjuncs[(chrom, strand)].update(set(introns))
            else:
                refseends[(chrom, strand)].add(exons[0])
    return refjuncs, refjuncchains, refseends


def categorize_isoforms_relative_to_reference(isoform_file, se_end_window, refjuncchains, refjuncs, refseends):
    fsm, ism, nic, nnc, sem, sen, tot = 0, 0, 0, 0, 0, 0, 0
    for line in open(isoform_file):  # bed isoforms file
        line = line.rstrip().split('\t')
        chrom, strand = line[0], line[5]
        iso, start, end = line[3], int(line[1]), int(line[2])
        esizes, estarts = [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
        exons = [(start+estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))]
        introns = tuple([(exons[x][1], exons[x+1][0]) for x in range(len(exons)-1)])
        tot += 1
        if len(introns) > 0:
            if introns in refjuncchains[(chrom, strand)]:
                fsm += 1
            else:
                isISM = False
                myjuncstring = str(introns)[1:-1]
                for juncchain in refjuncchains[(chrom, strand)]:
                    if myjuncstring in str(juncchain):
                        isISM = True
                        ism += 1
                        break
                if not isISM:
                    allFound = True
                    for j in introns:
                        if j not in refjuncs[(chrom, strand)]:
                            allFound = False
                            break
                    if allFound: nic += 1
                    else: nnc += 1
        else:
            endsMatch = False
            for refstart, refend in refseends[(chrom, strand)]:
                if abs(start-refstart) < se_end_window and abs(end-refend) < se_end_window:
                    endsMatch = True
                    break
            if endsMatch: sem += 1
            else: sen += 1
    return tot, fsm, ism, nic, nnc, sem, sen

def print_section_data(section_names, section_values, outfile):
    for i in range(len(section_names)):
        outfile.write(section_names[i] + '\t' + str(section_values[i]) + '\n')


if __name__ == '__main__':
    args = parse_args()
    # EVALUATE GENIC REGIONS FOUND RELATIVE TO THOSE REPRESENTED IN READS

    chromtoint, totreads = get_chromtoint(args.bed_reads) ##bed reads file
    bed_reads_intervals = args.output + '_bed_reads_intervals.bed'
    totregions = get_regions(chromtoint, bed_reads_intervals)
    
    chromtoint, _ = get_chromtoint(args.isoforms)
    isoform_intervals = args.output + '_isoform_intervals.bed'
    get_regions(chromtoint, isoform_intervals)
   
    foundregions = get_intersect_count(bed_reads_intervals, isoform_intervals)
    genicreads = get_intersect_count(args.bed_reads, isoform_intervals)

    out = open(args.output + '.tsv', 'w')
    out.write('\t'.join(['metric', args.output]) + '\n')
    section_names = ['total_read_regions', 'found_read_regions', 'total_reads', 'reads_overlapping_found_regions']
    section_values = [totregions, foundregions, totreads, genicreads]
    print_section_data(section_names, section_values, out)

    # EVALUATE SPLICE JUNCTION CHAINS FOUND RELATIVE TO THOSE REPRESENTED IN READS
    # EVALUATE SINGLE EXON READS/ISOFORMS SEPARATELY - round ends to ROUND_ENDS_WINDOW, check match

    read_sjc, read_se_ends = extract_sj_info(args.bed_reads, args.se_end_window)
    found_sjc, found_se_ends = extract_sj_info(args.isoforms, args.se_end_window)
    tot_sjc, sup_sjc, subset_sjc  = calculate_sjc_metrics(found_sjc, read_sjc)
    tot_se, sup_se = calculate_se_metrics(found_se_ends, read_se_ends)

    section_names = ['total_SJC_from_reads', 'SJC_reported_in_isoforms', 'reads_where_SJC_is_subset_of_isoform', 'total_SE_reads', 'SE_reported_in_isoforms']
    section_values = [tot_sjc, sup_sjc, subset_sjc, tot_se, sup_se]
    print_section_data(section_names, section_values, out)

    # EVALUATE TRANSCRIPTS COMPARED TO ANNOTATION
    # LABEL EACH TRANSCRIPT AS FSM, ISM, NIC, NNC
    # Single exon handled separately again, using SINGLE_EXON_END_WINDOW

    transcripttoexons = process_gtf(args.gtf)
    refjuncs, refjuncchains, refseends = categorize_reference_isoforms(transcripttoexons)

    tot, fsm, ism, nic, nnc, sem, sen = categorize_isoforms_relative_to_reference(args.isoforms, args.se_end_window, refjuncchains, refjuncs, refseends)
    section_names = ['total_isoforms', 'FSM_isoforms', 'ISM_isoforms', 'NIC_isoforms', 'NNC_isoforms', 'SE_match_annot_isoforms', 'SE_not_supported_isoforms']
    section_values = tot, fsm, ism, nic, nnc, sem, sen
    print_section_data(section_names, section_values, out)

    section_names = [
        'percent_of_regions_where_reads_align_with_isoforms_reported', 
        'percent_of_spliced_reads_represented_by_an_isoform',
        'percent_of_isoforms_that_are_single_exon',
        'percent_of_single_exon_isoforms_that_match_annot',
        'FSM_pct_of_spliced_isoforms',
        'ISM_pct_of_spliced_isoforms',
        'NIC_pct_of_spliced_isoforms',
        'NNC_pct_of_spliced_isoforms',
    ]
    section_values = [
        round(100*foundregions/totregions, 2),
        round(100*(sup_sjc+subset_sjc)/tot_sjc, 2),
        round(100*(sem+sen)/tot, 2),
        round(100*sem/(sem+sen), 2),
        round(100*fsm/(fsm+ism+nic+nnc), 2),
        round(100*ism/(fsm+ism+nic+nnc), 2),
        round(100*nic/(fsm+ism+nic+nnc), 2),
        round(100*nnc/(fsm+ism+nic+nnc), 2),
    ]
    print_section_data(section_names, section_values, out)

    out.close()
    pipettor.run([('rm', bed_reads_intervals, isoform_intervals)])


