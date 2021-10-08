import sys
import vcf
import pandas as pd
import csv
import argparse



parser = argparse.ArgumentParser( description='Select candidate SV events')
parser.add_argument("-v", "--VCF", help="VCF_file with SV events", nargs=1, type=str)
parser.add_argument("-o","--OUTFILE",help="Output tsv file",nargs=1,type=str)
parser.add_argument("-s","--SAMPLES_TABLE",help="table with sample info. First column: sample; second column: bam file; third column: population",nargs=1,type=str)

args = parser.parse_args()

reader = vcf.Reader(filename = args.VCF[0])

samples_table = pd.read_csv(args.SAMPLES_TABLE[0], names = ["sample", "bam", "population"])

csvfile = open(args.OUTFILE[0], "w")
outfile = csv.writer(csvfile, delimiter='\t')

header = ["CHR", "POS", "TYPE", "SVLEN"]
populations = samples_table.population.unique()

population_headers = [p + "_freq" for p in populations]

complete_header = header + population_headers + ["SYMBOL", "GENEID", "ANNOTATION", "DHBFC_1/1", "DHBFC_0/1", "DHBFC_0/0", "DHFFC_1/1", "DHFFC_0/1", "DHFFC_0/0" ]
outfile.writerow(complete_header)

infos = reader.infos.keys()


for record in reader:
    row_to_write = [record.CHROM, record.POS, record.INFO.get("SVTYPE"), record.INFO.get("SVLEN")[0]]
    
    pop_gt_dict = {}
    for sample in record.samples:
        population = samples_table.loc[samples_table["sample"] == sample.sample, "population"].item()
        if population in pop_gt_dict:
            pop_gt_dict[population].append(sample["GT"])
        else:
            pop_gt_dict[population] = [sample["GT"]]

    for key, value in pop_gt_dict.items():
        points = 0
        for el in value:
            if el == "1/1":
                points = points + 2
            elif el == "1/0" or el == "0/1":
                points = points + 1
            elif el == "0/0":
                continue
        freq = points / (len(value)*2) *100
        row_to_write.append(freq)

    csqs = record.INFO.get("CSQ")
    
    ann,symbol,geneid=[],[],[]
    for csq in csqs:
        csq=csq.split("|")
        ann.append(csq[1])
        if csq[3] == "":
            csq[3] = "-"
        if csq[3] not in symbol:
            symbol.append(csq[3])
        if csq[4] not in geneid:	
            geneid.append(csq[4])
            
    row_to_write.append(",".join(symbol))
    row_to_write.append(",".join(geneid))
    row_to_write.append(",".join(ann))
    
    if "DEPTH" in infos:
        depth = record.INFO.get("DEPTH")
        [row_to_write.append(d) for d in depth]
    else:
        row_to_write.append("-")
    outfile.writerow(row_to_write)
    
csvfile.close()
