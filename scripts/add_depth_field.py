#!/usr/bin/env python
## Adapted from a script by Martijn Derks

import sys
import gzip
import vcf
import argparse
import os
import numpy as np
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts

parser = argparse.ArgumentParser( description='Select candidate SV events')
parser.add_argument("-v", "--vcf_file", help="VCF_file with SV events", nargs=1)
parser.add_argument("-Q","--QUAL",help="Minimum Quality score",nargs=1,type=int)
parser.add_argument("-p","--Prefix",help="Prefix for output file",nargs=1,type=str)
# parser.add_argument("-bnd","--OutBND",help="Name of output VCF fir BNDs",nargs=1,type=str)


class Coverage():

	def CNV_VCF(self, vcf_file,Qual, Prefix):

		vcf_reader = vcf.Reader(filename = vcf_file, strict_whitespace=True)
		vcf_reader.infos['DEPTH'] = VcfInfo('DEPTH', vcf_field_counts['A'], 'String','Average coverage for Heteterozygous and Homozygous (reference) samples', source='devnull', version='devnull')
		vcf_writer = vcf.Writer(open(Prefix + "_DUP_DEL_INV.vcf", 'w'), vcf_reader)
		vcf_writer_BND = vcf.Writer(open(Prefix + "_BND.vcf", 'w'), vcf_reader)
		samples = vcf_reader.samples
		for record in vcf_reader:
			chr=record.CHROM
			type = record.INFO['SVTYPE']
			try: 
				imprecise = record.INFO['IMPRECISE']
				record.INFO['IMPRECISE'] = "TRUE"
			except:
				imprecise = False
				record.INFO['IMPRECISE'] = "FALSE"
			try: 
				event=record.INFO['EVENT']
			except:
				event="NA"
				record.INFO['EVENT'] = "NA"				

			if int(record.QUAL) < Qual:
				continue

			if type == "DUP" or type == "DEL" or type == "INV":
				het_ratios,hom_ref_ratios,hom_alt_ratios,het_dhd,hom_ref_dhd,hom_alt_dhd=[],[],[],[],[],[]

				for sample in record.get_hets():
					het_ratios.append(float(sample['DHBFC']))
					het_dhd.append(float(sample['DHFFC']))
				for sample in record.get_hom_refs():
					hom_ref_ratios.append(float(sample['DHBFC']))
					hom_ref_dhd.append(float(sample['DHFFC']))
				for sample in record.get_hom_alts():
					hom_alt_ratios.append(float(sample['DHBFC']))
					hom_alt_dhd.append(float(sample['DHFFC']))					
				avg_het,avg_hom_ref,avg_hom_alt,avg_het_dhd,avg_hom_ref_dhd,avg_hom_alt_dhd=np.mean(het_ratios), np.mean(hom_ref_ratios), np.mean(hom_alt_ratios),np.mean(het_dhd), np.mean(hom_ref_dhd), np.mean(hom_alt_dhd)		


				record.add_info('DEPTH', str(avg_hom_alt)+","+str(avg_het)+","+str(avg_hom_ref)+","+str(avg_hom_alt_dhd)+","+str(avg_het_dhd)+","+str(avg_hom_ref_dhd))
				vcf_writer.write_record(record)

			elif type == "BND":
				vcf_writer_BND.write_record(record)

			
				
if __name__ == '__main__':
	args = parser.parse_args()
	vcf_input = args.vcf_file[0]
	Qual = args.QUAL[0]
	Prefix = args.Prefix[0]
	C=Coverage()
	C.CNV_VCF(vcf_input,Qual,Prefix)			
	
