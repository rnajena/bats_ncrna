#!/usr/bin/env python3
# make_mitos_nice_again.py
# ML

import os

# species = ['dro', 'ehe', 'mbr', 'mda', 'mlu', 'pal', 'ppa', 'pva', 'rae', 'rfe', 'har']
species = ['ehe_rearranged']

# in
mitos_results = '/data/prostlocal2/projects/mh_bats_ncrna_annotation/2018/mitos/results'

# out
mitos_nice_dir = '/data/prostlocal2/projects/mh_bats_ncrna_annotation/2018/playground/mitos'

for spec in species:

	gff =  "{}/{}/{}.gff".format(mitos_results, spec, spec.split('_')[0])

	spec = spec.split('_')[0]


	out_dir = "{}/{}".format(mitos_nice_dir, spec)

	out_gtf = "{}/{}.gtf".format(out_dir, spec)

	if(os.path.isfile(out_gtf)):
		os.remove(out_gtf)

	rename_prefix = "{}_mito".format(spec)

	if(not os.path.isdir(out_dir)):
		os.makedirs(out_dir)

	gff_table = []
	with open(gff, 'r') as gf:
		for line in gf:
			line_tab = line.split("\t")
			line_tab[0] = rename_prefix
			gff_table.append(line_tab[:-1])

	features_types = [('gene', 'G', 'gene_id'), ('transcript', 'T', 'transcript_id'), ('exon', 'E', 'exon_id')]
	with open(out_gtf, 'a') as o_gt:
		for i, line in enumerate(gff_table):

			ids = ["{tag} \"{species}{type}O{num}\"".format(tag=f[2], species=spec.upper(), type=f[1], num=str(i+1).zfill(11)) for f in features_types]
			
			for j, f in enumerate(features_types):

				description_ids = '; '.join(ids[:j+1])
				if(line[2] == 'gene'):
					description_other = '; '.join(["gene_name \"{}\"".format(line[3]), "gene_source \"mitos\"", "gene_biotype \"{}\"".format('protein_coding'), "mitos_annotation_tool \"{}\"".format(line[1])])
				else:
					description_other = '; '.join(["gene_name \"{}\"".format(line[3]), "gene_source \"mitos\"", "gene_biotype \"{}\"".format(line[2]), "mitos_annotation_tool \"{}\"".format(line[1])])
				if (line[1] == 'mitos'):
					description_mitos_quality = '; '.join(["mitos_q_value \"{}\";".format(line[6])])
				else:
					description_mitos_quality = '; '.join(["mitos_e_value \"{}\";".format(line[6])])
				attributes = "\t".join([rename_prefix, 'mitos', f[0], line[4], line[5], '.', line[7], '.', '; '.join([description_ids, description_other, description_mitos_quality])])
				o_gt.write(attributes + '\n')
