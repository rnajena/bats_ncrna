import os
import sys
from itertools import compress

# expected blast output format
# 	0		1		2		3		4		5		6	7	8		9	10		11		12	13
# qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen

global min_query_coverage
global max_transcript_range

def analyze(current_lncRNA_dict, out):
	for strand in current_lncRNA_dict:
		for contig in current_lncRNA_dict[strand]:

			current_transcript_exons = []
			current_transcript_range_on_contig = [0, 0]
			current_transcript_coverage = 0
			current_transcript_counter = 0

			len_sorted_exons = sorted(current_lncRNA_dict[strand][contig], key=lambda x: x[0], reverse=True)
			free_exons = [True for x in range(len(len_sorted_exons))]

			for e, exon in enumerate(len_sorted_exons):
				if(free_exons[e]):
					free_exons[e] = False
					current_transcript_exons.append(exon)
					current_transcript_range_on_contig = [min(exon[4], exon[5]), max(exon[4], exon[5])]
					current_transcript_coverage += exon[0]

					while(current_transcript_range_on_contig[1]-current_transcript_range_on_contig[0] <= max_transcript_range):
						not_overlapping_exons = []
						for f in list(compress(len_sorted_exons, free_exons)):
							for t in current_transcript_exons:
								if(not (f[1] > t[2] or f[2] < t[1])):
									break
							else:
								not_overlapping_exons.append(f)

						exons_left = [i for i in not_overlapping_exons if current_transcript_range_on_contig[0] >= max(i[4], i[5])]
						diff_exons_left_to_start = list(map(lambda x: (current_transcript_range_on_contig[0] - max(x[4], x[5]), x), exons_left))
						exons_right = [i for i in not_overlapping_exons if min(i[4], i[5]) >= current_transcript_range_on_contig[1]]
						diff_exons_rigth_to_end = list(map(lambda x: (min(x[4], x[5]) - current_transcript_range_on_contig[1], x), exons_right))

						if(not exons_left and not exons_right):
							break
						elif(not exons_left):
							next_exon = min(diff_exons_rigth_to_end, key=lambda x: x[0])[-1]
						elif(not exons_right):
							next_exon = min(diff_exons_left_to_start, key=lambda x: x[0])[-1]
						else:
							next_exon = min(min(diff_exons_left_to_start, key=lambda x: x[0]), min(diff_exons_rigth_to_end, key=lambda x: x[0]), key=lambda x: x[0])[-1]

						current_transcript_range_on_contig = [min(current_transcript_range_on_contig[0], next_exon[4], next_exon[5]), max(current_transcript_range_on_contig[1], next_exon[4], next_exon[5])]
						if(current_transcript_range_on_contig[1]-current_transcript_range_on_contig[0] > max_transcript_range):
							break

						free_exons[len_sorted_exons.index(next_exon)] = False
						current_transcript_exons.append(next_exon)
						current_transcript_coverage += next_exon[0]

					if(current_transcript_coverage >= min_query_coverage*exon[3]):

						for a in current_transcript_exons:
							out_line = a[-1].split('\t')
							out_line[0] = out_line[0] + '_' + str(current_transcript_counter)
							out.write('\t'.join(out_line))


					current_transcript_exons = []
					current_transcript_range_on_contig = [0, 0]
					current_transcript_coverage = 0
					current_transcript_counter += 1


# expects tab separated, by query name sorted blast output
blast_file = os.path.abspath(sys.argv[1])

if(len(sys.argv) > 2):
	min_query_coverage = float(sys.argv[2])
else:
	min_query_coverage = 0.7
if(len(sys.argv) > 3):
	max_transcript_range = int(sys.argv[3])
else:
	max_transcript_range = 500000

blast_filtered_file = '{}.transcripts'.format(blast_file)
if(os.path.isfile(blast_filtered_file)):
	os.remove(blast_filtered_file)

current_lncRNA = ''
current_lncRNA_dict = {'+': {}, '-': {}}

with open(blast_file, 'r') as bf, open(blast_filtered_file, 'a') as out:
	for line in bf:
		line_tab = line.split('\t')

		if(current_lncRNA != line_tab[0]):
			# new lncRNA
			# analyze old, contigwise
			analyze(current_lncRNA_dict, out)
			# reset
			current_lncRNA = line_tab[0]
			current_lncRNA_dict = {'+': {}, '-': {}}

		if (line_tab[9] < line_tab[10]):
			# plus
			if(line_tab[1] not in current_lncRNA_dict['+']):
				current_lncRNA_dict['+'][line_tab[1]] = []
			current_lncRNA_dict['+'][line_tab[1]].append((int(line_tab[7])-int(line_tab[6]), int(line_tab[6]), int(line_tab[7]), int(line_tab[8]), int(line_tab[9]), int(line_tab[10]), line))
		else:
			# minus
			if(line_tab[1] not in current_lncRNA_dict['-']):
				current_lncRNA_dict['-'][line_tab[1]] = []
			current_lncRNA_dict['-'][line_tab[1]].append((int(line_tab[7])-int(line_tab[6]), int(line_tab[6]), int(line_tab[7]), int(line_tab[8]), int(line_tab[9]), int(line_tab[10]), line))

	analyze(current_lncRNA_dict, out)
