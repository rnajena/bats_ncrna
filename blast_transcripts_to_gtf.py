import os
import sys
import glob

# 	0		1		2		3		4		5		6	7	8		9	10		11		12	13
# qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen

gene_counter = 0
transcript_counter = 0
exon_counter = 0

def analyze(current_lncRNA, current_lncRNA_dict, out):
    global gene_counter
    global transcript_counter
    global exon_counter
    for strand in current_lncRNA_dict:

        for contig in current_lncRNA_dict[strand]:

            for transcript in current_lncRNA_dict[strand][contig]:

                for transcript_copy in current_lncRNA_dict[strand][contig][transcript]:

                    gene_counter += 1
                    transcript_counter += 1

                    min_exon_start = sys.maxsize
                    max_exon_end = 0
                    exon_list = []

                    for exon in sorted(current_lncRNA_dict[strand][contig][transcript][transcript_copy], key=lambda x: min(x.split('\t')[9], x.split('\t')[10])):
                        
                        exon_counter += 1

                        line_tab = exon.split('\t')

                        min_exon_start = min(min_exon_start, int(line_tab[9]), int(line_tab[10]))
                        max_exon_end = max(max_exon_end, int(line_tab[9]), int(line_tab[10]))
              
                        exon_id = "{tag} \"{species}{type}L{num}\"".format(tag='exon_id', species=spec.upper(), type='E', num=str(exon_counter).zfill(11))
                        exon_list.append((exon_id, '\t'.join([contig, 'blast', 'exon', str(min(int(line_tab[9]), int(line_tab[10]))), str(max(int(line_tab[9]), int(line_tab[10]))), '.', strand, '.'])))
            
                    description = '; '.join(["gene_name \"{}\"".format(current_lncRNA), "gene_source \"LNCipedia Version 5.2 high confidence set\"", "gene_biotype \"lncRNA\""])
                    gene_id = "{tag} \"{species}{type}L{num}\"".format(tag='gene_id', species=spec.upper(), type='G', num=str(gene_counter).zfill(11))
                    transcript_id = "{tag} \"{species}{type}L{num}\"".format(tag='transcript_id', species=spec.upper(), type='T', num=str(transcript_counter).zfill(11))
                    transcript_attributes =  '; '.join(["LNCipedia_name \"{}\"".format(line_tab[0].split('_')[0]), "LNCipedia_transcript_number \"{}\"".format(transcript), "transcript_copy \"{}\"".format(int(transcript_copy)+1)])

                    out.write('\t'.join([contig, 'blast', 'gene', str(min_exon_start), str(max_exon_end), '.', strand, '.']) + '\t' + '; '.join([gene_id, description]) + ';\n')
                    out.write('\t'.join([contig, 'blast', 'transcript', str(min_exon_start), str(max_exon_end), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, transcript_attributes, description]) + ';\n')

                    for e in exon_list:
                         out.write(e[1] + '\t' + '; '.join([gene_id, transcript_id, e[0], description]) + ';\n')


spec_blast = os.path.abspath(sys.argv[1])

spec = os.path.basename(spec_blast).split('.')[0]

out_gtf = spec_blast + '.gtf'
if(os.path.isfile(out_gtf)):
    os.remove(out_gtf)

current_lncRNA = ''
current_lncRNA_dict = {'+': {}, '-': {}}
gene_counter = 0

with open(spec_blast, 'r') as blast, open(out_gtf, 'a') as gtf:

    for line in blast:
        line_tab = line.split('\t')

        lncRNA = line_tab[0].split('_')[0].split(':')[0]
        transcript_number_lncipedia = line_tab[0].split('_')[0].split(':')[-1]
        exon_number_blast = line_tab[0].split('_')[-1]

        if(current_lncRNA != lncRNA):
            # new lncRNA
            # analyze old, contigwise
            analyze(current_lncRNA, current_lncRNA_dict, gtf)

            # reset
            current_lncRNA = lncRNA
            current_lncRNA_dict = {'+': {}, '-': {}}

        if (line_tab[9] < line_tab[10]):
            # plus
            if(line_tab[1] not in current_lncRNA_dict['+']):
                current_lncRNA_dict['+'][line_tab[1]] = {}
            if(transcript_number_lncipedia not in current_lncRNA_dict['+'][line_tab[1]]):
                current_lncRNA_dict['+'][line_tab[1]][transcript_number_lncipedia] = {}
            if(exon_number_blast not in current_lncRNA_dict['+'][line_tab[1]][transcript_number_lncipedia]):
                current_lncRNA_dict['+'][line_tab[1]][transcript_number_lncipedia][exon_number_blast] = []
            current_lncRNA_dict['+'][line_tab[1]][transcript_number_lncipedia][exon_number_blast].append(line)
        else:
            # minus
            if(line_tab[1] not in current_lncRNA_dict['-']):
                current_lncRNA_dict['-'][line_tab[1]] = {}
            if(transcript_number_lncipedia not in current_lncRNA_dict['-'][line_tab[1]]):
                current_lncRNA_dict['-'][line_tab[1]][transcript_number_lncipedia] = {}
            if(exon_number_blast not in current_lncRNA_dict['-'][line_tab[1]][transcript_number_lncipedia]):
                current_lncRNA_dict['-'][line_tab[1]][transcript_number_lncipedia][exon_number_blast] = []
            current_lncRNA_dict['-'][line_tab[1]][transcript_number_lncipedia][exon_number_blast].append(line)

    analyze(current_lncRNA, current_lncRNA_dict, gtf)