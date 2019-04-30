import sys
import os
import functools

gtf = os.path.abspath(sys.argv[1])
if(len(sys.argv) > 2):
	min_overlapp = int(sys.argv[2])
else:
	min_overlapp = 10

gtf_out = gtf + '.merged.gtf'
if(os.path.isfile(gtf_out)):
    os.remove(gtf_out)

spec = os.path.basename(gtf).split('.')[0]

############### create single and multiple exon dict (stand -> exon num -> contigs)

single_exon_gtf_dict = {'+': {}, '-': {}}
multiple_exon_gtf_dict =  {'+': {}, '-': {}}

single_exon_stack_dict = {'+': {}, '-': {}}
multiple_exon_stack_dict = {'+': {}, '-': {}}

with open(gtf, 'r') as inp:
    
    last_transcript_info = {}
    last_exons = []
    last_strand = ''
    last_contig = ''

    for line in inp:
        line_tab = line.split('\t')
        
        entry_type = line_tab[2]
        
        if(entry_type == 'exon'):
            last_exons.append((int(line_tab[3]), int(line_tab[4]), last_transcript_info))
            last_strand = line_tab[6]
            last_contig = line_tab[0]

        elif(entry_type == 'transcript'):
            description = line_tab[-1][:-1] # get rid of \n 
            tspl = description.split(';')
            if tspl[-1] == '':
                tspl = tspl[:-1]
            pairs = [t.strip().split(' ',1) for t in tspl]
            tags = {t: v.strip('"') for t, v in pairs}
            
            last_transcript_info['LNCipedia_name'] = tags['LNCipedia_name']
            last_transcript_info['gene_name'] = tags['gene_name']
            last_transcript_info['transcript_id'] = tags['transcript_id']
            last_transcript_info['start'] = int(line_tab[3])
            last_transcript_info['end'] = int(line_tab[4])
        else:
            # gene
            if(len(last_exons) > 0):
                if(len(last_exons) == 1):
                    ## single_exons

                    if(last_contig not in single_exon_gtf_dict[last_strand]):
                        single_exon_gtf_dict[last_strand][last_contig] = []
                    if(last_contig not in single_exon_stack_dict[last_strand]):
                        single_exon_stack_dict[last_strand][last_contig] = []

                    for e in last_exons:
                        single_exon_gtf_dict[last_strand][last_contig].append(e)
                else:
                    ## multiple_exons

                    if(last_contig not in multiple_exon_gtf_dict[last_strand]):
                        multiple_exon_gtf_dict[last_strand][last_contig] = []
                    if(last_contig not in multiple_exon_stack_dict[last_strand]):
                        multiple_exon_stack_dict[last_strand][last_contig] = {}

                    transcript_exon_dict = {}
                    for e in last_exons:
                        transcript = (e[-1]['start'], e[-1]['end'], e[-1]['transcript_id'], e[-1]['LNCipedia_name'], e[-1]['gene_name'])
                        if (transcript not in transcript_exon_dict):
                            transcript_exon_dict[transcript] = []
                        transcript_exon_dict[transcript].append((e[0], e[1]))
                    for transcript in transcript_exon_dict:
                        l = list(transcript)
                        l.append(transcript_exon_dict[transcript])
                        multiple_exon_gtf_dict[last_strand][last_contig].append(tuple(l))
                
                last_exons = []
                last_transcript_info = {}
                last_strand = ''
    
    ## for the last line(
    if(len(last_exons) == 1):
        ## single_exons

        if(last_contig not in single_exon_gtf_dict[last_strand]):
            single_exon_gtf_dict[last_strand][last_contig] = []
        if(last_contig not in single_exon_stack_dict[last_strand]):
            single_exon_stack_dict[last_strand][last_contig] = []

        for e in last_exons:
            single_exon_gtf_dict[last_strand][last_contig].append(e)
    else:
        ## multiple_exons

        if(last_contig not in multiple_exon_gtf_dict[last_strand]):
            multiple_exon_gtf_dict[last_strand][last_contig] = []
        if(last_contig not in multiple_exon_stack_dict[last_strand]):
            multiple_exon_stack_dict[last_strand][last_contig] = {}

        transcript_exon_dict = {}
        for e in last_exons:
            transcript = (e[-1]['start'], e[-1]['end'], e[-1]['transcript_id'], e[-1]['LNCipedia_name'], e[-1]['gene_name'])
            if (transcript not in transcript_exon_dict):
                transcript_exon_dict[transcript] = []
            transcript_exon_dict[transcript].append((e[0], e[1]))
        for transcript in transcript_exon_dict:
            l = list(transcript)
            l.append(transcript_exon_dict[transcript])
            multiple_exon_gtf_dict[last_strand][last_contig].append(tuple(l))

    ## ) last line
   
            
############### merge single dict (stand -> contigs -> stacks)

for strand in single_exon_gtf_dict:
    for contig in single_exon_gtf_dict[strand]:

        start_sorted_exons = sorted(single_exon_gtf_dict[strand][contig], key=lambda x: (x[0], x[1]))
        merged = [False for i in range(0, len(start_sorted_exons))]

        for id, exon in enumerate(start_sorted_exons):
            if (not merged[id]):

                start = exon[0]
                stack = [exon]
                merged[id] = True
                right_overlapping_exon_ends = [ex[1] for i, ex in enumerate(start_sorted_exons) if not merged[i] and ex[0] <= exon[1] - min_overlapp and ex[1] > exon[1]]
                
                if (right_overlapping_exon_ends):
                    end = max(right_overlapping_exon_ends)

                    if(id != len(start_sorted_exons)):
                        for i, ex in enumerate(start_sorted_exons[id+1:]):
                            if(ex[1] <= end and not merged[i+id+1]):
                                merged[i+id+1] = True
                                stack.append(ex)
                            elif(ex[0] <= end - min_overlapp and ex[1] > end and not merged[i+id+1]):
                                end = ex[1]
                                merged[i+id+1] = True
                                stack.append(ex)

                else:
                    end = exon[1]
                    if(id != len(start_sorted_exons)):
                        for i, ex in enumerate(start_sorted_exons[id+1:]):
                            if(ex[1] <= end and not merged[i+id+1]):
                                merged[i+id+1] = True
                                stack.append(ex)
                
                single_exon_stack_dict[strand][contig].append(stack)

############### merge multiple dict (strand -> contig -> stacks of exons (regardless the transcript origin) -> union of stacks of exons with common transcript origin)

for strand in multiple_exon_gtf_dict:
    for contig in multiple_exon_gtf_dict[strand]:

        ## stacks of exons (regardless the transcript origin)
        start_sorted_exons = sorted([(e[0], e[1], trsct) for trsct in multiple_exon_gtf_dict[strand][contig] for e in trsct[-1]], key=lambda x: (x[0], x[1]))
        merged = [False for i in range(0, len(start_sorted_exons))]

        list_of_merged_stacks = []
        
        for id, exon in enumerate(start_sorted_exons):
            if (not merged[id]):

                start = exon[0]
                stack = [exon]
                merged[id] = True
                right_overlapping_exon_ends = [ex[1] for i, ex in enumerate(start_sorted_exons) if not merged[i] and ex[0] <= exon[1] - min_overlapp and ex[1] > exon[1]]
                
                if (right_overlapping_exon_ends):
                    end = max(right_overlapping_exon_ends)

                    if(id != len(start_sorted_exons)):
                        for i, ex in enumerate(start_sorted_exons[id+1:]):
                            if(ex[1] <= end and not merged[i+id+1]):
                                merged[i+id+1] = True
                                stack.append(ex)
                            elif(ex[0] <= end - min_overlapp and ex[1] > end and not merged[i+id+1]):
                                end = ex[1]
                                merged[i+id+1] = True
                                stack.append(ex)

                else:
                    end = exon[1]
                    if(id != len(start_sorted_exons)):
                        for i, ex in enumerate(start_sorted_exons[id+1:]):
                            if(ex[1] <= end and not merged[i+id+1]):
                                merged[i+id+1] = True
                                stack.append(ex)
                list_of_merged_stacks.append(stack)

        ## find transcripts ids in each exon stack 
        list_of_sets_of_trscp_ids = []
        for i in list_of_merged_stacks:
            list_of_sets_of_trscp_ids.append(set([j[-1][2] for j in i]))

        ## find common transcript ids by comparing each stack with all other stacks 
        intersecting_stacks = [set((i,j)) for i, first_stack in enumerate(list_of_sets_of_trscp_ids) for j in range(i+1, len(list_of_sets_of_trscp_ids)) if first_stack.intersection(list_of_sets_of_trscp_ids[j])]

        ## resolve transitivity for easier union of exon stacks
        if(len(intersecting_stacks) != 0):
            transitive_intersecting_stacks = [intersecting_stacks.pop()]
            while len(intersecting_stacks) > 0:
                curr_set = intersecting_stacks.pop()
                for i, other_set in enumerate(transitive_intersecting_stacks):
                    if other_set.intersection(curr_set):
                        new_set = other_set.union(curr_set)
                        transitive_intersecting_stacks[i] = new_set

                        for j, to_check_set in enumerate(transitive_intersecting_stacks):
                            if to_check_set.intersection(new_set) and to_check_set != new_set:
                                new_new_set = to_check_set.union(new_set)
                                transitive_intersecting_stacks.pop(j)
                                transitive_intersecting_stacks[i] = new_new_set
                        break
                else:
                    transitive_intersecting_stacks.append(curr_set)

        else:
            transitive_intersecting_stacks = [set(range(0, len(list_of_merged_stacks)))]

        ## actual merging
        for union in transitive_intersecting_stacks:
            merged_exons = []
            for i in union:
                for j in list_of_merged_stacks[i]:
                    merged_exons.append(j[-1])

            genes = set([exon[4] for exon in merged_exons])
            if(len(genes) == 1):
                if('lnc' not in multiple_exon_stack_dict[strand][contig]):
                    multiple_exon_stack_dict[strand][contig]['lnc'] = []
                multiple_exon_stack_dict[strand][contig]['lnc'].append(merged_exons)
            else:
                if('hotspot' not in multiple_exon_stack_dict[strand][contig]):
                    multiple_exon_stack_dict[strand][contig]['hotspot'] = []
                multiple_exon_stack_dict[strand][contig]['hotspot'].append(merged_exons)

############### write merged gtf

normal_gene_counter = 0
normal_transcript_counter = 0
normal_exon_counter = 0
hot_spot_gene_counter = 0
hot_spot_transcript_counter = 0
hot_spot_exon_counter = 0

#### single exons

for strand in single_exon_stack_dict:
    for contig in  single_exon_stack_dict[strand]:
        for stack in single_exon_stack_dict[strand][contig]:
            
            start = min(stack, key=lambda x: x[0])[0]
            end = max(stack, key=lambda x: x[1])[1]
            
            genes = set([gene[-1]['gene_name'] for gene in stack])

            if(len(genes) == 1):
                # one gene, one or multiple transcripts
                nc_class = 'L'

                gene_description = '; '.join(["gene_name \""+stack[0][-1]['gene_name']+"\"", "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA\""])
                normal_gene_counter += 1
                gene_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='gene_id', species=spec.upper(), type='G', cls=nc_class, num=str(normal_gene_counter).zfill(11))
                
                with open(gtf_out, 'a') as outp:
                    outp.write('\t'.join([contig, 'blast', 'gene', str(start), str(end), '.', strand, '.']) + '\t' + '; '.join([gene_id, gene_description]) + ';\n')

                for s in stack:
                    transcript_description = '; '.join(["LNCipedia_name \""+s[-1]['LNCipedia_name']+"\"", "LNCipedia_transcript_number \""+s[-1]['LNCipedia_name'].split(':')[-1]+"\"", "gene_name \""+s[-1]['gene_name']+"\"", "gene_source \"LNCipedia Version 5.2 high confidence set\"", "gene_biotype \"lncRNA\""])
                    normal_transcript_counter += 1
                    transcript_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='transcript_id', species=spec.upper(), type='T', cls=nc_class, num=str(normal_transcript_counter).zfill(11))

                    exon_description = '; '.join(["gene_name \""+s[-1]['gene_name']+"\"", "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA\""])
                    normal_exon_counter += 1
                    exon_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='exon_id', species=spec.upper(), type='E', cls=nc_class, num=str(normal_exon_counter).zfill(11))

                    with open(gtf_out, 'a') as outp:
                        outp.write('\t'.join([contig, 'blast', 'transcript', str(s[0]), str(s[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, transcript_description]) + ';\n')
                        outp.write('\t'.join([contig, 'blast', 'exon', str(s[0]), str(s[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, exon_id, exon_description]) + ';\n')
            else:
                # hotspot, multiple genes
                nc_class = 'H'

                hot_spot_gene_counter += 1

                description = '; '.join(["gene_name \"lncRNAhotspot-{}\"".format(hot_spot_gene_counter), "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA_hot_spot\""])
                hot_spot_description = "hot_spot_genes \"" + ','.join([s[-1]['LNCipedia_name']+'('+str(s[0])+'-'+str(s[1])+')' for s in stack]) + "\""
                
                gene_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='gene_id', species=spec.upper(), type='G', cls=nc_class, num=str(hot_spot_gene_counter).zfill(11))
                with open(gtf_out, 'a') as outp:
                    outp.write('\t'.join([contig, 'blast', 'gene', str(start), str(end), '.', strand, '.']) + '\t' + '; '.join([gene_id, description, hot_spot_description]) + ';\n')
                
                for s in stack:
                    hot_spot_transcript_counter += 1
                    transcript_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='transcript_id', species=spec.upper(), type='T', cls=nc_class, num=str(hot_spot_transcript_counter).zfill(11))
                    
                    hot_spot_exon_counter += 1
                    exon_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='exon_id', species=spec.upper(), type='E', cls=nc_class, num=str(hot_spot_exon_counter).zfill(11))

                    with open(gtf_out, 'a') as outp:
                        outp.write('\t'.join([contig, 'blast', 'transcript', str(s[0]), str(s[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, description]) + ';\n')
                        outp.write('\t'.join([contig, 'blast', 'exon', str(s[0]), str(s[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, exon_id, description]) + ';\n')
### multiple exons

for strand in multiple_exon_stack_dict:
    for contig in multiple_exon_stack_dict[strand]:
        if('lnc' in multiple_exon_stack_dict[strand][contig]):
            for lnc_stack in multiple_exon_stack_dict[strand][contig]['lnc']:

                start = min(lnc_stack, key=lambda x: x[0])[0]
                end = max(lnc_stack, key=lambda x: x[1])[1]
                
                nc_class = 'L'

                gene_description = '; '.join(["gene_name \""+lnc_stack[0][4]+"\"", "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA\""])
                normal_gene_counter += 1
                gene_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='gene_id', species=spec.upper(), type='G', cls=nc_class, num=str(normal_gene_counter).zfill(11))
                
                with open(gtf_out, 'a') as outp:
                    outp.write('\t'.join([contig, 'blast', 'gene', str(start), str(end), '.', strand, '.']) + '\t' + '; '.join([gene_id, gene_description]) + ';\n')

                for s in lnc_stack:
                    transcript_description = '; '.join(["LNCipedia_name \""+s[3]+"\"", "LNCipedia_transcript_number \""+s[3].split(':')[-1]+"\"", "gene_name \""+s[4]+"\"", "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA\""])
                    normal_transcript_counter += 1
                    transcript_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='transcript_id', species=spec.upper(), type='T', cls=nc_class, num=str(normal_transcript_counter).zfill(11))
                    
                    with open(gtf_out, 'a') as outp:
                        outp.write('\t'.join([contig, 'blast', 'transcript', str(s[0]), str(s[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, transcript_description]) + ';\n')
                    
                    for exon in s[-1]:
                        exon_description = '; '.join(["gene_name \""+s[4]+"\"", "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA\""])
                        normal_exon_counter += 1
                        exon_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='exon_id', species=spec.upper(), type='E', cls=nc_class, num=str(normal_exon_counter).zfill(11))

                        with open(gtf_out, 'a') as outp:
                            outp.write('\t'.join([contig, 'blast', 'exon', str(exon[0]), str(exon[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, exon_id, exon_description]) + ';\n')
        
        if('hotspot' in multiple_exon_stack_dict[strand][contig]):
            for hs_stack in multiple_exon_stack_dict[strand][contig]['hotspot']:
                
                start = min(hs_stack, key=lambda x: x[0])[0]
                end = max(hs_stack, key=lambda x: x[1])[1]
                
                nc_class = 'H'

                hot_spot_gene_counter += 1

                description = '; '.join(["gene_name \"lncRNAhotspot-{}\"".format(hot_spot_gene_counter), "gene_source \"LNCipedia_Version_5.2_high_confidence_set\"", "gene_biotype \"lncRNA_hot_spot\""])
                hot_spot_description = "hot_spot_genes \"" + ','.join([s[3]+'('+','.join([str(e[0])+'-'+str(e[1]) for e in sorted(s[-1], key=lambda x: (x[0], x[1]))])+')' for s in hs_stack]) + "\""
                gene_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='gene_id', species=spec.upper(), type='G', cls=nc_class, num=str(hot_spot_gene_counter).zfill(11))

                with open(gtf_out, 'a') as outp:
                    outp.write('\t'.join([contig, 'blast', 'gene', str(start), str(end), '.', strand, '.']) + '\t' + '; '.join([gene_id, description, hot_spot_description]) + ';\n')

                for trsct in hs_stack:
                    hot_spot_transcript_counter += 1
                    transcript_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='transcript_id', species=spec.upper(), type='T', cls=nc_class, num=str(hot_spot_transcript_counter).zfill(11))

                    with open(gtf_out, 'a') as outp:
                        outp.write('\t'.join([contig, 'blast', 'transcript', str(trsct[0]), str(trsct[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, description]) + ';\n')

                    for exon in trsct[-1]:
                        hot_spot_exon_counter += 1
                        exon_id = "{tag} \"{species}{type}{cls}{num}\"".format(tag='exon_id', species=spec.upper(), type='E', cls=nc_class, num=str(hot_spot_exon_counter).zfill(11))

                        with open(gtf_out, 'a') as outp:
                            outp.write('\t'.join([contig, 'blast', 'exon', str(exon[0]), str(exon[1]), '.', strand, '.']) + '\t' + '; '.join([gene_id, transcript_id, exon_id, description]) + ';\n')
