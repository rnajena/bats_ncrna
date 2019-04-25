#!/usr/bin/env python3
# merge_gtf_ncbi.py
# SK

# same as merge_gtf_global_ids.py - but does not break on lines without gene_id and instead keeps them unchanged

import os
import sys


def error(string, error_type=1):
    sys.stderr.write('ERROR: ' + string + '\n')
    sys.exit(error_type)


def log(string, newline=False):
    if newline:
        sys.stderr.write('\n')
    sys.stderr.write('LOG: ' + string + '\n')

#####


class line_record:

    def __init__(self, gene_id, trans_id, exon_id, ftype, biotype, source, start, stop, strand, idx, line):
        self.gene_id = gene_id
        self.trans_id = trans_id
        self.exon_id = exon_id
        self.ftype = ftype
        self.biotype = biotype
        self.source = source
        self.start = start
        self.stop = stop
        self.strand = strand
        self.idx = idx
        self.line = line


class gene_record:

    def __init__(self, line_record):
        self.gene_id = line_record.gene_id
        self.source = line_record.source
        self.start = line_record.start
        self.stop = line_record.stop
        self.strand = line_record.strand
        self.lrecords = [line_record]
        self.biotype = line_record.biotype
        self.idxset = set([line_record.idx])

    def add_line(self, line_record):
        self.lrecords.append(line_record)
        assert line_record.idx not in self.idxset
        self.idxset.add(line_record.idx)


#####

log('==============================================================')
log(f'Started merging annotation files ...')


files_input = sys.argv[1:]

# check
files = []
for file in files_input:
    if not os.path.isfile(file):
        log(f'WARNING: Input file is missing: {file}')
    else:
        files.append(file)


sequence_names = {}
source_names = {}
feature_types = {}

# records data
# dict with all sequence names, each another dict with gene_ids
data = {}

gene_ids = []
trans_ids = []
exon_ids = []

warn_dupl = 0
warn_uids = 0

warn_overlaps = 0
three_plus_overlaps = 0
dele_overlaps_genes = 0
dele_overlaps_trans = 0
dele_overlaps_lines = 0


# collect lines that are not checked for overlaps and kept unchanged
keep_unchanged = []



#####################################


linectr = 0

# parse files
for file in files:

    with open(file) as infh:

        for line in infh:

            # skip comment / empty lines
            if line.startswith('#') or line == '\n':
                continue

            # status
            linectr += 1
            if linectr % 100 == 0:
                sys.stderr.write(f'LOG: Parsing lines: {linectr}\r')

            # global line index
            idx = linectr

            seqname, source, ftype, start, stop, _, strand, *_, desc = line.strip().split('\t')

            if seqname not in sequence_names:
                sequence_names[seqname] = 1
            else:
                sequence_names[seqname] += 1

            if ftype not in feature_types:
                feature_types[ftype] = 1
            else:
                feature_types[ftype] += 1

            # read tags
            try:
                tspl = desc.split(';')
                if tspl[-1] == '':
                    tspl = tspl[:-1]
                pairs = [t.strip().split(' ', 1) for t in tspl]
                tags = {t: v.strip('"') for t, v in pairs}
            except:
                error(f'Bad format of tags on line:\n{line}')


            # gather lines that have no gene_id and thus not checked and kept unchanged
            if 'gene_id' not in tags:
                keep_unchanged.append((seqname, int(start), line))
                continue


            # gather tags and check for duplicates
            gene_id = None
            trans_id = None
            exon_id = None

            # this is propagated to children transcripts/exons
            biotype = None

            if 'gene_id' in tags and 'transcript_id' not in tags and 'exon_id' not in tags:
                gene_id = tags['gene_id']

                if gene_id in gene_ids:
                    warn_dupl += 1
                    if warn_dupl <= 10:
                        log(f'WARNING: Duplicate gene_id: {gene_id} - in file: {file}')
                else:
                    gene_ids.append(gene_id)

                if 'gene_biotype' in tags:
                    biotype = tags['gene_biotype']


                if source not in source_names:
                    source_names[source] = 1
                else:
                    source_names[source] += 1



            elif 'gene_id' in tags and 'transcript_id' in tags and 'exon_id' not in tags:
                gene_id = tags['gene_id']
                trans_id = tags['transcript_id']

                if trans_id in trans_ids:
                    warn_dupl += 1
                    if warn_dupl <= 10:
                        log(f'WARNING: Duplicate transcript_id: {trans_id} - in file: {file}')
                else:
                    trans_ids.append(trans_id)



            elif 'gene_id' in tags and 'transcript_id' in tags and 'exon_id' in tags:
                gene_id = tags['gene_id']
                trans_id = tags['transcript_id']
                exon_id = tags['exon_id']

                if exon_id in exon_ids:
                    warn_dupl += 1
                    if warn_dupl <= 10:
                        log(f'WARNING: Duplicate exon_id: {exon_id} - in file: {file}')
                else:
                    exon_ids.append(exon_id)


            else:
                # other ftype, but has a gene_id
                # set gene_id and go on
                gene_id = tags['gene_id']

                log(f'ftype: {ftype} - gene_id: {gene_id}')



            # add this line record
            lrecord = line_record(gene_id, trans_id, exon_id, ftype, biotype, source, int(start), int(stop), strand, idx, line)

            if seqname not in data:
                data[seqname] = {gene_id: [lrecord]}
            else:
                # known sequence
                if gene_id in data[seqname]:
                    data[seqname][gene_id].append(lrecord)
                else:
                    data[seqname][gene_id] = [(lrecord)]

    # end parsing one file
sys.stderr.write(f'LOG: Parsing lines: {linectr}\r')
log('Parsed input gtf files.', newline=True)


###################################


# log stats
log(f'> Sequences: {len(sequence_names)}')
# for s, c in sorted(sequence_names.items(), key=lambda x: x[1], reverse=True):
#     log(f'{s}\t{c}')

log(f'> Feature types:')
for s, c in sorted(feature_types.items(), key=lambda x: x[1], reverse=True):
    log(f'{s}\t{c}')

log(f'> Annotation sources (on gene level):')
for s, c in sorted(source_names.items(), key=lambda x: x[1], reverse=True):
    log(f'{s}\t{c}')


log(f'> There were {len(keep_unchanged)} lines without gene_id. These are kept as is.')

#############################
# handle overlaps on all reference sequences

def check_if_overlap(rec1, rec2, thresh=0.5):
    ''' check if overlap as (% of shorter feature covered by larger feature) is > thresh (0.5)
    '''

    # make sure rec1 is 'left' of rec2
    if rec1.start > rec2.start:
        rec1, rec2 = rec2, rec1

    # check for any overlap
    if rec2.strand == rec1.strand and rec2.start <= rec1.stop:

        a1, a2 = rec1.start, rec1.stop
        b1, b2 = rec2.start, rec2.stop

        # overlap length
        # given: b1 >= a1
        ov = min(a2, b2) - b1 + 1

        alen = a2 - a1 + 1
        blen = b2 - b1 + 1

        # calculate overlap as (% of shorter feature covered by larger feature)
        ovpct = ov / min(alen, blen)

        #sys.stderr.write(f'ov check: {a1}-{a2}, {b1}-{b2}  =>  {ov} / min({alen}, {blen})  =>  {ovpct}\n')

        # decide by threshold
        if ovpct > thresh:
            return True

    # no overlap at all
    return False



def handle_overlap_exon(ovset):
    '''handle overlaps at exon level, remove transcripts
    '''

    global warn_overlaps, dele_overlaps_genes, dele_overlaps_trans, dele_overlaps_lines


    # trustworthiness index, later/higher is better, everything else gets -1337
    priority_of_source = ['mirdeep2', 'RNAmmer-1.2', 'gorap', 'tRNAscan-SE', 'mitos']


    # set of deleted lines to be returned
    delete_idxs = set()

    # gather all exons
    list_of_exons = []
    list_of_lrecs = []
    for gr in overlapset:
        for lr in gr.lrecords:
            list_of_lrecs.append(lr)
            if lr.ftype == 'exon':
                list_of_exons.append(lr)
    n_exons = len(list_of_exons)

    # get priorities
    prios = []
    for lr in list_of_exons:
        try:
            prio = priority_of_source.index(lr.source)
        except ValueError:
            prio = -1337
        prios.append(prio)

    # check exon line records for overlaps
    for nr1, lr1 in enumerate(list_of_exons):
        if lr1.idx in delete_idxs:
            continue
        for nr2 in range(nr1+1, n_exons):
            lr2 = list_of_exons[nr2]
            if lr2.idx in delete_idxs:
                continue
            if lr1.idx in delete_idxs:
                break
            
            # ignore internal overlaps of genes (multiple transcripts)
            if lr1.gene_id == lr2.gene_id:
                continue

            # handle overlapping exons by choosing best of two
            if check_if_overlap(lr1, lr2):
                
                warn_overlaps += 1
                decide_by_source = True


                # check biotype first
                if lr1.biotype == 'protein_coding':
                    if lr2.biotype != 'protein_coding':
                        del_lr = lr2
                        decide_by_source = False
                        
                elif lr2.biotype == 'protein_coding':
                    if lr1.biotype != 'protein_coding':
                        del_lr = lr1
                        decide_by_source = False


                # check tool priority second
                if decide_by_source:

                    p1, p2 = prios[nr1], prios[nr2]

                    # choose line record to delete
                    if p1 > p2:
                        del_lr = lr2
                    elif p2 > p1:
                        del_lr = lr1
                    else:
                        # log(f'Same priority exon overlap:\n{[lr1.gene_id, lr1.source, lr1.start, lr1.stop, lr1.strand]}' +
                        #     f'\n{[lr2.gene_id, lr2.source, lr2.start, lr2.stop, lr2.strand]}')
                        # choose longer feature:
                        if lr1.stop - lr1.start >= lr2.stop - lr2.start:
                            del_lr = lr2
                        else:
                            del_lr = lr1


                # add to deletes
                del_trans = del_lr.trans_id
                del_new_idxs = {lr.idx for lr in list_of_lrecs if lr.trans_id == del_trans}

                assert delete_idxs & del_new_idxs == set()
                delete_idxs |= del_new_idxs

                dele_overlaps_trans += 1
                dele_overlaps_lines += len(del_new_idxs)


    # check if any gene record is empty now
    for gr in overlapset:

        if len(gr.idxset) < 3:
            log(f'Warning: Gene record with < 3 lines: {gr.gene_id}')

        # check if only the gene line remains
        if (gr.idxset - delete_idxs) == set([gr.lrecords[0].idx]):
            # log(f'Gene record lost all transcripts: {[gr.gene_id, gr.source, gr.start, gr.stop, gr.strand]}')

            delete_idxs.add(gr.lrecords[0].idx)
            dele_overlaps_lines += 1
            dele_overlaps_genes += 1


    return delete_idxs


# UNUSED
def handle_overlap_gene(ovset):
    '''handle overlaps at gene level
    UNUSED FOR NOW
    '''

    global warn_overlaps, three_plus_overlaps, dele_overlaps_genes, dele_overlaps_lines

    if len(ovset) >= 3:
        three_plus_overlaps += 1

    warn_overlaps += 1
    if warn_overlaps <= 5:
        log(f'Overlap: {len(ovset)} gene records:')
        for gr in ovset:
            log(f'{[gr.gene_id, gr.source, gr.start, gr.stop, gr.strand]}')
    if warn_overlaps == 5:
        log(f'Only 5 overlaps are shown ...')

    # trustworthiness index, later/higher is better, everything else gets -1337
    priority_of_source = ['mirdeep2', 'gorap']

    prios = []
    for gr in ovset:
        try:
            prio = priority_of_source.index(gr.source)
        except ValueError:
            prio = -1337
        prios.append(prio)

    top = max(prios)
    bestl = [i for i, p in enumerate(prios) if p == top]

    if len(bestl) == 1:
        # one is the best, pick this
        picki = bestl[0]
    else:
        # choose longest of best
        lens = [gr.stop - gr.start + 1 for gr in ovset]
        picki = lens.index(max(lens))

    # return indices of unpicked rest
    delete_idxs = set()
    for gr in [gr for i, gr in enumerate(ovset) if i != picki]:
        dele_overlaps_genes += 1
        dele_overlaps_lines += len(gr.idxset)

        assert delete_idxs & gr.idxset == set()
        delete_idxs |= gr.idxset

    return delete_idxs


###########################

log('> Checking for overlaps ...')

deleted_idx_set = set()

for seq in data:
    # for all reference sequences

    dlist = data[seq]

    gene_records = []

    # # loop over records
    # # gather full gene_records
    # this_gr = None
    # for lrec in dlist:

    #     # same gene_id -> add to this gene record
    #     if this_gr is not None and this_gr.gene_id == lrec.gene_id:

    #         # propagate biotype:
    #         if this_gr.biotype is not None:
    #             lrec.biotype = this_gr.biotype

    #         this_gr.add_line(lrec)

    #     else:
    #         # new gene_id
    #         if this_gr is None:
    #             # start with first gene record
                
    #             this_gr = gene_record(lrec)

    #         else:
    #             # add finished gene record
    #             gene_records.append(this_gr)

    #             # then continue with next one
                
    #             this_gr = gene_record(lrec)

    # # end for all line records
    # # append last gene record
    # gene_records.append(this_gr)



    # build gene records
    for gid in dlist:

        lines_list = dlist[gid]

        grec = None
        gidx = None
        biotype = None
        for i, lrec in enumerate(lines_list):
            if lrec.ftype == 'gene':
                if grec is not None:
                    error(f'Double gene line in one gene_id set: {gene_id}')
                # start gene record
                grec = gene_record(lrec)
                gidx = i
                biotype = lrec.biotype
        
        if grec is None:
            error(f'gene_id set without gene line: {gene_id}')
        
        # add rest
        for i, lrec in enumerate(lines_list):
            if i == gidx:
                continue
            
            if biotype is not None:
                lrec.biotype = biotype
            grec.add_line(lrec)

        # add finished gene record
        gene_records.append(grec)


    # sort by start pos
    gene_records.sort(key=lambda gr: gr.start)


    # handle overlaps
    base = None
    overlapset = None


    # gather sets of gene records that overlap at gene level
    for gr in gene_records:

        # log(f'{len(gr.lrecords)}')

        if base is None:
            overlapset = [gr]
            base = gr
            continue

        # find all overlapping records
        # check strand and position
        if check_if_overlap(base, gr, 0):

            overlapset.append(gr)

        # not overlapping
        else:
            # handle finished set
            if len(overlapset) > 1:
                # handle overlaps at exon level
                deletes = handle_overlap_exon(overlapset)
                assert deleted_idx_set & deletes == set()
                deleted_idx_set |= deletes


            # start new set
            base = gr
            overlapset = [gr]

    # end for all gene records
    # handle final set
    if len(overlapset) > 1:
        deletes = handle_overlap_exon(overlapset)
        assert deleted_idx_set & deletes == set()
        deleted_idx_set |= deletes


    # delete overlapping line records
    # data[seq] = [lrec for lrec in dlist if lrec.idx not in deleted_idx_set]
    # NOT NEEDED delete idxs are global now


log(f'> There were {warn_overlaps} overlapping exons between gene records.')
if warn_overlaps > 0:
    log(f'> Removed {dele_overlaps_genes} gene records (lost all transcripts).')
    log(f'> Removed {dele_overlaps_trans} transcript records.')
    log(f'> Removed {dele_overlaps_lines} total lines.')
    


log(f'Sorting kept lines without gene_id ...')

keep_unchanged_dict = {}
for lrc in keep_unchanged:
    seq, start, line = lrc
    if seq in keep_unchanged_dict:
        keep_unchanged_dict[seq].append(lrc)
    else:
        keep_unchanged_dict[seq] = [lrc]
for seq in keep_unchanged_dict:
    keep_unchanged_dict[seq] = sorted(keep_unchanged_dict[seq], key=lambda x: x[1])

keep_seqs = list(keep_unchanged_dict.keys())


############################
# output
log('Outputting records ...')
for seq in sorted(data.keys()):
    for gid in data[seq]:
        for rec in data[seq][gid]:
            if not rec.idx in deleted_idx_set:
                sys.stdout.write(rec.line)
    if seq in keep_unchanged_dict:
        for rec in keep_unchanged_dict[seq]:
            sys.stdout.write(rec[2])
        keep_seqs.remove(seq)

# write out rest that may be on new contigs
if keep_seqs != []:
    for seq in keep_seqs:
        for rec in keep_unchanged_dict[seq]:
            sys.stdout.write(rec[2])


if warn_dupl > 0:
    log(f'> WARNING: There were {warn_dupl} duplicate ids! Check input files!')
if warn_uids > 0:
    log(f'> WARNING: There were {warn_uids} nonmatching ids! Check gene_id, transcript_id and exon_id of records!')
log('Done.')
