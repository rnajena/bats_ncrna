#!/usr/bin/env python3
# format_ncbi.py
# SK

import sys


def error(string, error_type=1):
    sys.stderr.write(f'ERROR: {string}\n')
    sys.exit(error_type)


def log(string, newline_before=False):
    if newline_before:
        sys.stderr.write('\n')
    sys.stderr.write(f'LOG: {string}\n')

#####

infile = sys.argv[1]
log('==============================================================')
log(f'Started formatting NCBI annotation {infile}')

gfftoggle = False
gtftoggle = False


abbr = infile.rsplit('/',1)[-1].split('.',1)[0].upper()
NID = '_NCBI_ID_'

if len(abbr) != 3:
    error(f'Could not infer bat abbreviation: {infile}')

num_lines = 0
del_lines = {}
warn_types = []

####

def print_line(lr):
    # print line in gtf format
    global num_lines
    field9 = '; '.join([f'{l} "{r}"' for l, r in lr[-1]]) + ';'
    print(*lr[:-1], field9, sep='\t')
    num_lines += 1


def get_ID(lr):
    tags = lr[-1]
    ID = [v for t, v in tags if t == 'ID']
    if ID == []:
        error(f'No ID in line: {lr}')
    return ID[0]


def get_Parent(lr):
    tags = lr[-1]
    par = [v for t, v in tags if t == 'Parent']
    if par == []:
        return None
    return par[0]


def handle_gene_rec(gr):

    global del_lines

    maintype = gr[0][2]

    ID = get_ID(gr[0])
    gene_id = abbr + 'G' + NID + ID


    if maintype in ['gene']:

        # add tags to all lines
        # gene line
        assert get_Parent(gr[0]) is None
        gr[0][-1].append(('gene_id', gene_id))
        print_line(gr[0])

        # rest
        tid = None
        TID = None
        for line in gr[1:]:
            if line[2] != 'exon':
                # transcript of some sort
                try:
                    assert get_Parent(line) == ID
                except AssertionError:
                    log(f'ID: {ID} - Parent: {get_Parent(line)}')
                    for l in gr:
                        log(l)
                    error('')
                TID = get_ID(line)
                tid = abbr + 'T' + NID + get_ID(line)
                line[-1] += [('gene_id', gene_id), ('transcript_id', tid)]
                print_line(line)
            else:
                # exon
                if tid is None:
                    # add transcript line
                    assert get_Parent(line) == ID
                    TID = ID
                    tid = abbr + 'T' + NID + ID
                    tline = gr[0].copy()
                    tline[2] = 'transcript'
                    tline[-1] = gr[0][-1].copy() + [('gene_id', gene_id), ('transcript_id', tid)]
                    print_line(tline)

                assert get_Parent(line) == TID
                eid = abbr + 'E' + NID + get_ID(line)
                line[-1] += [('gene_id', gene_id), ('transcript_id', tid), ('exon_id', eid)]
                print_line(line)


    elif maintype in ['tRNA', 'rRNA']:

        # add gene line
        gline = gr[0].copy()
        gline[2] = 'gene'
        gline[-1] = gr[0][-1].copy() + [('gene_id', gene_id)]
        print_line(gline)

        # fix transcript line
        assert get_Parent(gr[0]) is None
        tid = abbr + 'T' + NID + ID
        gr[0][-1] += [('gene_id', gene_id), ('transcript_id', tid)]
        print_line(gr[0])

        # exon line(s)
        for line in gr[1:]:
            assert get_Parent(line) == ID
            eid = abbr + 'E' + NID + get_ID(line)
            line[-1] += [('gene_id', gene_id), ('transcript_id', tid), ('exon_id', eid)]
            print_line(line)


    elif maintype == 'pseudogene':

        # keep as is
        for line in gr:
            print_line(line)


        # # add a 'gene' line
        # gline = gr[0].copy()
        # gline[2] = 'gene'
        # gline[-1] = gline[-1].copy() + [('gene_id', gene_id)]
        # print_line(gline)

        # # add tags to rest of lines
        # tline = gr[0]
        # tid = abbr + 'T' + NID + get_ID(tline)
        # tline[-1] += [('gene_id', gene_id), ('transcript_id', tid)]
        # print_line(tline)

        # for line in gr[1:]:
        #     try:
        #         assert get_Parent(line) == ID
        #     except AssertionError:
        #         log(f'{gr}')
        #         exit()
        #     eid = abbr + 'E' + NID + get_ID(line)
        #     line[-1] += [('gene_id', gene_id), ('transcript_id', tid), ('exon_id', eid)]
        #     print_line(line)

    # all other types
    else:
        if maintype not in warn_types:
            log(f'Warning: Found unhandled main feature type: {maintype}')
            warn_types.append(maintype)

        # keep as is
        for line in gr:
            print_line(line)




#### main loop over file

linectr = 0

gene_records = {}
parents = {}

main_ftypes = {}

# iterate over file
with open(infile) as infh:
    for ln, line in enumerate(infh):


        if line.startswith('##gff-version'):
            gfftoggle = True
            log('gff format detected.')
            continue

        if line.startswith('#') or line == '\n':
            continue

        # status
        linectr += 1
        if linectr % 100 == 0:
            sys.stderr.write(f'LOG: Parsing lines: {linectr}\r')


        # normal record, split line
        try:
            contig, source, ftype, start, stop, _, strand, _, description = line.strip().split('\t')
        except:
            error(f'Line {ln} does not have 9 fields:\n{line}')


        # ignore some feature types
        if ftype in ['region', 'CDS', 'cDNA_match']:
            if ftype in del_lines:
                del_lines[ftype] += 1
            else:
                del_lines[ftype] = 1
            continue


        # read description
        if not gfftoggle and not gtftoggle:
            tag = description.split(';')[0]
            if '=' in tag:
                gfftoggle = True
                log('gff format detected.')
            elif ' "' in tag:
                gtftoggle = True
                log('gtf format detected.')
            else:
                error(f'Could not determine format from field 9: {description}')

        if gfftoggle:
            tags = [(f.split('=')) for f in description.split(';')]

        elif gtftoggle:
            # read in to check format
            tspl = description.split(';')
            if tspl[-1] == '':
                tspl = tspl[:-1]
            pairs = [t.strip().split(' ',1) for t in tspl]

            tags = [(t, v.strip('"')) for t, v in pairs]
        else:
            error(f'Could not determine format from field 9: {description}')

        # add contig
        if contig not in gene_records:
            gene_records[contig] = {}


        # line record
        lrec = [contig, source, ftype, int(start), int(stop),
                            '.', strand, '.', tags]
        ID = get_ID(lrec)
        if get_Parent(lrec) is None:

            # new gene record
            try:
                assert ID not in gene_records[contig]
            except AssertionError:
                error(f'ID {ID} already in gene records.')
            gene_records[contig][ID] = {'gline': lrec, 'trans': {}}

            # track main types
            if ftype in main_ftypes:
                main_ftypes[ftype] += 1
            else:
                main_ftypes[ftype] = 1

        else:
            par = get_Parent(lrec)
            if par in gene_records[contig]:
                # add as transcript
                gene_records[contig][par]['trans'][ID] = {'tline': lrec, 'exons': []}
                try:
                    assert ID not in parents
                except AssertionError:
                    # WTF?
                    log(f'WARNING: Duplicate ID {ID} in file!?')
                    continue
                
                # set parent
                parents[ID] = par
            
            elif par in parents:
                # add as exon
                gene = parents[par]
                gene_records[contig][gene]['trans'][par]['exons'].append(lrec)

            else:
                error(f'Line parent not found: {par}\n{lrec}')
        
        
sys.stderr.write(f'LOG: Parsing lines: {linectr}\r')
log('Parsed input file.', newline_before=True)


log(f'Main feature types: \n{sorted([(k, main_ftypes[k]) for k in main_ftypes], key=lambda x: x[1], reverse=True)}')

#####
# handle all records
log('Printing records ...')

for contig in gene_records:

    data = gene_records[contig]
    
    # sort by start pos on contig
    for ID in sorted(data.keys(), key=lambda x: data[x]['gline'][3]):

        # start with gene line
        gr = []
        gr.append(data[ID]['gline'])

        # add transcripts
        tras = data[ID]['trans']
        for tra in sorted(tras.keys(), key=lambda x: tras[x]['tline'][3]):

            gr.append(tras[tra]['tline'])
            # add exons
            for exon in sorted(tras[tra]['exons'], key=lambda x: x[3]):
                gr.append(exon)

        # print this gene record
        handle_gene_rec(gr)

log(f'Records printed: {num_lines} lines.')

log(f'Deleted {sum([del_lines[k] for k in del_lines])} lines:\n{sorted([(k, del_lines[k]) for k in del_lines], key=lambda x: x[1], reverse=True)}')

log('Done.')
