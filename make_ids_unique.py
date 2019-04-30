import sys
import os

anno_gff = '/data/prostlocal2/projects/mh_bats_ncrna_annotation/2018/annotations/abbr/efu.gff'
anno_gff_new = '/home/ji57pog/Projects/mh_bats_ncrna_annotation/util/efu.gff'

root_dict = {}
root_id_list = []
child_dict = {}

with open(anno_gff, 'r') as inp:
    for line in inp:
        if not line.startswith('#'):

            # normal record, split line
            try:
                contig, source, ftype, start, stop, _, strand, _, description = line.strip().split('\t')
            except:
                print(f'Line {line} does not have 9 fields:\n{line}')
                sys.exit(1)

            tags = [(f.split('=')) for f in description.split(';')]

            # get id
            ID = [v for t, v in tags if t == 'ID']
            if ID == []:
                print(f'No ID in line: {line}')
                sys.exit(1)
            else:
                ID = ID[0]

            # get parent
            parent = [v for t, v in tags if t == 'Parent']
            if parent == []:
                # no parnet -> top level
                root_id_list.append(ID)
                if ID not in root_dict:
                    root_dict[ID] = [line]
                else:
                    # print(f'Duplicated id {ID} at root level:\n{line}')
                    root_dict[ID].append(line)
            else:

                parent = parent[0]
                    
                if parent in child_dict:
                    child_dict[parent].append((ID, line))
                else:
                    child_dict[parent] = [(ID, line)]

def change_id_in_record(old_ID, new_ID, record):
    first_part = record.split(old_ID, 1)[0]
    second_part = record.split(old_ID, 1)[1]

    return old_ID + '_' + str(new_ID), first_part + old_ID + '_' + str(new_ID) + second_part

def change_parent_id_in_record(old_parent_ID, new_parent_ID, record):
    first_part = record.split(old_parent_ID, 1)[0]
    second_part = record.split(old_parent_ID, 1)[1]

    return first_part + str(new_parent_ID) + second_part

def traverse(ID, start, end, out_file, new_parent_id=None):
    global counter
    if ID in child_dict:
        # ID has one or more children

        next_rec = child_dict[ID]
        next_ids = [i for i, l in next_rec]

        # print('-------')
        # print(ID, next_ids, start, end)

        if len(set(next_ids)) == len(next_ids):
            # unique child ids

            for child in next_rec:
                start_child = int(child[1].split('\t')[3])
                end_child = int(child[1].split('\t')[4])

                if start <= start_child and end >= end_child:
                    if new_parent_id:
                        out_file.write(change_parent_id_in_record(ID, new_parent_id, child[1]))
                    else:
                        out_file.write(child[1])
            for child in next_rec:
                start_child = int(child[1].split('\t')[3])
                end_child = int(child[1].split('\t')[4])

                if start <= start_child and end >= end_child:
                    traverse(child[0], int(child[1].split('\t')[3]), int(child[1].split('\t')[4]), out_file)

        else:
            # one or more ambiguous child ids
            for child in next_rec:
                start_child = int(child[1].split('\t')[3])
                end_child = int(child[1].split('\t')[4])

                if start <= start_child and end >= end_child:

                    if next_ids.count(child[0]) == 1:
                        # unique
                        if new_parent_id:
                            out_file.write(change_parent_id_in_record(ID, new_parent_id, child[1]))
                        else:
                            out_file.write(child[1])
                        traverse(child[0], int(child[1].split('\t')[3]), int(child[1].split('\t')[4]), out_file)
                    else:
                        # ambiguous
                        new_ID, new_rec = change_id_in_record(child[0], counter, child[1])
                        
                        if new_parent_id:
                            out_file.write(change_parent_id_in_record(ID, new_parent_id, new_rec))
                        else:
                            out_file.write(new_rec)

                        counter += 1

                        traverse(child[0], int(child[1].split('\t')[3]), int(child[1].split('\t')[4]), out_file, new_ID)

counter = 1

if os.path.isfile(anno_gff_new):
    os.remove(anno_gff_new)

with open(anno_gff_new, 'a') as outf:

    for r in root_id_list:
        root_rec = root_dict[r]

        if len(root_rec) == 1:
            # unique id
            outf.write(root_rec[0])

            traverse(r, int(root_rec[0].split('\t')[3]), int(root_rec[0].split('\t')[4]), outf)
        
        else:
            # ambiguous id
            outf.write(change_id_in_record(r, counter, root_rec[0])[1])
            counter += 1

            if r in child_dict:
                # TODO get children of ambiguous roots with new parent ids
                print('Skipping children of ambiguous roots: ' + root_rec[0].strip())