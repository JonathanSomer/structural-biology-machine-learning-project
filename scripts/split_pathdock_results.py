# compatibility for python 2 and 3
from builtins import str

import os
import traceback

from Constants import BASE_DATA_PATH

overwrite_existing = False


def parse_pdb_line(pdb_line):
    pdb_line_tup = pdb_line.split()
    # ATOM NUM parse
    line_atom = int(pdb_line_tup[1])
    # CHAIN parse
    line_chain = pdb_line_tup[4] if len(pdb_line_tup[4]) == 1 and pdb_line_tup[4].isalpha() else None
    # RESIDUE parse
    if line_chain is None:
        maybe_chain = pdb_line_tup[4][0]
        maybe_res = pdb_line_tup[4][1:]
        if maybe_chain.isalpha() and str(maybe_res).isdecimal() and int(maybe_res) >= 1000:
            line_res = int(maybe_res)
        else:
            line_res = int(pdb_line_tup[4])
    else:
        # sometimes there is an 'A' in the residue id for some dumb reason
        if pdb_line_tup[5][-1].isalpha():
            line_res = int(pdb_line_tup[5][:-1])
        else:
            line_res = int(pdb_line_tup[5])
    return line_atom, line_chain, line_res


for root, dirs, files in os.walk(BASE_DATA_PATH):
    for f in files:
        filename, file_extension = os.path.splitext(f)
        fullpath = os.path.join(root, f)
        direct_dir = os.path.basename(os.path.dirname(fullpath))
        if direct_dir != 'patch_dock' or file_extension != '.pdb':
            continue
        new_dirpath = os.path.join(root, 'processed')
        if not os.path.exists(new_dirpath):
            print("Creating dir {}".format(new_dirpath))
            os.makedirs(new_dirpath)

        receptor_file = os.path.join(new_dirpath, filename + '.r' + file_extension)
        ligand_file = os.path.join(new_dirpath, filename + '.l' + file_extension)
        # skip if exist and shouldn't overwrite
        if not overwrite_existing and os.path.exists(receptor_file) and os.path.exists(ligand_file):
            continue

        receptor_lines = list()
        ligand_lines = list()
        is_receptor = True
        prev = None
        cur = None
        try:
            with open(fullpath) as buf:
                for line in buf:
                    cur = line
                    # check if we got to the ligand
                    if is_receptor and prev is not None:
                        cur_atom, cur_chain, cur_res = parse_pdb_line(line)
                        prev_atom, prev_chain, prev_res = parse_pdb_line(prev)
                        # ATOM and RESIDUE should be monotonically increasing
                        if cur_atom < prev_atom or cur_res < prev_res:
                            is_receptor = False
                    prev = line
                    # append line to receptor or ligand
                    if is_receptor:
                        receptor_lines.append(line)
                    else:
                        ligand_lines.append(line)
        except Exception as e:
            print ("Error parsing pdb file {}:{}".format(f, str(e)))

            continue

        with open(receptor_file, 'w') as rec_buf:
            rec_buf.writelines(receptor_lines)

        with open(ligand_file, 'w') as lig_buf:
            lig_buf.writelines(ligand_lines)
