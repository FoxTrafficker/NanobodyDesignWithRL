import subprocess
from Bio.PDB import PDBIO, Structure
import os
from pathlib import Path
script_dir = 'src/docking/TMalign/TMalign'
PDB1_dir = "src/docking/TMalign/.tmp/PDB1.pdb"
PDB2_dir = "src/docking/TMalign/.tmp/PDB2.pdb"
matrix_dir = "src/docking/TMalign/.tmp/matrix"


def TMalign(structure1: Structure, structure2: Structure, use_fast_mode=True, output_matrix=True, use_split=True):
    prepare_file(structure1, structure2)
    cmd = [script_dir, PDB1_dir, PDB2_dir, '-outfmt', '2']

    if use_fast_mode:
        cmd.append('-fast')
    if use_split:
        cmd.extend(['-split', '2', '-ter', '1'])
    if output_matrix:
        cmd.extend(['-m', matrix_dir])

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    output, _ = process.communicate()
    print(output)
    TMscore = output.split()[13]
    for i in os.listdir(Path(matrix_dir).parent):
        if i.startswith(Path(matrix_dir).name):
            with open(Path(matrix_dir).parent / i, 'r') as f:
                matrix = [i.split()[1:] for i in f.readlines() if i[0].isdigit()]

    return TMscore, matrix


def prepare_file(structure1: Structure, structure2: Structure):
    assert isinstance(structure1, Structure.Structure) and isinstance(structure2, Structure.Structure), "输入必须是Structure类型"
    
    io = PDBIO()

    io.set_structure(structure1)
    io.save(PDB1_dir)

    io.set_structure(structure2)
    io.save(PDB2_dir)


# todo
def TMalign_batch():
    pass


# todo
def prepare_batch_file():
    pass