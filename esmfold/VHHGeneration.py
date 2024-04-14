from Protein import CDR, VHH_9G8, VHH_7D12,VHH_7D12_shortcdr3
from Protein import VHH_EGB4,VHH_EGB4_shortcdr3
import random
import subprocess


def randomly_generate_cdr(cdr_length):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(random.choice(amino_acids) for _ in range(cdr_length))


def assemble_random_CDR_to_VHH(VHH):
    for i in range(3):
        new_cdr = CDR(randomly_generate_cdr(cdr_length=VHH.CDRs[f'CDR{i + 1}'].__len__))
        VHH.change_cdr(cdr_idx=i + 1, new_cdr=new_cdr)
    return VHH.seq


def generate_dataset(dataset_size, VHH, filename):
    num_digits = 6
    sequences = set()
    originalVHHseq = VHH.seq
    while len(sequences) < dataset_size:
        sequence = assemble_random_CDR_to_VHH(VHH)
        sequences.add(sequence)

    sequences = list(sequences)

    with open(filename, 'w') as fasta_file:
        header = '>originalVHH\n'
        fasta_file.write(header)
        fasta_file.write(originalVHHseq + '\n')
        for i, sequence in enumerate(sequences):
            formatted_index = str(i).zfill(num_digits)
            header = f'>{formatted_index}\n'
            fasta_file.write(header)
            fasta_file.write(sequence + '\n')


def esmfold(fasta_file_dir, output_dir):
    conda_env_name = 'esmfold'
    command = 'esm-fold'

    cmd = ['conda', 'run', '-n', conda_env_name, command, '-i', fasta_file_dir, '-o', output_dir]
    subprocess.run(cmd, text=True, capture_output=True)


if __name__ == '__main__':
    fasta_file_dir = 'src/esmfold/random_cdr_VHH_EGB4.fasta'
    structures_dir = 'data/1_Structures/EGFR/antibody/EGB4'
    generate_dataset(200000, VHH=VHH_EGB4_shortcdr3, filename=fasta_file_dir)
    esmfold(fasta_file_dir,structures_dir)
