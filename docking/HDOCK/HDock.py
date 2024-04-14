import subprocess
from pathlib import Path
import os
import shutil

def HDockDocking(antibody_dir, antigen_dir, script_path='src/docking/HDOCK/hdock', angle='15', restr='restr.txt', out='Hdock.out'):
    cmd = f"{script_path} {antigen_dir} {antibody_dir} -angle {angle} -restr {restr} -out {out}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

def HDockCreatepl(Hdock_res, script_path='src/docking/HDOCK/createpl', nmax='100', restr='restr.txt'):
    output = f'{Hdock_res}-top{nmax}.pdb'
    cmd=f'{script_path} {Hdock_res} {output} -nmax {nmax} -restr {restr}  -complex -chid -models'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

def batchHDockDocking(antibody_dirs, antigen_dirs, HDock_res_dir, models_dir, angle = 15, nmax=100):
    for antibody_dir in antibody_dirs:
        for antigen_dir in antigen_dirs:

            HDock_name = f'Hdock-{antibody_dir.stem}-{antigen_dir.stem}'
            HDockDocking(antibody_dir, antigen_dir, out=HDock_name,angle=angle)
            HDockCreatepl(HDock_name, nmax=nmax)
            
            # 移动结果
            os.makedirs(HDock_res_dir, exist_ok=True)
            shutil.move(HDock_name,HDock_res_dir)
            for i in os.listdir('./'):
                if i.startswith('model'):
                    os.makedirs(Path(models_dir) / HDock_name, exist_ok=True)
                    shutil.move(i,Path(models_dir) / HDock_name)

def batchHDockCreatepl(HDock_names, models_dir,nmax=100):
    for HDock_name in HDock_names:
        HDockCreatepl(HDock_name, nmax=nmax)
        for i in os.listdir('./'):
                if i.startswith('model'):
                    os.makedirs(Path(models_dir) / HDock_name, exist_ok=True)
                    shutil.move(i,Path(models_dir) / HDock_name)


if __name__ == '__main__':
    # batchHDockCreatepl(['Hdock-VHH1832-3ZSJ','Hdock-VHH1834-3ZSJ','Hdock-VHH1835-3ZSJ','Hdock-VHH1836-3ZSJ','Hdock-VHH1837-3ZSJ','Hdock-WT-3ZSJ'],
    #                    models_dir='data/3_DockingResults/HDock/models')
    
    # raise
    root_dir=Path('data/1_Structures/galectin-3')

    antibody_folder = root_dir / 'antibody/VHH'
    antigen_folder= root_dir / 'antigen'
    
    antibody_dirs = [antibody_folder / file_name for file_name in os.listdir(antibody_folder)]
    antigen_dirs = [antigen_folder / file_name for file_name in os.listdir(antigen_folder)]

    batchHDockDocking(antibody_dirs, antigen_dirs, angle=15,nmax=5000,
                      HDock_res_dir='data/3_DockingResults/HDock/docking_result',
                      models_dir='data/3_DockingResults/HDock/models')

