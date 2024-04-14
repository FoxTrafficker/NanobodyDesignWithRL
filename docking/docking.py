from TMalign.TMalign import *
from Bio.PDB import PDBParser, Structure, Chain, Model
from Bio.PDB import Superimposer
import numpy as np
import copy
import warnings
warnings.filterwarnings("ignore")


def TMDock(antibody_moving: Structure,antibody_target: Structure, antigen: Structure, output_file):
    TMscore, matrix = TMalign(antibody_moving, antibody_target, use_split=True)
    new_antibody = perform_trans(antibody_moving,matrix)

    merge_chains(new_antibody,antigen,output_file=output_file)
    return TMscore, matrix

# todo
def HDock():
    pass

def perform_trans(structure, trans_matrix, ):
    # 将输入的变换矩阵转换为NumPy数组
    trans_matrix = np.array(trans_matrix, dtype=float)
    
    # 从变换矩阵中提取平移向量和旋转矩阵
    t = trans_matrix[:, 0]
    u = trans_matrix[:, 1:]

    # 创建structure的深拷贝，以便修改而不改变原始对象
    transformed_structure = copy.deepcopy(structure)

    for model in transformed_structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # 读取原子坐标
                    coord = np.array(atom.get_coord(), dtype=float)
                    # 应用旋转和平移
                    transformed_coord = np.dot(u, coord) + t
                    # 更新原子坐标
                    atom.set_coord(transformed_coord)
    
    return transformed_structure


# todo
def merge_chains(structure1, structure2, output_file):
    # 创建一个新的Structure对象来保存合并后的结果
    merged_structure = Structure.Structure("merged")
    
    # 添加一个模型到merged_structure中，ID设置为0
    merged_model = Model.Model(0)
    
    # 遍历第一个结构的所有链，将它们的ID设置为"A"，然后添加到新模型中
    for chain in structure1.get_chains():
        chain.id = "A"
        merged_model.add(chain)
    
    # 遍历第二个结构的所有链，将它们的ID设置为"B"，然后添加到新模型中
    for chain in structure2.get_chains():
        chain.id = "B"
        merged_model.add(chain)
    
    merged_structure.add(merged_model)
    # 使用PDBIO将合并后的结构写入到一个文件中
    io = PDBIO()
    io.set_structure(merged_structure)
    io.save(output_file)



def get_chain(PDB_dir, modelid=0, chainid='A', new_chain_id='A', return_as='structure'):
    parser = PDBParser()
    chain = parser.get_structure('A', PDB_dir)[modelid][chainid]
    if return_as == 'structure':
        new_structure = Structure.Structure(new_chain_id)
        new_model = Model.Model(0)
        new_model.add(chain)
        new_structure.add(new_model)
        return new_structure

    elif return_as == 'chain':
        return chain

    else:
        raise 'return as structure or chain'


if __name__ == "__main__":
    antibody_root_dir = Path('data/1_Structures/EGFR/antibody/EGB4')
    antigen_root_dir = Path('data/1_Structures/EGFR/antigen')
    docking_res_root = Path('data/3_DockingResults/TMDock/EGFR/EGB4')

    for antibody_pdb_file in os.listdir(antibody_root_dir):
        antigen_pdb_file = '7om4.pdb'
        TMscore, matrix = TMDock(get_chain(antibody_root_dir / antibody_pdb_file, return_as='structure'), 
                                    get_chain(antigen_root_dir /antigen_pdb_file, chainid="B", return_as='structure'),
                                    get_chain(antigen_root_dir /antigen_pdb_file, chainid="A", return_as='structure'),
                                    output_file = str(docking_res_root / f'-{antibody_pdb_file[:-4]}-.pdb'))
