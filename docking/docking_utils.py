from Bio.PDB import PDBParser, PDBIO, Select
import torch


def kabsch(A, B):
    a_mean = A.mean(dim=1, keepdims=True)
    b_mean = B.mean(dim=1, keepdims=True)
    A_c = A - a_mean
    B_c = B - b_mean
    # Covariance matrix
    H = torch.bmm(A_c.transpose(1,2), B_c)  # [B, 3, 3]
    U, S, V = torch.svd(H)
    # Flip
    sign = (torch.det(U) * torch.det(V) < 0.0)
    if sign.any():
        S[sign] = S[sign] * (-1)
        U[sign,:] = U[sign,:] * (-1)
    # Rotation matrix
    R = torch.bmm(V, U.transpose(1,2))  # [B, 3, 3]
    # Translation vector
    t = b_mean - torch.bmm(R, a_mean.transpose(1,2)).transpose(1,2)
    A_aligned = torch.bmm(R, A.transpose(1,2)).transpose(1,2) + t
    return A_aligned, R, t
def prepare_file_for_TMalign(pdb1,chainid1,pdb2,chainid2,):
    # 使用 PDBParser 读取原始 PDB 文件
    structure = parser.get_structure("A", "/home/sding/Project/Nanobody_Design_with_RL/data/1_Structures/EGFR/antigen/4krl.pdb")

    # 初始化 PDBIO 实例
    io = PDBIO()

    # 设置要处理的结构
    io.set_structure(structure)

    # 使用 ChainSelect 选择器来选择特定的链，比如链 'A'
    io.save("output_chain_file.pdb", ChainSelect("A"))

def 