import numpy as np


class protein:
    def __init__(self, seq) -> None:
        self.seq = seq
        self.__len__ = self.__len__()

    def __len__(self):
        return len(self.seq)


class CDR(protein):
    def __init__(self, seq) -> None:
        super().__init__(seq)

    def change_one_aminoacid(self, pos, new_aminoacid):
        self.seq = self.seq[:pos] + new_aminoacid + self.seq[pos + 1:]


class FR(protein):
    def __init__(self, seq) -> None:
        super().__init__(seq)


class VHH(protein):
    def __init__(self, VHH_seq, antibody_cdr) -> None:
        super().__init__(VHH_seq)
        self.antibody_cdr = antibody_cdr
        self.FRs = {"FR1": None, "FR2": None, "FR3": None, "FR4": None}
        self.CDRs = {"CDR1": None, "CDR2": None, "CDR3": None}
        self.FRs_pos = []
        self.CDRs_pos = []
        self.partition()

    def partition(self, cdr_type='123'):
        if len(self.seq) != len(self.antibody_cdr):
            raise "length of protein and cdr not the same"

        positions = [0] + [item for sublist in [self.find_first_last(self.antibody_cdr, i) for i in cdr_type] for item in sublist] + [self.__len__]
        self.FRs_pos = [(positions[i], positions[i + 1]) for i in range(len(positions) - 1) if i % 2 == 0]
        self.CDRs_pos = [(positions[i], positions[i + 1]) for i in range(len(positions) - 1) if i % 2 == 1]

        parts = [self.seq[positions[i]:positions[i + 1]] for i in range(len(positions) - 1)]
        self.FRs = {f"FR{i + 1}": FR(parts[2 * i]) for i in range(4)}
        self.CDRs = {f"CDR{i + 1}": CDR(parts[2 * i + 1]) for i in range(3)}

    @staticmethod
    def find_first_last(cdr_str, target):
        np_arr = np.array(list(cdr_str))
        indices = np.where(np_arr == target)[0]
        if indices.size > 0:
            return indices[0], indices[-1] + 1
        else:
            return -1, -1

    def change_one_aminoacid_on_cdr(self, cdr_idx, pos, new_aminoacid):
        self.CDRs[f'CDR{cdr_idx}'].change_one_aminoacid(pos, new_aminoacid)
        self.seq = self.FRs['FR1'].seq + self.CDRs['CDR1'].seq + \
                   self.FRs['FR2'].seq + self.CDRs['CDR2'].seq + \
                   self.FRs['FR3'].seq + self.CDRs['CDR3'].seq + self.FRs['FR4'].seq

    def change_cdr(self, cdr_idx, new_cdr):
        self.CDRs[f'CDR{cdr_idx}'] = new_cdr
        self.seq = self.FRs['FR1'].seq + self.CDRs['CDR1'].seq + \
                   self.FRs['FR2'].seq + self.CDRs['CDR2'].seq + \
                   self.FRs['FR3'].seq + self.CDRs['CDR3'].seq + self.FRs['FR4'].seq


VHH_9G8 = VHH("EVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVVAINWSSGSTYYADSVKGRFTISRDNAKNTMYLQMNSLKPEDTAVYYCAAGYQINSGNYNFKDYEYDYWGQGTQVT",
              "0000000000000000000000000011111111100000000000000222222220000000000000000000000000000000000000000333333333333333333300000000")
VHH_7D12 = VHH("QVKLEESGGGSVQTGGSLRLTCAASGRTSRSYGMGWFRQAPGKEREFVSGISWRGDSTGYADSVKGRFTISRDNAKNTVDLQMNSLKPEDTAIYYCAAAAGSAWYGTLYEYDYWGQGTQVTV",
               "00000000000000000000000001111111111000000000000002222222222222222200000000000000000000000000000000333333333333333000000000")
VHH_7D12_cdr1 = VHH("QVKLEESGGGSVQTGGSLRLTCAASGRTSRSYGMGWFRQAPGKEREFVSGISWRGDSTGYADSVKGRFTISRDNAKNTVDLQMNSLKPEDTAIYYCAAAAGSAWYGTLYEYDYWGQGTQVTV",
                    "00000000000000000000000001111111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000")
VHH_7D12_shortcdr3 = VHH("QVKLEESGGGSVQTGGSLRLTCAASGRTSRSYGMGWFRQAPGKEREFVSGISWRGDSTGYADSVKGRFTISRDNAKNTVDLQMNSLKPEDTAIYYCAAAAGYDYWGQGTQVTV",
                         "00000000000000000000000001111111111000000000000002222222222222222200000000000000000000000000000000333333000000000")

VHH_EGB4 =           VHH("QVQLQESGGGSVQAGGSLKLSCAASGRSFSTYAMGWFRQAPGQDREFVATISWTDSTDYADSVKGRFTISRDNAKNTGYLQMNSLKPEDTAVYYCAADRWASSRRNVDYDYWGQGTQVTVSS",
                         "00000000000000000000000000111111111000000000000002222222000000000000000000000000000000000000000033333333333333300000000000")
VHH_EGB4_shortcdr3 = VHH("QVQLQESGGGSVQAGGSLKLSCAASGRSFSTYAMGWFRQAPGQDREFVATISWTDSTDYADSVKGRFTISRDNAKNTGYLQMNSLKPEDTAVYYCAADRWDYDYWGQGTQVTVSS",
                         "0000000000000000000000000011111111100000000000000222222200000000000000000000000000000000000000003333333300000000000")

