import pandas as pd

class Check_Point():
    """This class if for creating a heckpoint object wich methots check that the secuences have certaing requiraments, if they pass the check
    they return Flase as tere no problem, if dont they return True"""
    def check_1(self, aa_pos: int, wt_prot: str):
        """1º checkpoint, makes sure that the anotation dosent put the mutation in a aminoacid beyond the protein length"""
        if int(aa_pos.split("-")[-1]) > len(wt_prot):
            return True
        return False
        
    def check_2(self, mut_seq: str, wt_seq: str):
        """2º checkpoint, makes sure that the secuence contains only known aminoacids and the secuences arent the same"""
        AAs = "ACDEFGHIKLMNPQRSTVWY"
        
        if mut_seq != wt_seq:
            for i in range(len(mut_seq)):
                if mut_seq[i] not in AAs or wt_seq[i] not in AAs:
                    return True
        else:
            return True
        return False
    
    def check_3(self, mut_seq: str, sum: pd.DataFrame, enst: str, aa_tag: str, gen: str):
        """3º checkpoint, using the summary df as cache, checks that the sequence (mut or neopep) didn't appear before. 
        It starts at the end of the df since repetitions are more likely there."""
        repetition = False
        for idx in reversed(sum.index):
            if sum.at[idx, "mut_seq"] == mut_seq and sum.at[idx, "gen"] == gen:
                sum.at[idx, "ENST"] += "/ " + enst
                if len(mut_seq) == 25:  
                    sum.at[idx, "aa_mut"] += " " + aa_tag
                repetition = True
                break
            elif sum.at[idx, "gen"] != gen:
                break

                    
        return repetition
