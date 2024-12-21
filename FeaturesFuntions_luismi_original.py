import pandas as pd # ver 2.2.2
import os
import math
import heapq
import numpy as np
import threading

wd = os.getcwd()
ann_dir = os.path.join(wd + "/annotate")
patients_dir = os.path.join(wd + "/patients")
mat_dir = os.path.join(wd + "/material")
tmp_dir = os.path.join(wd + "/tmp")

                                        
def make_input_file(seqs: list[str], soft: str):
    """Crea un archivo de entrada para NetMHC en el formato esperado."""
    
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        
    file_path = f"{tmp_dir}/{soft}_l.pep"
    

    # Escribe las secuencias de mutado y silvestre en líneas separadas
    mt = seqs[0]
    wt = seqs[1]
        
    with open(file_path, "w") as f:
        f.write(f"{mt}\n{wt}")
    return file_path

def netMHC(seqs:list[str], haplotypes:list[str], pep_features: dict, mut_features: dict):
    haplotypes_str = ",".join(haplotypes)
    input_path = make_input_file(seqs, "NetMHCpan")
    
    xls_path = f'{tmp_dir}/NetMHCpan_out.xls'  # Reemplaza esto con la ruta correcta
    os.system(f"netMHCpan -p {input_path} -a {haplotypes_str} -BA -inptype 1 -xls -xlsfile {xls_path}  > /dev/null")
    
    # Leer resultados del archivo .xls
    try:
        df = pd.read_csv(xls_path, sep='\t', header=1)
    except FileNotFoundError:
        raise FileNotFoundError(f"No se encontró la salida esperada: {xls_path}")
    
    # Extraer las afinidades del mutante y wild-type
    mut_affinity = df.loc[0, 'BA-score']  # Primera fila: mutante
    wt_affinity = df.loc[1, 'BA-score']   # Segunda fila: wild-type
    mut_rank = df.loc[0, 'BA_Rank']
    wt_rank = df.loc[1, 'BA_Rank']

    # Calcular la amplitud A
    A = (wt_affinity / mut_affinity) if mut_affinity > 0 else 0

    # Actualizar el diccionario de características
    pep_features.update({
        "mutant_BA_score_netMHCpan": mut_affinity,
        "wildtype_BA_score_netMHCpan": wt_affinity,
        "mutant_BA_Rank_netMHCpan": mut_rank,
        "wildtype_BA_Rank_netMHCpan": wt_rank,
        "amplitude_A": A
    })

    # Salida
    print(f"Mutante BA_score: {mut_affinity}, Wild-type BA_score: {wt_affinity}, Amplitud A: {A}")
    print(f"Rank MT: {mut_rank}, Rank wt: {wt_rank}")
    return mut_affinity, wt_affinity, A, mut_rank, wt_rank



def iniciate_dict(sum: dict, pep: bool = False):
    """Inicializa un diccionario de características basado en la mutación."""
    required_keys = ["patient", "nc_mut", "aa_mut", "gen", "consequence", "ENST", "wt_seq", "mut_seq", "mut_pos"]
    for key in required_keys:
        if key not in sum:
            raise KeyError(f"Clave faltante en la entrada: {key}")
    
    features = {
        "patient": sum["patient"],
        "chr": sum["nc_mut"].split(".")[0],
        "genomic_coord": sum["nc_mut"].split(".")[1],
        "protein_coord": sum["aa_mut"].split(".")[0],
        "ref": sum["nc_mut"].split(".")[-1].split("/")[0],
        "alt": sum["nc_mut"].split(".")[-1].split("/")[-1],
        "gen": sum["gen"],
        "consequence": sum["consequence"],
        "transcript": sum["ENST"],
        "wt_seq": sum["wt_seq"],
        "mut_seq": sum["mut_seq"],
        "mut_pos": sum["mut_pos"],
    }
    if pep:
        features["seq_len"] = sum["len"]

    return features

def transfer_features(pep_features: dict, mut_features: dict, mode: str):
    if mode == "pep2mut":
        # Transferir características relacionadas con la afinidad
        if "mutant_score_netMHCpan" in pep_features and "wildtype_score_netMHCpan" in pep_features:
            mut_features["mutant_score_netMHCpan"] = pep_features["mutant_score_netMHCpan"]
            mut_features["wildtype_score_netMHCpan"] = pep_features["wildtype_score_netMHCpan"]
        
        if "mutant_BA_Rank_netMHCpan" in pep_features and "wildtype_BA_Rank_netMHCpan" in pep_features:
            mut_features["mutant_BA_Rank_netMHCpan"] = pep_features["mutant_BA_Rank_netMHCpan"]
            mut_features["wildtype_BA_Rank_netMHCpan"] = pep_features["wildtype_BA_Rank_netMHCpan"]
        
        if "amplitude_A" in pep_features:
            mut_features["amplitude_A"] = pep_features["amplitude_A"]

    elif mode == "mut2pep":
        # Transferir características de la mutación al péptido
        if "mutant_score_netMHCpan" in mut_features and "wildtype_score_netMHCpan" in mut_features:
            pep_features["mutant_score_netMHCpan"] = mut_features["mutant_score_netMHCpan"]
            pep_features["wildtype_score_netMHCpan"] = mut_features["wildtype_score_netMHCpan"]
        
        if "mutant_BA_Rank_netMHCpan" in mut_features and "wildtype_BA_Rank_netMHCpan" in mut_features:
            pep_features["mutant_BA_Rank_netMHCpan"] = mut_features["mutant_BA_Rank_netMHCpan"]
            pep_features["wildtype_BA_Rank_netMHCpan"] = mut_features["wildtype_BA_Rank_netMHCpan"]
        
        if "amplitude_A" in mut_features:
            pep_features["amplitude_A"] = mut_features["amplitude_A"]


def process_peptide(seqs, mut_neo_ext, mut_neo_pos, lengt, haplotypes, pep_row, mut_row_features, min_values):
    pep_row_features = iniciate_dict(pep_row, pep=True)
    
    # Verificar que haplotypes tenga al menos un elemento
    if len(haplotypes) < 1:
        raise ValueError("No se proporcionaron haplotipos válidos para la predicción.")

    # Llamada a netMHC para obtener el valor de afinidad
    # netMHC debe devolver un valor, necesitamos actualizar este valor después de join
    mut_affinity, wt_affinity, A, mut_rank, wt_rank = netMHC(seqs, haplotypes[1], pep_row_features, mut_row_features)

    # Asegurarse de que pep_row_features esté actualizado con los valores de afinidad y amplitud
    pep_row_features["mutant_score_netMHCpan"] = mut_affinity
    pep_row_features["wildtype_score_netMHCpan"] = wt_affinity
    pep_row_features["mutant_BA_Rank_netMHCpan"] = mut_rank  # mut_rank debe ser devuelto por netMHC
    pep_row_features["wildtype_BA_Rank_netMHCpan"] = wt_rank  # wt_rank debe ser devuelto por netMHC
    pep_row_features["amplitude_A"] = A

    # Guardar los valores para los features derivados y mantener los 3 más pequeños
    for key in min_values.keys():
# Verifica si la clave existe en pep_row_features
        if key in pep_row_features:
            if key == "mutant_BA_Rank_netMHCpan":
                heapq.heappush(min_values[key], [pep_row_features[key], wt_affinity])
            else:
                heapq.heappush(min_values[key], pep_row_features[key])
            # Limitar la lista a los 3 valores más pequeños
            if len(min_values[key]) > 3:
                min_values[key] = min_values[key][:3]
        else:
            print(f"Warning: Key '{key}' not found in pep_row_features")

    return pep_row_features, min_values