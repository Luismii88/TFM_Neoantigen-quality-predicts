import pandas as pd # ver 2.2.2
import os
from pysam import VariantFile # ver 0.22.1
import datetime
import time
import argparse
from CheckClass import Check_Point
from SecuenceFuntions_luismi import *
from FeaturesFuntions_luismi_original import *
import numpy as np

parser = argparse.ArgumentParser(
    prog="Feature_extractor"
)
parser.add_argument("-f", "--force", action="store_true", help="overwrite ouput", default=False)
args = parser.parse_args()

### mut dict
mut_dict = {
    "patient": str,
    "chromosome": int,
    "genomic_coord": float,
    "ref": str,
    "alt": str,
    "gene": str,
    "transcript": str,
    "protein_coord": int,
    "aa_mutant": str,
    "aa_wt": str,
    "mutant_seq": str,
    "wt_seq": str,
    "consequence": str,
    "mutant_BA_score_netMHCpan": float,
    "wildtype_BA_score_netMHCpan": float,
    "mutant_BA_Rank_netMHCpan": float,
    "wildtype_BA_Rank_netMHCpan": float,
    "amplitude_A": float
}

### pep dict
pep_dict = {
    "patient": str,
    "chromosome": int,
    "genomic_coord": float,
    "ref": str,
    "alt": str,
    "gene": str,
    "transcript": str,
    "protein_coord": int,
    "aa_mutant": str,
    "aa_wt": str,
    "mutant_seq": str,
    "wt_seq": str,
    "consequence": str,
    "mutant_BA_score_netMHCpan": float,
    "wildtype_BA_score_netMHCpan": float,
    "mutant_BA_Rank_netMHCpan": float,
    "wildtype_BA_Rank_netMHCpan": float,
    "amplitude_A": float
}

wd = os.getcwd()
ann_dir = os.path.join(wd + "/annotate")
patients_dir = os.path.join(wd + "/patients")
mat_dir = os.path.join(wd + "/material")


with open(f"{mat_dir}/MHC_pseudo.txt", "r") as f:
        # Leer todas las líneas y quitar espacios en blanco y saltos de línea
        aviables_hlas = [line.strip() for line in f.readlines()]


def get_HLA(patient: str):
    file = f"{patients_dir}/{patient}/{patient}_HLA.txt"
    Mix_haplo = []
    net_haplo = []
    with open(file, "r") as f:
        l = f.readline().split(",")
        for hla in l:
            if hla == "-" or hla in ["Not","typed"]:
                continue
            else:
                parts = hla.split("*")
                hla = parts[0]
                f_digits = parts[1].split(":")[0]
                s_digits = parts[1].split(":")[1]

                net = hla + f_digits + ":" + s_digits

                if net in aviables_hlas:
                    Mix_haplo.append(hla[-1] + f_digits + s_digits) #A0101
                    net_haplo.append(net) #HLA-A01:01
                 
        return None if net_haplo == [] else (Mix_haplo, net_haplo) # si no se logra optener un HLA, pues no y pasamos del paciente



def main(vcf_file: str, patient: str, haplotypes: list[str], log_path: str):
    targets = ["missense_variant"]

    all_mut_sum = pd.DataFrame(columns=["patient", "ENST", "gen", "consequence", "nc_mut", "aa_mut", "wt_seq", "mut_seq", "mut_pos"])
    all_peps_sum = pd.DataFrame(columns=["patient", "ENST", "gen", "consequence", "nc_mut", "aa_mut", "wt_seq", "mut_seq", "mut_pos", "len"])

    all_mut_features = pd.DataFrame(columns=mut_dict)
    all_pep_features = pd.DataFrame(columns=pep_dict)  # Asegúrate de que empieza como un DataFrame
    
    check_point = Check_Point()

    # Abrir el archivo VCF
    vcf = VariantFile(vcf_file)
    n_mut = 0

    for header_line in vcf.header.records:
            header_line = str(header_line)
            if '##tumor_sample=' in header_line:
                tumor = header_line.split("=")[1][:-1]
                break
            
    for record in vcf.fetch():
        # La variante pertenece solo a la muestra de cancer
        if record.samples[tumor]["GT"] in [(0, 1), (1, 1), (1, 2)]:
            if 'CSQ' in record.info:
                # Las anotaciones de CSQ están separadas por comas
                for i in range(len(record.info['CSQ'])):
                    annotations = record.info['CSQ'][i].split(',')
                    for annotation in annotations:
                        fields = annotation.split('|')
                        # El índice 1 corresponde a "consecuencia"
                        type = fields[1] 
                        if type in targets:  # Solo missense_variant
                                                    
                            aa_pos = int(fields[14].split("-")[0])  # Convertir posición a entero
                            gen = fields[3]
                            wt_prot = fields[28]
                            
                            # Validación de posición con check_1
                            if check_point.check_1(aa_pos=str(aa_pos), wt_prot=wt_prot):
                                break
                            
                            ref = fields[15].split("/")[0]
                            alt = fields[15].split("/")[-1]

                            # Generación de secuencia mutada para missense_variant
                            mutated_mer, mutated_mer_ext, wt_mer, mutation_position = Make_MissInser(wt_prot, aa_pos, alt)
                                                    
                            chr = record.chrom
                            genomic_coord = record.pos
                            nc_ref = record.ref
                            nc_alt = record.alts[0]  # assuming single alt allele
                            
                            nc_tag = f"{chr}.{genomic_coord}.{nc_ref}/{nc_alt}"
                            aa_tag = f"{aa_pos}.{ref}/{alt}"
                            enst = fields[6]

                            # Checkpoints 2 y 3
                            check_2, check_3 = check_point.check_2(mut_seq=mutated_mer, wt_seq=wt_mer), check_point.check_3(mut_seq=mutated_mer, sum=all_mut_sum, gen=gen, enst=enst, aa_tag=aa_tag)
                            if check_2 or check_3:
                                break

                            n_mut += 1
                            print(f"procesando {mutated_mer}\nnc_tag {nc_tag}\naa_tag {aa_tag}")

                            mut_row = {"patient": patient, "ENST": enst, "gen": gen, "consequence": type, 
                                    "nc_mut": nc_tag, "aa_mut": aa_tag, "wt_seq": wt_mer, "mut_seq": mutated_mer,  
                                    "mut_pos": mutation_position}
                            
                            mut_row_features = iniciate_dict(mut_row)
                            
                            # Generación de neopeptidos y procesamiento con NetMHC
                            neo_gen = NeoGenerator(mutated_mer, wt_mer, mutation_position, mutated_mer_ext)
                            min_values = {"mutant_score_netMHCpan": []}
                            
                            try:
                                while True:
                                    mut_neo, mut_neo_ext, wt_neo, mut_neo_pos, lengt = next(neo_gen)
                                    # Verificación de secuencia mutada y de tipo salvaje
                                    if check_point.check_2(mut_seq=mut_neo, wt_seq=wt_neo):
                                        continue
                                    if check_point.check_3(mut_seq=mut_neo, sum=all_peps_sum, gen=gen, enst=enst, aa_tag=aa_tag):
                                        continue

                                    print(mut_neo)

                                    # Asegúrate de que all_pep_features sea un DataFrame antes de concatenar
                                    if isinstance(all_pep_features, list):
                                        all_pep_features = pd.DataFrame(all_pep_features)

                                    # Definir pep_row dentro del ciclo
                                    pep_row = {"patient": patient, "ENST": enst, "gen": gen, "consequence": type, 
                                            "nc_mut": nc_tag, "aa_mut": aa_tag, "wt_seq": wt_neo, "mut_seq": mut_neo, 
                                            "mut_pos": mut_neo_pos, "len": lengt}
                                    
                                    # Convertir pep_row a un DataFrame
                                    pep_row_df = pd.DataFrame([pep_row])

                                    # **NO** eliminar columnas vacías
                                    # Concatenar pep_row_df al DataFrame all_peps_sum sin usar dropna
                                    all_peps_sum = pd.concat([all_peps_sum, pep_row_df], ignore_index=True)

                                    seqs = [mut_neo, wt_neo]
                                    pep_row_features, min_values = process_peptide(seqs, mut_neo_ext, mut_neo_pos, lengt, haplotypes, pep_row, mut_row_features, min_values)
                                    all_pep_features = pd.concat([all_pep_features, pd.DataFrame([pep_row_features])], ignore_index=True)

                            except StopIteration:
                                print(f"neopeptidos de {enst} generados")
                        
    vcf.close()
    print(f"{n_mut} mutaciones procesadas")
    return all_mut_sum, all_peps_sum, all_mut_features, all_pep_features



if __name__ == "__main__":
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H:%M:%S")
    log_path = f"{wd}/logs/failed_muts_log_{timestamp}.txt"
    start_time = time.time()
    
    # Asegúrate de que los directorios existen
    log_dir = f"{wd}/logs"
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)

    for patient in os.listdir(patients_dir):
        print(f"Procesando paciente: {patient}")
        
        patient_dir = os.path.join(patients_dir, patient)
        print(f"Directorio del paciente: {patient_dir}")
        if not os.path.exists(patient_dir):
            os.makedirs(patient_dir)  # Crear el directorio si no existe

        if not args.force and os.path.isfile(f"{patient_dir}/all_mut_sum.csv"):
            print(f"{patient} ya ha sido procesado")
            continue

        vcf_file = os.path.join(ann_dir, f"{patient}_hc_gatk-mutect2_ann_fil.vcf.gz")

        haplotypes = get_HLA(patient)
        if haplotypes:
            all_mut_sum, all_pep_sum, all_mut_features, all_pep_features = main(vcf_file, patient, haplotypes, log_path)
            
            # Verifica si los DataFrames tienen datos antes de intentar guardarlos
            print("Guardando resultados para:", patient)
            print("all_mut_sum:", all_mut_sum.head())
            print("all_pep_sum:", all_pep_sum.head())

            all_mut_sum.to_csv(f"{patient_dir}/all_mut_sum.csv", index=False)
            all_pep_sum.to_csv(f"{patient_dir}/all_pep_sum.csv", index=False)
            all_mut_features.to_csv(f"{patient_dir}/all_mut_features.txt", sep="\t", index=False)
            all_pep_features.to_csv(f"{patient_dir}/all_pep_features.txt", sep="\t", index=False)

        else:
            print(f"{patient} no tiene haplotipos reconocibles")
            with open(f"{wd}/logs/failed_patient_log_{timestamp}.txt", "a") as f:
                f.write(f"{patient}\n")
            
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Tiempo transcurrido: {elapsed_time} segundos")
else:
    print(f"{__name__}, no estas en main")


