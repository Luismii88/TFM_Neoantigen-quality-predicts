def Make_MissInser(wt_protein: str, aa_position: int, alt: str, flank_length: int= 12):
    """
    Función para obtener una secuencia de 25 aa, modelo para snvs
    """
    total_len = 1 + flank_length*2
    
    aa_position -= 1

        # Calcular los límites del mer
    start = max(aa_position - flank_length, 0)
    end = min(aa_position + flank_length + 1, len(wt_protein))

    # Construir la secuencia mutada (alt es un solo aminoácido)
    mutated_mer = wt_protein[start:aa_position] + alt + wt_protein[aa_position + 1:end]
    mutated_mer_ext = mutated_mer + wt_protein[end:end + 10]  # Extensión de 10 aa adicionales a la derecha
    wt_mer = wt_protein[start:end]
    new_mutation_position = aa_position - start + 1  # Posición relativa de la mutación en la secuencia final

    return mutated_mer, mutated_mer_ext, wt_mer, new_mutation_position



def NeoGenerator(mut_mer: str, wt_mer: str, mut_pos: int, mer_ext: str):
    mut_pos -= 1
    
    for length in range(8, 13):
        l = max(mut_pos - length + 1, 0)
        r = mut_pos + 1 if l > 0 else length
        
        while l <= mut_pos and r <= len(mut_mer):
            mut_neo = mut_mer[l:r]
            mut_neo_ext = mer_ext[l:r+11] if r + 11 <= len(mer_ext) else mer_ext[l:]
            wt_neo = wt_mer[l:r]

            mut_neo_pos = int(mut_pos - l + 1)

            yield mut_neo, mut_neo_ext, wt_neo, mut_neo_pos, length

            l += 1
            r += 1