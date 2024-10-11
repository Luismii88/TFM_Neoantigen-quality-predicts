#!/bin/bash
# Define el intérprete de bash para ejecutar el script.

start=`date +%H:%M:%S`
# Guarda la hora de inicio del script para calcular el tiempo de ejecución.

work_dir=$(pwd)
# Define la variable "work_dir" que será el directorio de trabajo actual.

p_dir=$work_dir/vcf
# Define la variable "p_dir" para almacenar la ruta de los archivos VCF de entrada.

p2_dir=$work_dir/vcf_hc
# Define la variable "p2_dir" para almacenar los archivos VCF procesados.

m_dir=$work_dir/material
# Define la variable "m_dir" para almacenar archivos de materiales de referencia (como el genoma).

a_dir=$work_dir/annotate
# Define la variable "a_dir" para almacenar los archivos VCF anotados.

raw_files=$(ls $p_dir | grep -v 'tbi' | grep 'indel')
# Encuentra los archivos VCF de inserciones/eliminaciones (indel) y excluye los archivos .tbi.

for file in $raw_files
do
    patient_id=$(echo "$file" | cut -d'.' -f1,2 | tr '.' '_')
    # Extrae el ID del paciente del nombre del archivo y reemplaza los puntos con guiones bajos.

    method=$(echo "$file" | cut -d '.' -f 6)
    # Extrae el método (por ejemplo, mutect2 o sanger) del nombre del archivo.

    mut=$(echo "$file" | cut -d '.' -f 8)
    # Extrae el tipo de mutación (indel o snv) del nombre del archivo.

    echo "---------------------------------------------------------------------------------"
    output_vcf="${patient_id}_${method}"
    # Define el nombre del archivo VCF de salida usando el ID del paciente y el método.

    echo "Creating: ${output_vcf}.vcf.gz"
    # Muestra un mensaje indicando que se está creando el archivo de salida.

    snv=$(echo $file | sed 's/indel/snv/')
    # Cambia "indel" por "snv" en el nombre del archivo para encontrar el archivo de variantes de un solo nucleótido (SNV).

    bcftools concat "$p_dir/$file" "$p_dir/$snv" -a -D --threads 9 -Oz -o "$p2_dir/${output_vcf}.vcf.gz"
    # Usa bcftools para concatenar los archivos VCF de indel y snv, comprime y guarda el resultado en el directorio de salida.

    echo "VCF files concatenated"
    # Muestra un mensaje indicando que los archivos VCF han sido concatenados.

    if [[ $file == *mutect2* ]]
    then
        echo "File from ${method}... Correcting"
        # Si el archivo es de mutect2, se corrige el formato.

        zcat $p2_dir/${output_vcf}.vcf.gz | python3 $work_dir/bin/correct_mutect.py | bgzip > $p2_dir/${output_vcf}.vcf && tabix -p vcf $p2_dir/${output_vcf}.vcf
        # Descomprime el archivo, lo pasa por un script de corrección y lo vuelve a comprimir.

        echo "Recompressing"
        bcftools view -I "$p2_dir/${output_vcf}.vcf" -O z -o "$p2_dir/${output_vcf}.vcf.gz"
        # Vuelve a comprimir el archivo VCF corregido.

        rm "$p2_dir/${output_vcf}.vcf"
        # Elimina el archivo temporal descomprimido.

    fi

    echo "Indexing"
    bcftools index -f -t -o "${p2_dir}/${output_vcf}.vcf.tbi" "${p2_dir}/${output_vcf}.vcf.gz"
    # Indexa el archivo VCF comprimido para su uso futuro.

    if [[ ${output_vcf}.vcf.gz == *sanger* ]]
    then
        echo "Filtering common variants separately for mutect2 and sanger"
        # Filtra las variantes comunes entre los métodos Mutect2 y Sanger, pero las procesa por separado.

        bcftools isec -p $work_dir/tmp/isec $p2_dir/${patient_id}_gatk-mutect2.vcf.gz $p2_dir/${patient_id}_sanger-wgs.vcf.gz
        # Usa bcftools isec para encontrar las variantes comunes y separarlas.

        echo "Saving mutect2 filtered variants"
        bcftools view -I "$work_dir/tmp/isec/0000.vcf" -O z -o "$p2_dir/${patient_id}_filtered_gatk-mutect2.vcf.gz"
        # Guarda las variantes filtradas de Mutect2.

        echo "Saving sanger filtered variants"
        bcftools view -I "$work_dir/tmp/isec/0001.vcf" -O z -o "$p2_dir/${patient_id}_filtered_sanger.vcf.gz"
        # Guarda las variantes filtradas de Sanger.

    fi
    echo "---------------------------------------------------------------------------------"
done

# Annotation and filtration of the VCF files
# A partir de aquí, se realiza la anotación y filtrado de los archivos VCF.

echo "VEP annotation"
# Anotación de variantes utilizando el software VEP.

p_files=$(ls $p2_dir | grep 'filtered_' | grep -v 'tbi' | sed 's/\.vcf\.gz//')
# Encuentra los archivos filtrados para anotarlos.

for file in $p_files
do
    patient=$(echo "$file" | cut -d'_' -f2)
    # Extrae el ID del paciente del nombre del archivo.

    echo "Annotating $file"
    # Muestra un mensaje indicando que se está anotando el archivo VCF.

    ~/ensembl-vep/vep -i "$p2_dir/${file}.vcf.gz" --cache --force --offline --vcf --compress_output gzip --fork 4 --protein -o "$a_dir/${file}_ann.vcf.gz" --sf $a_dir/${file}_ann.html --fa ${m_dir}/GRCh38_full_analysis_set_plus_decoy_hla.fa --gtf ${m_dir}/gencode.v45.annotation.gtf.gz --force_overwrite --terms SO --tsl --biotype --hgvs --plugin Frameshift --plugin Wildtype
    # Ejecuta la anotación VEP en los archivos VCF usando referencias de genoma y diferentes plugins para variantes de genes.

    echo "Filtering: variants of non-expressed or non-coding genes will be removed"
    # Filtra las variantes de genes que no son codificantes o que no se expresan.

    cohort=$(echo "$file" | cut -d '_' -f 1)
    # Determina a qué cohorte pertenece el archivo (AU o CA).

    if [[ $cohort == *APIG_AU* ]]
    then
        genes_to_remove="Exp0_AU.txt" 
    else
        genes_to_remove="Exp0_CA.txt"
    fi
    # Define qué archivo de genes no expresados se utilizará para cada cohorte (AU o CA).

    ~/ensembl-vep/filter_vep -i $a_dir/${file}_ann.vcf.gz -f "(SYMBOL != $genes_to_remove) and (BIOTYPE is protein_coding)" --format vcf --gz -o "$a_dir/${file}_ann_fil.vcf"
    # Filtra los genes no expresados o que no son codificantes y guarda el resultado.

    echo "Recompressing"
    bcftools view -I "$a_dir/${file}_ann_fil.vcf" -O z -o "$a_dir/${file}_ann_fil.vcf.gz"
    # Vuelve a comprimir los archivos filtrados.

    rm "$a_dir/${file}_ann_fil.vcf"
    # Elimina el archivo temporal descomprimido.

    echo "Indexing the filtered file"
    bcftools index -f -t -o "$a_dir/${file}_ann_fil.vcf.tbi" "$a_dir/${file}_ann_fil.vcf.gz"
    # Indexa el archivo VCF filtrado y comprimido.

    echo "Filtering completed"
    # Indica que el filtrado ha terminado.

done

end=`date +%H:%M:%S`
# Guarda la hora de finalización del script.

echo "Execution time from $start to $end"
# Muestra el tiempo de ejecución total del script.
