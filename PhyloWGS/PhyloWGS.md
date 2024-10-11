Este software se utiliza para creal filogenias (clones), para poder utilizar 
el software para predecir la calidad de los neoantígenos.

Como input en nuestro caso vamos a utilizar archivos VCF que deben ser
formateados con el script parser (https://github.com/morrislab/phylowgs/tree/master/parser)
Para esto necesitamos un VCF del hospital

# INSTALACION

- Hay que instalar Python 2, NumPy, SciPy y PyVCF (este ultimo puede ser instalado via pip:
```
pip2 install --user pyvcf
```

** Crear ssm_data.txt desde archivo sanger**
```
./create_phylowgs_inputs.py -s 5000 --vcf-type sample1=sanger sample1=sample.vcf
```
**-s 5000: Esto especifica que se está usando el argumento --sample-size 
para submuestrear variantes somáticas. En este caso, se tomarán un máximo
de 5000 variantes en lugar de procesar todas las variantes del archivo VCF.
Esto es útil para reducir el tiempo de ejecución, especialmente si se trabaja
con un gran número de variantes.**

