# TFM-Bioinformatica

Repositorio asociado al Trabajo de Fin de Máster:

"Análisis transcriptómico de la respuesta inmune y molecular en un modelo murino de sepsis neonatal"

## Descripción

Este repositorio contiene el script utilizado para el análisis de expresión diferencial y enriquecimiento 
funcional de datos transcriptómicos procedentes de los conjuntos GEO GSE220048 y GSE268868.

## Análisis realizados

- Transformación log2.
- Filtrado de genes con baja expresión.
- Normalización por cuantiles.
- Análisis de componentes principales (PCA).
- Expresión diferencial mediante limma.
- Volcano plots.
- Mapas de calor.
- Identificación de genes compartidos.
- Enriquecimiento funcional Gene Ontology (GO).

## Requisitos

- R 4.5.1
- Bioconductor 3.22

Paquetes principales:

- limma
- clusterProfiler
- org.Mm.eg.db
- enrichplot
- ggplot2
- ggrepel
- pheatmap
- dplyr
- patchwork
- png
- gridExtra
- VennDiagram

## Archivos de entrada

El script requiere los siguientes archivos de expresión:

- GSE220048_spleen_RPKM.txt
- GSE268868_RNAseq_RPKM.txt

Estos archivos deben encontrarse en el directorio de trabajo.

## Ejecución

1. Instalar R y los paquetes necesarios.
2. Colocar el script y los archivos de entrada en el mismo directorio.
3. Establecer ese directorio como directorio de trabajo.
4. Ejecutar el script en R o RStudio.

## Resultados

El script genera archivos de resultados de limma,gráficos, mapas de calor, resultados de enriquecimiento GO
y tablas de genes compartidos.
