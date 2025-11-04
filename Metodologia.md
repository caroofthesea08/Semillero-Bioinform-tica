# **Metodología**
### **1. Identificar y descargar el gen de la escualeno sintasa en una especie de tiburón**

```
# Buscar la secuencia de proteínas en la base de datos de la ncbi

!pip install biopython

from Bio import Entrez
from Bio import SeqIO

Entrez.email = "xxx@xxx.xx"

with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_041071437.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

print("La longitud de la secuencia es de", len(seq_record.seq), "aminoácidos")
print("La secuencia", seq_record.description, " es: ")
print(seq_record.seq)
SeqIO.write(seq_record, "whiteshark_sqs.fa", "fasta")
```

### **2. Comparación del gen de la escualeno sintasa entre especies de tiburones**

Esto se realiza porque la secuencia elegida está establecida como "squalene synthase-like", por lo que se debe confirmar si otros tiburones sí tienen una secuencia similar que codifica para lo mismo.

- Realizar un blast para buscar las secuencias de otras especies de tiburones, que deben ser las más similares al query elegido.

<img width="887" height="340" alt="image" src="https://github.com/user-attachments/assets/77b5a3a7-4d5c-4fd5-9a33-68476a5724e4" />

Figura 3. Secuencias de tiburones más similares a la secuencia query, a las que se les hará alineamiento

- Realizar alineamiento de la secuencia query junto con las 6 secuencias elegidas de tiburones

```
#Secuencia squalene synthase de Scyliorhinus canicula
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_038673545.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "scanicula_sqs.fa", "fasta")

#Secuencia squalene synthase isoform X2 de Chiloscyllium plagiosum
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_043551139.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "cplagiosum_sqs.fa", "fasta")

#Secuencia squalene synthase isoform X1 de Stegostoma tigrinum
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_048392098.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "stigrinum_sqs.fa", "fasta")

#Secuencia squalene synthase isoform X2 de Hemiscyllium ocellatum
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_060686698.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "hocellatum.fa", "fasta")

#Secuencia squalene synthase isoform X1 de Scyliorhinus torazame
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_072371610.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "storazame_sqs.fa", "fasta")

#Secuencia squalene synthase isoform X2 de Chiloscyllium punctatum
with Entrez.efetch(
    db="protein", rettype="fasta", retmode="text", id="XP_072436270.1") as handle:
    seq_record = SeqIO.read(handle, "fasta")

SeqIO.write(seq_record, "cpunctatum_sqs.fa", "fasta")
```
- Agrupar las secuencias en un archivo de texto y realizar un alineamiento con clustal

```
files = ["whiteshark_sqs.fa", "scanicula_sqs.fa", "cplagiosum_sqs.fa", "stigrinum_sqs.fa", "hocellatum.fa", "storazame_sqs.fa", "cpunctatum_sqs.fa"]

output = "sqs_seqs.txt"

with open(output, "w") as out_handle:
    for file in files:
        for record in SeqIO.parse(file, "fasta"):
            SeqIO.write(record, out_handle, "fasta")

print(f"Archivo combinado creado: {output}")
```

```
!sudo apt-get install clustalw

!clustalw -INFILE=sqs_seqs.txt -TYPE=PROTEIN
```

- Observar el alineamiento y analizar si hay mucha diferencia entre especies similares (si no tienes descargado ClustalX, puedes importar tu archivo sqs_seqs.txt en https://www.ebi.ac.uk/jdispatcher/msa/clustalo)

```
clustalx sqs_seqs.aln
```

<img width="1803" height="670" alt="image" src="https://github.com/user-attachments/assets/c34839ea-250a-446a-b423-8e778ca6aee9" />

Figura 4. Interfaz web de clustal


- Entrar a https://www.ebi.ac.uk/interpro/search/sequence/ y pegar en la casilla la secuencia de aminoácidos de las proteínas elegidas para determinar si los dominios de estas corresponden con el verdadero funcionamiento de la proteína, y es el mismo entre las 6 especies elegidas.

<img width="1672" height="744" alt="image" src="https://github.com/user-attachments/assets/b87d34d2-5e92-430e-bbdf-421a76a3adc1" />


Figura 5. Interfaz web de interpro

### **3.	Alineamiento de la escualeno sintasa en diversos organismos**

Se realiza el mismo flujo de trabajo que en la sección 2, pero con organismos diferentes a tiburones (plantas, algas, levaduras y bacterias)

### **4.	Realizar árbol filogenético para determinar la distancia entre los organismos**

```
!sudo apt-get install -y iqtree

!iqtree2 -s new_sqs_seqs.aln -bb 1000
```
- Para observar el árbol filogenético:
```
from Bio import Phylo

tree = Phylo.read("new_sqs_seqs.aln.treefile", "newick")
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
Phylo.draw(tree, axes=ax)
```

- Para observar la matriz de distancias
```
!cat new_sqs_seqs.aln.mldist
```

### **5.	Selección de un organismo con una proteína similar**
- Modelar las proteínas en https://swissmodel.expasy.org/
- Comparar si su estructura su plegamiento es similar

