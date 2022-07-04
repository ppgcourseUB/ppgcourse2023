# Bioinformatics primer and data formats

Instructor: **Cinta Pegueroles Queralt**


In this session we will learn some basic bash commands. First we will create the folders architecture that we will use in the two sessions:

```
mkdir session1
cd session1
```

### Downloading data:
We will play parsing a fasta and gtf files. First we will download the data from the genome of the Giant panda (<em>Ailuropoda melanoleuca/em>):

```
wget http://ftp.ensembl.org/pub/release-106/fasta/ailuropoda_melanoleuca/pep/Ailuropoda_melanoleuca.ASM200744v2.pep.all.fa.gz

wget http://ftp.ensembl.org/pub/release-106/gtf/ailuropoda_melanoleuca/Ailuropoda_melanoleuca.ASM200744v2.106.gtf.gz
```
    
### Creating and deleting files:
With command line we can create a file, read it and remove it:

```
touch test.txt
gedit test.txt &
cat test.txt
rm test.txt
```

### Example 1: fasta file
We will uncompress the fasta file and we will store the first 1000 lines in a new file (WARNING: last sequence may be truncated!):

```
gunzip -c Ailuropoda_melanoleuca.ASM200744v2.pep.all.fa.gz | head -n 1000 > testProt.fa
```
How many proteins contains this file?

```
grep '>' testProt.fa |wc -l
```
Obtain a file with the list of all identifiers

```
grep '>' testProt.fa |cut -d ' ' -f1 |sort|uniq >  Ids.list
```
Can you compress the new file obtained?

```
gzip testProt.fa 
```
Can you uncompress it again?

```
gunzip testProt.fa 
```

### Example 1: gtf file
    
Copy the first 1000 lines of the GTF file in a folder 
```
gunzip -c Ailuropoda_melanoleuca.ASM200744v2.106.gtf.gz | head -1000 > test.gtf 
```
How many genes contains this file?
```
cut -f3 test.gtf |sort|uniq -c
```
How many exons contain the ENSAMEG00000005298 gene?
```
grep 'ENSAMEG00000005298' test.gtf |grep 'exon\t'
grep 'ENSAMEG00000005298' test.gtf |grep --perl-regex 'exon\t'
grep 'ENSAMEG00000005298' test.gtf |grep --perl-regex 'exon\t'|wc -l
```

### My first script: bash and loops:
We will create your first bash script. Open a text file (for instance gedit or nano in the terminal) and copy the following lines (IMPORTANT: mind the intendents!).

```
#!/bin/bash

for i in $(ls)
	do echo $i
done
```
    
to run the script use the following command:
```
bash myScript.sh
```

Alternatively you could run this command directly on the terminal:
    
```
for i in $(ls); do echo "my file is: "$i; done
```
### AWK    
    
Awk is a programming language designed for text processing. Here some examples:

Print lines within a range:
```
awk 'NR>528 && NR <537' testProt.fa
```
Convert a multi-line fasta to a singleline fasta:
```
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' testProt.fa > testProt_singleline.fa
```
Check if it worked:
```
head -2 testProt*
```
