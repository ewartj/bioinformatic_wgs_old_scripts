# made at APHA#

import os


#
# read 1 = The forward read of the DNA Seq run.
# read 2 = The reverse read of the DNA Seq run.
# reference = The reference genome to align the reads to, this needs to be a nucleotide fasta file.

# need to get list of files done in alphabetical order

fastq_R1 = [f for f in os.listdir('.') if f.endswith('_R1_001.fastq.gz')]
fastq_R1 = sorted(fastq_R1)
fastq_R2 = [f for f in os.listdir('.') if f.endswith('_R2_001.fastq.gz')]
fastq_R2 = sorted(fastq_R2)

count = 0
names = []
total = len(fastq_R1)
for i in fastq_R2:
    R1 = fastq_R1[count]
    run = R1.split('_')[0]
    names.append(run)
    count = count + 1

    if count > total:
        end

print("Total isolates:")
print(total)
print("Sequeces:")
print(names, fastq_R1, fastq_R2)

# Snippy

reference = 'WA1.gb'

# sequences input into snippy

count = 0
names = []
for i in fastq_R2:
    one = count
    R1 = fastq_R1[count]
    R2 = fastq_R2[count]
    run = R1.split('_')[0]
    os.system("snippy --cpus 4 --reference " + reference + " --outdir " + run + "_snippy" + " --pe1 " + R1 + " --pe2 " + R2 + " --force")
    count = count + 1
print(names)
