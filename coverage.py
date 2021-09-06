import glob, os
import pandas as pd
from datetime import datetime


"""
Author: EJS
NB: will need to specify if fatsqs have 'L001' (MiSeq) or not
"""
###############
##FOLDERS etc##
###############

path = os.getcwd()
path = path + "/"

#getting the date
time = datetime.now().strftime('%d_%m_%Y_%H')

#make all the folders that are needed, if they don't already exist
folders = ['bam','fasta','qualimap','other','sam','assembly_map','kraken']
[os.mkdir(x) for x in folders if not os.path.isdir(x)] 
#############
##FUNCTIONS##
#############

#kraken/braken

def krkan(fastq_R1):
    R1 = str(fastq_R1)
    run = R1.split('_')[0]
    print('Working on: ' + run)
    os.system("kraken --threads 10 --db /mnt/VM0516/minikraken_20171101_8GB_dustmasked " + R1 + " > " + run + ".output")
    os.system("kraken-report --db /mnt/VM0516/minikraken_20171101_8GB_dustmasked " + run + ".output" + " > " + run + ".report" )
    os.system("python /home/p938497/Programs/Bracken-master/src/est_abundance.py -i " + run + ".report -k /mnt/VM0516/minikraken_20171101_8GB_dustmasked/minikraken_8GB_100mers_distrib.txt -o " + run + ".braken -t 17")

#assembling the genome
#If you want spades change the # but you will need to modify list of contig lengths to '	NODE_'
def spades(R1):
    run = str(R1).split('_')[0]
    S = str(R1).split('_')[1]
    name = run + '_' + S
    read1 = str(name) + '_L001_R1_001.fastq.gz'
    read2 = str(name) + '_L001_R2_001.fastq.gz'
    os.system('shovill --outdir ' + run + ' --force --trim --cpus 10 --R1 ' + read1 + ' --R2 '  + read2)
    os.system("mv " + run + "/contigs.fa " + path + run + ".fasta")
    os.system("mv " + run + "/contigs.gfa " +  path + 'assembly_map/' + run + ".gfa") 

#ALL fasta files need to be in same folder as raw reads

def cov(i):
    ID = str(i).split('_')[0]
    S = str(i).split('_')[1]
    name = ID + '_' + S
    read1 = str(name) + '_L001_R1_001.fastq.gz'
    read2 = str(name) + '_L001_R2_001.fastq.gz'
    print(ID)
    os.system("bwa index " + ID + ".fasta")
    print("bwa index " + ID + ".fasta")
    os.system("bwa mem -t 10 " + ID + ".fasta " +  read1 + " " +  read2 + " > " +  ID + ".sam")
    os.system("samtools faidx " + ID + ".fasta")
    os.system("samtools import " + ID + ".fai " + ID + ".sam " +  ID + ".bam")
    os.system("samtools sort -T /tmp/" + ID + ".sorted.bam -o " + ID + ".sorted.bam " + ID + ".bam")
    os.system("samtools index " + ID + ".sorted.bam")
    os.system("qualimap bamqc -bam " + ID + ".sorted.bam -outdir " + path + ID + ".qualimap_results")
    os.system("mv " + path + ID + ".qualimap_results/genome_results.txt " + path + ID + "_results.txt")

# rename all genome_results.txt to ID + '_results.txt')



###############
##LISTS/LOOPS##
###############

#making a list of fastq files
fastq_R1 = [f for f in os.listdir('.') if f.endswith('_L001_R1_001.fastq.gz')]
fastq_R1 = sorted(fastq_R1)
print(fastq_R1)

#running kraken
names = []
for i in fastq_R1:
    krkan(i)

#making database of top kraken hit
Kraken_list = []
percent_list = []
names_list = []
krk = [f for f in os.listdir('.') if f.endswith('.braken')]
krk = sorted(krk) 

for i in krk: 
        file_ = path + str(i)
        name = str(i).split('.braken')[0]
        df = pd.read_csv(file_,sep = '\t')
        df.sort_values(by=['fraction_total_reads'],ascending = False, inplace=True)
        df = df.reset_index(drop=True)
        kraken = df.at[0,'name']
        percentage= df.at[0,'fraction_total_reads']
        Kraken_list.append(kraken)
        names_list.append(name)
        percent_list.append(percentage)





#ASSEMBLY
#making a list of fastq files
for i in fastq_R1:
    spades(i)

#coverage
for i in fastq_R1:
    cov(i)
#making a list of fqualimapped files
quali = [f for f in os.listdir('.') if f.endswith('_results.txt')]
quali = sorted(quali)

cover = []
gcs = []
conts = []
reads = []
lst = []
name_list = []
contig = []
ref=[]
N50 = []

for i in quali:
    name = str(i).split('_')[0]
    name_list.append(name)
    print(name_list)
    with open(i,'rt') as f:
        for line in f:
        #coverage
            if line.startswith('     mean coverageData'):
                tx=float(line.strip()[:-1].split(' = ')[1])
                cover.append(tx)
        #gc content
            if line.startswith('     GC percentage'):
                gc=float(line.strip()[:-1].split(' = ')[1])
                gcs.append(gc)
        #total contigs
            if line.startswith('     number of contigs'):
                cont=int(line.strip().split(' = ')[1])
                conts.append(cont)
        ##total reads
            if line.startswith('     number of reads'):
                read=int(line.strip().split(' = ')[1].replace(',',''))
                lst.append(read)
        #reference genome size for n50
            if line.startswith('     number of bases'):
                #tot = line.decode("utf-8")
                tot=line.strip().split(' = ')[1]
                tot=int(tot.strip().split(' bp')[0].replace(',',''))
                ref.append(tot)
                print('reference size:' + str(tot))
        #list of contig lengths
            if line.startswith('	contig0') == True:
                ln=line.strip().split('\t')[1]
                ln=int(ln.strip().split('\t')[0])
                contig.append(ln)
                contig.sort(reverse = True)
#getting the N50
    nfifty = tot/2
    summ = 0
    for i in contig:
        summ = summ + i
        if summ >= nfifty:
            print('n50:' + str(i))
            n50 = int(i)
            N50.append(n50)
            break

# take all the lists (names, coverage and kraken/braken results) and put them into a csv file
        
df_k = pd.DataFrame(list(zip(names_list,Kraken_list,percent_list)),columns=['Name','Species','Percent Reads Mapped']) 
df_k.set_index('Name')
#df_k.to_csv('kraken' + time + '.csv')
# coverage csv file
df = pd.DataFrame(list(zip(name_list,cover,gcs,conts,lst,ref,N50)),columns=['Name','Coverage','GC content', 'Total Contigs','total reads','Reference genome size','N50']) 
df.set_index('Name')
print(df_k)
print(df) 
#df.to_csv('quality' + time + '.csv')


############
##CLEAN UP##
############

os.system("mv *.fasta fasta/")
os.system("mv *.txt qualimap/")
os.system("mv *.bam bam/")
os.system("mv *.sam sam/")
os.system("mv *.amb other/") 
os.system("mv *.ann other/") 
os.system("mv *.bwt other/") 
os.system("mv *.fai other/") 
os.system("mv *.pac other/") 
os.system("mv *.sa other/")
os.system("mv *.report kraken/")
os.system("mv *.braken kraken/")
os.system("mv *.output kraken/")


