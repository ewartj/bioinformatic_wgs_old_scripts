import os

genome = [f for f in os.listdir('.') if f.endswith('.fasta')]
genome = sorted(genome)

#making sure theres a record of whats been annotated

with open('fastas_annotated.tab', 'wb') as f:
	writer = csv.writer(f, dialect = 'excel-tab')
	writer.writerows(izip(names,genome))


count = 0
names = []
for i in genome:
    G1 = genome[count]
    run = G1.split('.')[0]
    print(G1)

    os.system("prokka --outdir prokka/" + run + " --prefix " + run + " " + G1 + " --force")
    count = count + 1

    if count > total:
        print("finished annotation")
