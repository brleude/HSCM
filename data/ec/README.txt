Ec Key Proteins.tsv：Escherichia coli key protein data.
The data format is as follows:
Database	Accession	Data Type	Species	Title	Url	Description	Update Date	Attributes
deg	DEG10480311	gene	[Escherichia coli ST131 strain EC958]	DEG10480311	http://tubic.org/deg/public/index.php/information/bacteria/DEG10480311	50S ribosomal subunit protein L11	2019-08-26T08:18:18.000Z	{"GeneName":"acrF","refseq":"HG941718"}
deg	DEG10480305	gene	[Escherichia coli ST131 strain EC958]	DEG10480305	http://tubic.org/deg/public/index.php/information/bacteria/DEG10480305	hypothetical protein	2019-08-26T08:18:18.000Z	{"GeneName":"yjgQ","refseq":"HG941718"}
deg	DEG10480304	gene	[Escherichia coli ST131 strain EC958]	DEG10480304	http://tubic.org/deg/public/index.php/information/bacteria/DEG10480304	hypothetical protein	2019-08-26T08:18:18.000Z	{"GeneName":"yjgP","refseq":"HG941718"}
deg	DEG10480303	gene	[Escherichia coli ST131 strain EC958]	DEG10480303	http://tubic.org/deg/public/index.php/information/bacteria/DEG10480303	valyl-tRNA synthetase	2019-08-26T08:18:18.000Z	{"GeneName":"valS","refseq":"HG941718"}
......
We extracted the GeneName data from the Attributes column.


Ec PPI Network.csv：Escherichia coli PPI network data.
The data format is as follows:
Protein_A	Protein_B
carB	carA
PRNP	groEL
malE	malG
xerC	xerD
......
Each row of data indicates that the two proteins in the row have an interaction.