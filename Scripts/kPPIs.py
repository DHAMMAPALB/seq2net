#!/usr/bin/python

import sys
import requests
from Bio import Entrez

ref=sys.argv[1]
inname=sys.argv[2]
outname=sys.argv[3]
indir=sys.argv[4]
infile="/".join([indir,inname])
ref1=ref.replace(".faa","")
ref2=ref1.replace("_"," ")

Entrez.email="dhammapalb@gmail.com"
handle=Entrez.esearch(db="Taxonomy",term=ref2)
record=Entrez.read(handle)
species_id=record["IdList"][0]

names=[]
for line in open(infile):
	fields=line.split("\t")
	names.append(fields[0])
	names.append(fields[1])

my_genes=set(names)

string_api_url="https://string-db.org/api"
output_format="tsv-no-header"
method="interaction_partners"
request_url="/".join([string_api_url,output_format,method])
params={
	"identifiers" : "%0d".join(my_genes),
	"species" : species_id
}
response=requests.post(request_url, data=params)
ppis=[]
for line in response.text.strip().split("\n"):
	l=line.strip().split("\t")
	if l[2] in my_genes and l[3] in my_genes and float(l[5])>=0.8:
		query_name=l[2]
		partner_name=l[3]
		combined_score=l[5]
		ppis.append("\t".join([query_name,partner_name,combined_score]))
outfile="/".join([indir,outname])
o=open(outfile,"w")
i=0
for line in ppis:
	fields=line.split("\t")
	line2="\t".join([fields[1],fields[0],str(fields[2])])
	if line2 in ppis[i:]:
		print("",end='')
	else:
		print("%s\t%s\t%s\n" % (fields[0],fields[1],str(fields[2])),end='',file=o)
	i=i+1

o.close()

raise SystemExit
