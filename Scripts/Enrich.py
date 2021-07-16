#!/usr/bin/python

import sys
import requests
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez

ref=sys.argv[1]
inname=sys.argv[2]
indir=sys.argv[3]
infile="/".join([indir,inname])
ref1=ref.replace(".faa","")
ref2=ref1.replace("_"," ")

Entrez.email="dhammapalb@gmail.com"
handle=Entrez.esearch(db="Taxonomy", term=ref2)
record=Entrez.read(handle)
species_id=record["IdList"][0]

names=[]
o_in=open(infile)
next(o_in)
for line in o_in:
	fields=line.split("\t")
	names.append(fields[0])
	names.append(fields[1].replace('\n',''))

my_genes=set(names)

request_url="http://pantherdb.org/services/oai/pantherdb/enrich/overrep?-"
params={
	"geneInputList" : ",".join(my_genes),
	"organism" : species_id,
	"annotDataSet" : "ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP",
	"enrichmentTestType" : "FISHER",
	"correction" : "FDR",
	"-H" : "accept: application/json"
}
response=requests.post(request_url, data=params)
outbpj="/".join([indir,"Enriched_BP.json"])
o=open(outbpj,"w")
for line in response.text.strip().split("\n"):
	print(line,file=o)

o.close()

df=pd.read_json(outbpj)
table=df['results']['result']
res=pd.DataFrame(table,columns=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','term'])
f_res=res[res['pValue']<=0.05]
f_res2=pd.DataFrame.from_records(f_res.term)
f_res3=pd.concat([f_res['number_in_list'],f_res['fold_enrichment'],f_res['fdr'],f_res['expected'],f_res['number_in_reference'],f_res['pValue'],f_res2['id'],f_res2['label']],keys=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','GO_id','GO_label'],axis=1)
outbpc="/".join([indir,"Enriched_BP.csv"])
f_res3.to_csv(outbpc,sep="\t",index=False)
df2=pd.DataFrame(res,columns=['number_in_list','term'])
df3=pd.DataFrame.from_records(df2.term)
df4=pd.concat([df2['number_in_list'],df3['label']],keys=['number','GO'],axis=1)
x=df4.number[0:10]
y=df4.GO[0:10]
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots()
fig.set_size_inches(11,7)
ax.invert_yaxis()
plt.barh(y, x, align='center', alpha=0.5,height=0.9)
plt.subplots_adjust(left=0.45, right=0.97, bottom=0.1,top=0.9)
plt.xlabel("Gene count",fontsize=10)
plt.ylabel("Terms",fontsize=10)
plt.title("Enriched_BP",fontsize=12)
outbpp="/".join([indir,"Enriched_BP.pdf"])
plt.savefig(outbpp)

request_url="http://pantherdb.org/services/oai/pantherdb/enrich/overrep?-"
params={
	"geneInputList" : ",".join(my_genes),
	"organism" : species_id,
	"annotDataSet" : "ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF",
	"enrichmentTestType" : "FISHER",
	"correction" : "FDR",
	"-H" : "accept: application/json"
}
response=requests.post(request_url, data=params)
outmfj="/".join([indir,"Enriched_MF.json"])
o=open(outmfj,"w")
for line in response.text.strip().split("\n"):
	print(line,file=o)

o.close()

df=pd.read_json(outmfj)
table=df['results']['result']
res=pd.DataFrame(table,columns=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','term'])
f_res=res[res['pValue']<=0.05]
f_res2=pd.DataFrame.from_records(f_res.term)
f_res3=pd.concat([f_res['number_in_list'],f_res['fold_enrichment'],f_res['fdr'],f_res['expected'],f_res['number_in_reference'],f_res['pValue'],f_res2['id'],f_res2['label']],keys=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','GO_id','GO_label'],axis=1)
outmfc="/".join([indir,"Enriched_MF.csv"])
f_res3.to_csv(outmfc,sep="\t",index=False)
df2=pd.DataFrame(res,columns=['number_in_list','term'])
df3=pd.DataFrame.from_records(df2.term)
df4=pd.concat([df2['number_in_list'],df3['label']],keys=['number','GO'],axis=1)
x=df4.number[0:10]
y=df4.GO[0:10]
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots()
fig.set_size_inches(11,7)
ax.invert_yaxis()
plt.barh(y, x, align='center', alpha=0.5,height=0.9)
plt.subplots_adjust(left=0.45, right=0.97, bottom=0.1,top=0.9)
plt.xlabel("Gene count",fontsize=10)
plt.ylabel("Terms",fontsize=10)
plt.title("Enriched_MF",fontsize=12)
outmfp="/".join([indir,"Enriched_MF.pdf"])
plt.savefig(outmfp)

request_url="http://pantherdb.org/services/oai/pantherdb/enrich/overrep?-"
params={
	"geneInputList" : ",".join(my_genes),
	"organism" : species_id,
	"annotDataSet" : "ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC",
	"enrichmentTestType" : "FISHER",
	"correction" : "FDR",
	"-H" : "accept: application/json"
}
response=requests.post(request_url, data=params)
outccj="/".join([indir,"Enriched_CC.json"])
o=open(outccj,"w")
for line in response.text.strip().split("\n"):
	print(line,file=o)

o.close()

df=pd.read_json(outccj)
table=df['results']['result']
res=pd.DataFrame(table,columns=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','term'])
f_res=res[res['pValue']<=0.05]
f_res2=pd.DataFrame.from_records(f_res.term)
f_res3=pd.concat([f_res['number_in_list'],f_res['fold_enrichment'],f_res['fdr'],f_res['expected'],f_res['number_in_reference'],f_res['pValue'],f_res2['id'],f_res2['label']],keys=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','GO_id','GO_label'],axis=1)
outccc="/".join([indir,"Enriched_CC.csv"])
f_res3.to_csv(outccc,sep="\t",index=False)
df2=pd.DataFrame(res,columns=['number_in_list','term'])
df3=pd.DataFrame.from_records(df2.term)
df4=pd.concat([df2['number_in_list'],df3['label']],keys=['number','GO'],axis=1)
x=df4.number[0:10]
y=df4.GO[0:10]
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots()
fig.set_size_inches(11,7)
ax.invert_yaxis()
plt.barh(y, x, align='center', alpha=0.5,height=0.9)
plt.subplots_adjust(left=0.45, right=0.97, bottom=0.1,top=0.9)
plt.xlabel("Gene count",fontsize=10)
plt.ylabel("Terms",fontsize=10)
plt.title("Enriched_CC",fontsize=12)
outccp="/".join([indir,"Enriched_CC.pdf"])
plt.savefig(outccp)

request_url="http://pantherdb.org/services/oai/pantherdb/enrich/overrep?-"
params={
	"geneInputList" : ",".join(my_genes),
	"organism" : species_id,
	"annotDataSet" : "ANNOT_TYPE_ID_PANTHER_PATHWAY",
	"enrichmentTestType" : "FISHER",
	"correction" : "FDR",
	"-H" : "accept: application/json"
}
response=requests.post(request_url, data=params)
outpathj="/".join([indir,"Enriched_Pathways.json"])
o=open(outpathj,"w")
for line in response.text.strip().split("\n"):
	print(line,file=o)

o.close()

df=pd.read_json(outpathj)
table=df['results']['result']
res=pd.DataFrame(table,columns=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','term'])
f_res=res[res['pValue']<=0.05]
f_res2=pd.DataFrame.from_records(f_res.term)
f_res3=pd.concat([f_res['number_in_list'],f_res['fold_enrichment'],f_res['fdr'],f_res['expected'],f_res['number_in_reference'],f_res['pValue'],f_res2['id'],f_res2['label']],keys=['number_in_list','fold_enrichment','fdr','expected','number_in_reference','pValue','GO_id','GO_label'],axis=1)
outpathc="/".join([indir,"Enriched_Pathways.csv"])
f_res3.to_csv(outpathc,sep="\t",index=False)
df2=pd.DataFrame(res,columns=['number_in_list','term'])
df3=pd.DataFrame.from_records(df2.term)
df4=pd.concat([df2['number_in_list'],df3['label']],keys=['number','GO'],axis=1)
x=df4.number[0:10]
y=df4.GO[0:10]
plt.rcParams.update({'font.size': 7})
fig, ax = plt.subplots()
fig.set_size_inches(11,7)
ax.invert_yaxis()
plt.barh(y, x, align='center', alpha=0.5,height=0.9)
plt.subplots_adjust(left=0.45, right=0.97, bottom=0.1,top=0.9)
plt.xlabel("Gene count",fontsize=10)
plt.ylabel("Terms",fontsize=10)
plt.title("Enriched_Pathways",fontsize=12)
outpathp="/".join([indir,"Enriched_pathways.pdf"])
plt.savefig(outpathp)



