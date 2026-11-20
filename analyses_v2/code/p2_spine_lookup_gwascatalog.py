import pysam, urllib.request
def gcst_url(g):
    n=int(g[4:]); lo=(n-1)//1000*1000+1; hi=lo+999
    w=8 if len(g)>9 else 6
    rng=f"GCST{lo:0{len(g)-4}d}-GCST{hi:0{len(g)-4}d}"
    return f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{rng}/{g}/harmonised/{g}.h.tsv.gz"
studies={"GCST90129430":"Chronic back pain (EUR)","GCST90018797":"Back pain (EUR, Sakaue 2021)","GCST90018577":"Back pain (EAS, BBJ)","GCST90859374":"Sciatica (EUR)","GCST90077887":"Sciatica (UKB)","GCST90077923":"Back pain (UKB)"}
for g,n in studies.items():
    u=gcst_url(g)
    try:
        t=pysam.TabixFile(u,index=u+".tbi")
        hdr=t.header[-1] if t.header else ""
        rows=[r for r in t.fetch("3",37453000,37454000)]
        hit=[r for r in rows if "rs1549596" in r or "\t37453519\t" in r]
        print(g,n,"|",hdr[:300],"|",hit,flush=True)
    except Exception as e: print(g,n,"ERR",u,e,flush=True)
