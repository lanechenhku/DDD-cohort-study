import pysam
B="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90566001-GCST90567000"
names={795:"All OA",796:"Finger",797:"Hand",798:"Hip",799:"Hip/knee",800:"Knee",801:"Spine",802:"THR",803:"Thumb",804:"TJR",805:"TKR"}
out=open("go2_itga9_region.tsv","w")
for i,n in names.items():
    g=f"GCST90566{i}"; u=f"{B}/{g}/harmonised/{g}.h.tsv.gz"
    t=pysam.TabixFile(u, index=u+".tbi")
    best=None; k=0
    for r in t.fetch("3",37451000,37825000):   # ITGA9 GRCh38 approx
        f=r.split("\t")
        try: p=float(f[7])
        except: continue
        k+=1; out.write(n+"\t"+r+"\n")
        if best is None or p<best[0]: best=(p,f[15],f[2],f[4],f[6])
    print(n, k, best, flush=True)
