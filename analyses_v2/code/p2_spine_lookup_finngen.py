import pysam
B="https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_"
eps=["M13_INTERVERTEB","M13_DORSOPATHY","M13_LOWBACKPAIN","M13_SCIATICA","M13_SPONDYLOSIS","M13_SPONDYLOPATHY","M13_LOWBACKPAINORANDSCIATICA"]
for e in eps:
    u=B+e+".gz"
    try:
        t=pysam.TabixFile(u,index=u+".tbi")
        h=t.header[-1] if t.header else ""
        hit=[r for r in t.fetch("3",37453000,37454000) if "37453519" in r]
        print(e,"|",h,"|",hit,flush=True)
    except Exception as ex: print(e,"ERR",ex,flush=True)
