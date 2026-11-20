import json, urllib.request, time
loci={"MC.avg":["3:37495010","6:140217584","10:83902678","3:119107610","22:33653104"],
"Core.slope":["11:119276600","10:121281837","9:16252807","6:129715289","3:110087264"],
"HIZ.avg":["9:138154270","9:93597535","7:118571442","7:52334147","17:9889258"],
"DH.avg":["8:5449701","1:92142744","6:124194354","7:54177773","5:91698010"]}
res={}
for ph,ls in loci.items():
    for l in ls:
        c,p=l.split(":"); p=int(p)
        for w in (250000,1000000):
            u=f"https://grch37.rest.ensembl.org/overlap/region/human/{c}:{max(1,p-w)}-{p+w}?feature=gene;content-type=application/json"
            g=json.load(urllib.request.urlopen(u,timeout=60))
            g=[x for x in g if x.get("biotype")=="protein_coding" and x.get("external_name")]
            if g: break
        def dist(x): return 0 if x["start"]<=p<=x["end"] else min(abs(x["start"]-p),abs(x["end"]-p))
        g.sort(key=dist)
        res[l]=[(x["external_name"],dist(x)) for x in g[:2]]
        print(ph,l,res[l],flush=True); time.sleep(0.2)
json.dump(res,open("nearest.json","w"))
