# Check numbered references against PubMed (NCBI E-utilities: ecitmatch + esummary).
import re, sys, json, time, urllib.request, urllib.parse
md, out = sys.argv[1], sys.argv[2]
s = open(md).read(); refs = s[s.index('## References'):]
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(url):
    for k in range(4):
        try: return urllib.request.urlopen(url, timeout=30).read().decode()
        except Exception as e: time.sleep(2 * (k + 1))
    return ''
rows = []
for m in re.finditer(r'^(\d+)\. (.+)$', refs, re.M):
    n, txt = int(m.group(1)), m.group(2)
    au = re.match(r"\s*([A-Za-zÀ-ÿ'\- ]+?)\s+[A-Z]{1,3}[,.]", txt)
    au = au.group(1).split()[-1] if au else ''
    jm = re.search(r'\*([^*]+)\*\.?\s*(\d{4})[^;]*;\s*(\d+)[^:]*:\s*([A-Za-z]?\d+)', txt)
    title = re.sub(r'\*.*', '', txt.split('. ', 1)[1] if '. ' in txt else txt).strip()
    rec = dict(n=n, ref=txt, author=au)
    pmid = ''
    if jm:
        j, y, v, p = jm.groups(); rec.update(journal=j, year=y, vol=v, page=p)
        q = f'{j}|{y}|{v}|{p}|{au}|k{n}|'
        r = get(E + 'ecitmatch.cgi?db=pubmed&retmode=xml&bdata=' + urllib.parse.quote(q))
        pm = r.strip().split('|')[-1] if r else ''
        if pm.isdigit(): pmid = pm
    if not pmid:  # fall back on a title search
        t = re.sub(r'[^A-Za-z0-9 ]', ' ', title)[:200]
        r = get(E + 'esearch.fcgi?db=pubmed&retmode=json&term=' + urllib.parse.quote(t + '[ti]'))
        try:
            ids = json.loads(r)['esearchresult']['idlist']
            if len(ids) == 1: pmid = ids[0]; rec['via'] = 'title'
        except Exception: pass
    rec['pmid'] = pmid
    if pmid:
        r = json.loads(get(E + f'esummary.fcgi?db=pubmed&retmode=json&id={pmid}'))['result'][pmid]
        rec.update(pm_title=r['title'], pm_first=r['authors'][0]['name'] if r['authors'] else '',
                   pm_journal=r['source'], pm_year=r['pubdate'][:4], pm_vol=r['volume'], pm_pages=r['pages'],
                   doi=next((a['value'] for a in r['articleids'] if a['idtype'] == 'doi'), ''))
    rows.append(rec); time.sleep(0.4)
json.dump(rows, open(out, 'w'), indent=1, ensure_ascii=False)
for r in rows:
    flag = []
    if not r['pmid']: flag.append('NO PMID')
    else:
        if r.get('year') and r['year'] != r['pm_year']: flag.append(f"year {r['year']}≠{r['pm_year']}")
        if r.get('vol') and r['vol'] != r['pm_vol']: flag.append(f"vol {r['vol']}≠{r['pm_vol']}")
        if r.get('page') and not r['pm_pages'].startswith(r['page']): flag.append(f"pages {r['page']}≠{r['pm_pages']}")
        if r['author'] and r['author'].lower() not in r['pm_first'].lower(): flag.append(f"1st author {r['author']}≠{r['pm_first']}")
    print(r['n'], r['pmid'] or '-', '; '.join(flag) or 'ok', '|', r.get('via', ''))
