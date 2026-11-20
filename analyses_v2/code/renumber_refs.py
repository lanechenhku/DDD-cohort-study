# Renumber Vancouver-style ^n^ citations in order of first citation; drop uncited references.
# Citations are only read after the "## Abstract" heading (so affiliation superscripts are untouched).
import re, sys
p = sys.argv[1]; s = open(p).read()
a = s.index('## Abstract'); r = s.index('## References')
head, body, refs = s[:a], s[a:r], s[r:]
cit = re.compile(r'\^([0-9][0-9,–]*)\^')
def expand(g):
    out = []
    for part in g.split(','):
        if '–' in part: x, y = part.split('–'); out += list(range(int(x), int(y) + 1))
        else: out.append(int(part))
    return out
order = []
for m in cit.finditer(body):
    for n in expand(m.group(1)):
        if n not in order: order.append(n)
new = {old: i + 1 for i, old in enumerate(order)}
def compress(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        out.append(str(nums[i]) if j == i else (f'{nums[i]},{nums[j]}' if j == i + 1 else f'{nums[i]}–{nums[j]}'))
        i = j + 1
    return ','.join(out)
body = cit.sub(lambda m: '^' + compress([new[n] for n in expand(m.group(1))]) + '^', body)
entries = {int(m.group(1)): m.group(2) for m in re.finditer(r'^(\d+)\. (.+)$', refs, re.M)}
missing = [n for n in order if n not in entries]
assert not missing, f'cited but not listed: {missing}'
dropped = sorted(set(entries) - set(order))
lines = [l for l in refs.split('\n') if not re.match(r'^\d+\. ', l)]
while lines and lines[-1].strip() == '': lines.pop()
newrefs = '\n'.join(lines) + '\n\n' + '\n'.join(f'{new[o]}. {entries[o]}' for o in order) + '\n'
open(p, 'w').write(head + body + newrefs)
print(p, 'refs', len(order), 'dropped (uncited):', [(d, entries[d][:60]) for d in dropped])
