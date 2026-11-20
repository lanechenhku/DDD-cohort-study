"""Build docx + pdf from a living md file with every [ ... TODO ... ] item highlighted in yellow (Lane, 2026-11-07).
The md file itself is NOT changed: a temporary copy wraps each bracketed TODO in a pandoc span with the character style
"TODO" (bold, yellow highlight) defined in todo_reference.docx (pandoc's default reference.docx + that one style).
Usage (from /mnt/project-files/DDD_papers):  python3 analyses_v2/code/build/build_doc.py md/X.md 2026-11-07
Writes docx/X_<date>.docx and pdf/X_<date>.pdf. Also importable: todo_items(path) lists the TODOs for the audit doc."""
import os, re, subprocess, sys, tempfile, shutil
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "todo_reference.docx")

def _brackets(text):
    """Yield (start, end) of outermost balanced [...] groups whose content is a TODO item."""
    i, n = 0, len(text)
    while i < n:
        if text[i] == "[" and (i == 0 or text[i - 1] not in "\\!"):
            depth, j = 0, i
            while j < n:
                if text[j] == "\\": j += 2; continue
                if text[j] == "[": depth += 1
                elif text[j] == "]":
                    depth -= 1
                    if depth == 0: break
                j += 1
            if j < n:
                inner = text[i + 1:j]
                nxt = text[j + 1:j + 2]
                if nxt not in ("(", "{") and re.search(r"(^TODO\b)|(\bTODO(\b|:|\s*\())", inner) and not re.search(r"with TODO$", inner):
                    yield i, j + 1; i = j + 1; continue
        i += 1

def highlight(text):
    out, last = [], 0
    for a, b in _brackets(text):
        out.append(text[last:a]); out.append("[\\[" + text[a + 1:b - 1] + "\\]]{custom-style=\"TODO\"}"); last = b
    out.append(text[last:]); return "".join(out)

def todo_items(path):
    """[(section heading, TODO text)] in document order."""
    text = open(path, encoding="utf-8").read(); items = []
    for a, b in _brackets(text):
        head = [l for l in text[:a].splitlines() if l.startswith("#")]
        sec = head[-1].lstrip("#").strip() if head else "(front matter / draft notes)"
        items.append((sec, re.sub(r"\s+", " ", text[a + 1:b - 1]).strip()))
    return items

def build(md, date):
    name = os.path.splitext(os.path.basename(md))[0]
    tmpd = tempfile.mkdtemp()
    tmp = os.path.join(tmpd, name + ".md"); open(tmp, "w", encoding="utf-8").write(highlight(open(md, encoding="utf-8").read()))
    docx = os.path.join(ROOT, "docx", f"{name}_{date}.docx")
    subprocess.run(["pandoc", tmp, f"--resource-path={ROOT}", f"--reference-doc={REF}", "-o", docx], check=True)
    subprocess.run(["soffice", "--headless", "--convert-to", "pdf", "--outdir", tmpd, docx], check=True, capture_output=True)
    shutil.copy(os.path.join(tmpd, f"{name}_{date}.pdf"), os.path.join(ROOT, "pdf", f"{name}_{date}.pdf"))
    shutil.rmtree(tmpd); print("built", docx)

if __name__ == "__main__":
    build(sys.argv[1], sys.argv[2])
