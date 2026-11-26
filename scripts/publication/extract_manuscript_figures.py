"""Recover exact PNG bytes from the privately held author-review DOCX."""
from pathlib import Path
import argparse,csv,hashlib,zipfile

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manuscript',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();root=Path(__file__).resolve().parents[2]
    out=Path(a.output).resolve();source=Path(a.manuscript)
    if out.exists() or out==root or root in out.parents:p.error('Use a new directory outside the repository.')
    rows=list(csv.DictReader((root/'figures/publication/FIGURE_SOURCE_MAP.csv').open()))
    assert hashlib.sha256(source.read_bytes()).hexdigest()==rows[0]['source_document_sha256'],'Manuscript identity mismatch'
    with zipfile.ZipFile(source) as z:
        data=[(Path(r['file']).name,z.read(r['source_member']),r['sha256']) for r in rows]
    assert all(hashlib.sha256(b).hexdigest()==h for _,b,h in data)
    out.mkdir()
    for name,b,_ in data:(out/name).write_bytes(b)
    print('9 exact manuscript PNG objects recovered. No artwork redrawn.')

if __name__=='__main__':main()
