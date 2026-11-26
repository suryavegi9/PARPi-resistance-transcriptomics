"""Download manifest-listed GEO inputs into a NEW external directory; no analysis."""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,urllib.request

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    target=Path(args.output).resolve()
    if target==root or root in target.parents or target.exists():
        parser.error('Output must be new and outside the repository.')
    if not target.parent.is_dir(): parser.error('Output parent must exist.')
    rows=list(csv.DictReader((root/'metadata/DATA_ACQUISITION_MANIFEST.csv').open(encoding='utf-8')))
    for r in rows:
        dest=target/r['destination']
        if target not in dest.resolve().parents: raise ValueError('Unsafe destination')
        if not r['url'].startswith('https://ftp.ncbi.nlm.nih.gov/geo/series/'): raise ValueError('Unexpected source')
    target.mkdir()
    report=[]
    for r in rows:
        request=urllib.request.Request(r['url'],headers={'User-Agent':'PARPi-release-acquisition/1.0'})
        with urllib.request.urlopen(request,timeout=120) as response: data=response.read()
        if r['transformation']=='gunzip': data=gzip.decompress(data)
        elif r['transformation']!='none': raise ValueError('Unknown transformation')
        digest=hashlib.sha256(data).hexdigest()
        passed=digest==r['sha256'] and len(data)==int(r['bytes'])
        report.append(dict(accession=r['accession'],destination=r['destination'],sha256=digest,bytes=len(data),pass_checksum=passed))
        (target/'ACQUISITION_REPORT.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        if not passed: raise ValueError('Frozen input checksum mismatch: '+r['filename'])
        dest=target/r['destination'];dest.parent.mkdir(parents=True,exist_ok=True)
        with dest.open('xb') as out: out.write(data)
        print(r['accession'],r['filename'],'PASS',flush=True)
    print('All inputs verified; no analyses executed.')

if __name__=='__main__': main()
