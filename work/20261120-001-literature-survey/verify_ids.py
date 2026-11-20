"""Re-check every arXiv ID cited in literature/*.md against the arXiv API.

Usage (from repo root): python tools/pin.py python work/20261120-001-literature-survey/verify_ids.py
Prints one line per ID: id | title | first authors. IDs the API does not return are listed as MISSING.
Network only, no compute.
"""
import pathlib
import re
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[2]
NEW = r"\d{4}\.\d{4,5}"
OLD = r"(?:hep-th|hep-lat|hep-ph|quant-ph|cond-mat|math-ph)/\d{7}"
PAT = re.compile(rf"arXiv:({NEW}|{OLD})")


def main() -> None:
    ids = set()
    for f in sorted((ROOT / "literature").glob("*.md")):
        ids.update(PAT.findall(f.read_text(encoding="utf8")))
    ids = sorted(ids)
    found = {}
    for k in range(0, len(ids), 25):
        chunk = ids[k:k + 25]
        url = ("https://export.arxiv.org/api/query?max_results=50&id_list=" + ",".join(chunk))
        xml = urllib.request.urlopen(url, timeout=60).read().decode("utf8")
        for e in xml.split("<entry>")[1:]:
            m = re.search(r"<id>http://arxiv.org/abs/(.*?)v\d+</id>", e)
            if not m:
                continue
            title = " ".join(re.search(r"<title>(.*?)</title>", e, re.S).group(1).split())
            authors = re.findall(r"<name>(.*?)</name>", e)[:3]
            found[m.group(1)] = f"{title} | {', '.join(authors)}"
        time.sleep(3)  # arXiv API etiquette
    for i in ids:
        print(f"{i} | {found.get(i, 'MISSING')}")
    print(f"\n{len(ids)} IDs cited, {sum(i in found for i in ids)} returned by the API")


if __name__ == "__main__":
    main()
