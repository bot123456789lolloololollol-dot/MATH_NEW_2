#!/usr/bin/env python3
"""Query the OEIS JSON search API and print compact results.

Usage: python oeis_search.py "query words" [max_results]
"""
import json
import sys
import urllib.parse
import urllib.request

UA = "Mozilla/5.0 research-session-S4"


def main() -> None:
    query = sys.argv[1]
    max_results = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    url = ("https://oeis.org/search?q="
           + urllib.parse.quote(query)
           + "&fmt=json")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=30) as r:
        results = json.load(r)
    if isinstance(results, dict):
        print("NO RESULTS:", results.get("error", results))
        return
    for s in results[:max_results]:
        print(f"--- A{s['number']:06d} | {s.get('name','')}")
        data = s.get("data", "")
        terms = data.split(",") if data else []
        head = ", ".join(terms[:24])
        tail = ", ".join(terms[-12:]) if len(terms) > 30 else ""
        print(f"    data[{len(terms)}]: {head}" + (f" ... {tail}" if tail else ""))
        for key in ("comment",):
            for c in s.get(key, [])[:2]:
                print(f"    {key}: {c[:280]}")
        print(f"    revtime: {s.get('time','')}  refs:{s.get('references','')}")


if __name__ == "__main__":
    main()
