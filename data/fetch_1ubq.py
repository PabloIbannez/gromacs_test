#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen


PDB_ID = "1UBQ"
URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"


def main() -> None:
    output = Path(__file__).resolve().with_name(f"{PDB_ID.lower()}.pdb")
    with urlopen(URL, timeout=30) as response:
        content = response.read()
    output.write_bytes(content)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
