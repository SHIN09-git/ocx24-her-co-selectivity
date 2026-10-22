"""Maintainer-only release freeze. Never use this to hide a failed integrity check."""
from pathlib import Path
import hashlib, json
ROOT=Path(__file__).resolve().parents[1]
records=[]
for p in sorted(ROOT.rglob("*")):
    relative=p.relative_to(ROOT)
    if any(x in (".git",".venv","reproductions","__pycache__") for x in relative.parts):
        continue
    if not p.is_file() or relative.as_posix()=="RELEASE_MANIFEST.json":
        continue
    records.append({"path":relative.as_posix(),"bytes":p.stat().st_size,
                    "sha256":hashlib.sha256(p.read_bytes()).hexdigest()})
(ROOT/"RELEASE_MANIFEST.json").write_text(
    json.dumps({"schema_version":1,"scope":"All release files except this self-referential manifest.",
                "files":records},ensure_ascii=False,indent=2)+"\n",encoding="utf-8",newline="\n")
print(f"Recorded {len(records)} release files.")

