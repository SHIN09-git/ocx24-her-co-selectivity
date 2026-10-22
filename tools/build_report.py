"""Optional presentation build: pip install -r requirements-report.txt."""
from pathlib import Path
import base64
import html
import re
import markdown

ROOT=Path(__file__).resolve().parents[1]
source=ROOT/"05_reports/首期计算研究报告.md"
body=markdown.markdown(source.read_text(encoding="utf-8"),extensions=["tables","fenced_code","toc"])
def image(match):
    relative=html.unescape(match.group(1))
    path=(source.parent/relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise ValueError("Image outside release or missing: "+relative)
    return 'src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'"'
body=re.sub(r'src="([^"]+\.png)"',image,body)
body=body.replace("<table>","<div class=\"table-scroll\"><table>").replace("</table>","</table></div>")
css="""
:root{color-scheme:light}*{box-sizing:border-box}
body{margin:0;background:#f5f7fa;color:#172b3a;font:17px/1.8 system-ui,"Microsoft YaHei","Noto Sans CJK SC",sans-serif}
main{max-width:1060px;margin:40px auto;padding:44px 54px;background:white;border:1px solid #dce3e9;border-radius:8px}
h1{font-size:31px;line-height:1.4;margin:0 0 24px}
h2{font-size:24px;margin:40px 0 16px;border-top:1px solid #dde4eb;padding-top:24px}
h3{font-size:20px}p{margin:16px 0}a{color:#155e91;text-underline-offset:3px}
img{display:block;max-width:100%;height:auto;margin:24px auto}
table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.6;white-space:normal}
th,td{border-bottom:1px solid #dbe3e9;padding:10px 12px;text-align:left;vertical-align:top}
th{background:#eef3f7}tbody tr:nth-child(even){background:#fafbfd}.table-scroll{overflow-x:auto}
code{background:#edf2f6;padding:2px 5px;border-radius:3px;font-size:.9em;overflow-wrap:anywhere}
pre{background:#edf2f6;padding:16px;overflow:auto}li{margin:8px 0}
footer{font-size:14px;color:#526573;margin-top:40px;border-top:1px solid #dce3e9;padding-top:18px}
@media(max-width:700px){main{margin:0;padding:26px 18px;border:0;border-radius:0}body{font-size:16px}h1{font-size:26px}}
@media print{body{background:white}main{max-width:none;border:0;margin:0;padding:0}h2,h3{break-after:avoid}img,table{break-inside:avoid}a{color:inherit}}
"""
out='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>OCx24：HER 前测对 CO 选择性预测的增量检验</title><style>'+css+'</style></head><body><main>'+body+'<footer>JIN HANZHE · 公开版整理 2026-10-22 · CC BY 4.0。图像已内嵌；证据文件链接相对本发布包解析。</footer></main></body></html>\n'
(ROOT/"05_reports/report.html").write_text(out,encoding="utf-8",newline="\n")
print("Generated report.html with seven embedded archived figures.")

