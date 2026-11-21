"""Joins panel A and panel B into Figure 5.

Usage: python3 figure5_join.py A.png B.png A.pdf B.pdf OUT_PREFIX

PNG: the two images stacked (1950 px wide, 300 dpi). PDF: the two panel pages scaled to a width of 468 pt and placed on
one page. Writes OUT_PREFIX.png and OUT_PREFIX_vector.pdf by exclusive create."""
import sys
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation, PageObject
A_PNG, B_PNG, A_PDF, B_PDF, OUT = sys.argv[1:6]
a, b = (Image.open(p).convert("RGB") for p in (A_PNG, B_PNG))
assert a.size[0] == b.size[0] == 1950, (a.size, b.size)
full = Image.new("RGB", (1950, a.size[1] + b.size[1]), "white")
full.paste(a, (0, 0)); full.paste(b, (0, a.size[1]))
with open(OUT + ".png", "xb") as f:
    full.save(f, "PNG", dpi=(300, 300))
pa, pb = (PdfReader(p).pages[0] for p in (A_PDF, B_PDF))
(wa, ha), (wb, hb) = ((float(p.mediabox.width), float(p.mediabox.height)) for p in (pa, pb))
sa, sb = 468 / wa, 468 / wb
page = PageObject.create_blank_page(width=468, height=ha * sa + hb * sb)
page.merge_transformed_page(pa, Transformation().scale(sa).translate(0, hb * sb))
page.merge_transformed_page(pb, Transformation().scale(sb).translate(0, 0))
w = PdfWriter(); w.add_page(page)
with open(OUT + "_vector.pdf", "xb") as f:
    w.write(f)
print("written:", OUT + ".png", full.size, "|", OUT + "_vector.pdf", "468 x %.1f pt" % (ha * sa + hb * sb))
