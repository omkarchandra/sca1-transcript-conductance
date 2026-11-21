"""Writes Figure 6 at its final size from the panel file.

Usage: python3 figure6_page.py PANELS.png PANELS.pdf OUT_PREFIX

PNG: the panel image with its resolution set to 300 dpi. PDF: the panel page scaled to a width of 468 pt.
Writes OUT_PREFIX.png and OUT_PREFIX_vector.pdf by exclusive create."""
import sys
from PIL import Image
from pypdf import PdfReader, PdfWriter, Transformation, PageObject
P_PNG, P_PDF, OUT = sys.argv[1:4]
im = Image.open(P_PNG).convert("RGB")
assert im.size[0] == 1950, im.size
with open(OUT + ".png", "xb") as f:
    im.save(f, "PNG", dpi=(300, 300))
p = PdfReader(P_PDF).pages[0]
w0, h0 = float(p.mediabox.width), float(p.mediabox.height)
s = 468 / w0
page = PageObject.create_blank_page(width=468, height=h0 * s)
page.merge_transformed_page(p, Transformation().scale(s))
w = PdfWriter(); w.add_page(page)
with open(OUT + "_vector.pdf", "xb") as f:
    w.write(f)
print("written:", OUT + ".png", im.size, "|", OUT + "_vector.pdf", "468 x %.1f pt" % (h0 * s))
