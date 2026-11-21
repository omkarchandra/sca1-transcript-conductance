"""Build p16_score_all_v01.py from p15_score_all_v05.py (bfa035d6...) by counted replacements.
Usage: python3 -I -B build_p16_score_all_v03.py AREA_TSV OUT_PY P15_SCORE_ALL_V05_PY P16_VARIANT_SCORE_V02_PY
AREA_TSV: six lines 'A<TAB>factor sha256<TAB>name of the wave-8 result folder' for the membrane-area versions, A = 0.9, 0.85,
0.8, 0.75, 0.7, 0.6 in this order; each hash must be the one of the wave-9 factor file of that A. OUT_PY is written by
exclusive create. v03 = v02 with these checks of AREA_TSV and of the two source files; the script it builds is unchanged."""
import hashlib, re, sys
if len(sys.argv) != 5:
    sys.exit(__doc__)
src = open(sys.argv[3]).read()
if hashlib.sha256(src.encode()).hexdigest() != "bfa035d60feb56440311a5a3390ef2553ca79e6bd61900cb9ee159fa437b34c1":
    sys.exit("STOP: P15_SCORE_ALL_V05_PY does not have its recorded hash")
scorer_sha = hashlib.sha256(open(sys.argv[4], "rb").read()).hexdigest()
if scorer_sha != "1e0a200d7d6a6e066e6d48796b84fec47f527a5b37a4185274a950fd184d9a5e":
    sys.exit("STOP: P16_VARIANT_SCORE_V02_PY does not have its recorded hash")
PIN = {"0.9": "25264b8faa433761c2468e04a2947351c3a438d50c7d3670f8801eb1b8d9f5c7", "0.85": "7253210c3c233de005824a6374235c6052925c200e045f5c857ca12cf8da5c26",
       "0.8": "4d83728d4925efb313a2b69f7fe74d64f49007284e5286beea99f3912c100426", "0.75": "0d45a06cb93b656289c2281b599252cc9771f93cbce576228f48f4400c039a41",
       "0.7": "a52848b5320cafc7fcd7316ef3f94e98ebbe0aa75a7a0040215e4f1d0800e247", "0.6": "f2566bff023b75b0b20ce7ef8096ac7abcb3e3537ff04a4d3c94653f527def4a"}
ns = {}
exec(src[src.index("GENES = ("):src.index("EXTRA_NAMES")], {}, ns)
V15 = {v[0]: v for v in ns["VARIANTS"]}
G1 = ["KI_18w_mid", "KI_1c_restore_Cacna1a", "KI_4_protein_Kcnma1_Cacna1g", "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85"]
G2 = ["KI_3_only_Kcnc3", "KI_3_only_Kcnma1", "KI_3_only_Kcnd3", "KI_3_only_Cacna1a", "KI_3_only_Cacna1g",
      "KI_1c_restore_Kcnc3", "KI_1c_restore_Kcnma1", "KI_1c_restore_Kcnd3", "KI_1c_restore_Cacna1g",
      "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3", "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3"]
G3 = [v for v in V15 if (v.startswith(("KI_3_only_", "KI_2b_", "KI_4_", "KI_5_")) and v not in G1 + G2)]
assert len(G1) == 4 and len(G2) == 11 and len(G3) == 24, (len(G1), len(G2), len(G3))
assert len(set(G1 + G2 + G3)) == 39
base = V15["KI_18w_mid"][3]
area = [l.rstrip("\n").split("\t") for l in open(sys.argv[1]) if l.strip()]
if any(len(a) != 3 for a in area) or [a[0] for a in area] != ["0.9", "0.85", "0.8", "0.75", "0.7", "0.6"]:
    sys.exit("STOP: AREA_TSV must hold six lines 'A<TAB>sha256<TAB>folder', A = 0.9, 0.85, 0.8, 0.75, 0.7, 0.6 in this order")
for A, sha, root in area:
    if sha != PIN[A]:
        sys.exit("STOP: the hash given for A %s is not that of its wave-9 factor file" % A)
    if not re.fullmatch(r"[A-Za-z0-9_.]+", root):
        sys.exit("STOP: folder name of A %s has a character outside letters, digits, underscore and point" % A)
taken = ["wave8_syn_" + v for v in G1 + G2 + G3] + ["wave8_syn_WT"]
names = [a[2] for a in area]
if len(set(names)) != 6 or any(n in taken or n.endswith("_copy") for n in names):
    sys.exit("STOP: the six folder names must differ from each other and from the folders of the other versions, and must not end in _copy")
rows = []
def row(name, sha, root, factors, group):
    return '    ("%s", "%s", "%s", (%s), "%s"),' % (name, sha, root, ", ".join(repr(float(x)) for x in factors), group)
for v in G1:
    rows.append(row(v, V15[v][1], "wave8_syn_" + v, V15[v][3], "1"))
for A, sha, root in area:
    rows.append(row("P8_atrophy_A" + A, sha, root, [x / float(A) for x in base], "1"))
for v in G2:
    rows.append(row(v, V15[v][1], "wave8_syn_" + v, V15[v][3], "2"))
for v in G3:
    rows.append(row(v, V15[v][1], "wave8_syn_" + v, V15[v][3], "3"))
out = src
def rep(old, new, count=1):
    global out
    assert out.count(old) == count, (old[:70], out.count(old))
    out = out.replace(old, new)
a, b = out.index("VARIANTS = ["), out.index("EXTRA_NAMES")
out = out[:a] + "VARIANTS = [\n" + "\n".join(rows) + "\n]\n" + out[b:]
rep('EXTRA_NAMES = {"KI_18w_mid": ["ki_copy_20261113"]}', "EXTRA_NAMES = {}")
rep('SCORER_SHA = "cbf5189d549b25c689dc7f390b0a05511f0a125cd742a6a902d3819eae632dd7"', 'SCORER_SHA = "%s"' % scorer_sha)
rep('"""P15: score every finished variant root and print one line per variant (14 Nov 2026; v02 = v01 with the six\nwave-5 variants added to the table and a column aux).\n',
    '"""P16 (synaptic-input runs, wave 8): score every finished variant root and print one line per variant (14 Nov 2026).\n'
    'This is p15_score_all_v05.py (bfa035d6...) with the list of the wave-8 versions, 131 sets per root, the scoring script\n'
    'p16_variant_score_v02.py, and the last entry of each list line used as the submission group (column "group", in the\n'
    'place of the column aux). Everything else is unchanged.\n')
rep("Usage: python3 -I -B p15_score_all_v05.py SCORER_PY WT_ROOT COPY_PARENT OUT_DIR TABLE_OUT", "Usage: python3 -I -B p16_score_all_v01.py SCORER_PY WT_ROOT COPY_PARENT OUT_DIR TABLE_OUT")
rep("SCORER_PY    p15_variant_score_v01.py (its sha256 is checked against the value below)", "SCORER_PY    p16_variant_score_v02.py (its sha256 is checked against the value below)")
rep("             table below); KI_18w_mid is also looked for as ki_copy_20261113\n", "             table below)\n")
rep("the copy is used only if it holds 130 set folders that all", "the copy is used only if it holds 131 set folders that all")
rep("    if len(sets) != 130 or len(done) != 130:", "    if len(sets) != 131 or len(done) != 131:")
rep('    cols = ["variant"] + list(GENES) + ["aux", "status",', '    cols = ["variant"] + list(GENES) + ["group", "status",')
rep("level); aux (b4, a2d2, b4+a2d2 or -); status;", "level); group (submission group 1, 2 or 3); status;")
rep("# variant, sha256 of its factor file, name of its result folder on the cluster, factors of the five genes, auxiliary-subunit\n# changes (b4 = beta4_shift_r 0.38; a2d2 = newCaP.pcabar 0.85; - = none)\n",
    "# variant, sha256 of its factor file, name of its result folder on the cluster, factors of the five genes, submission group\n")
open(sys.argv[2], "x").write(out)
print("variants", len(rows), "scorer", scorer_sha[:12], "written", sys.argv[2], hashlib.sha256(out.encode()).hexdigest())
