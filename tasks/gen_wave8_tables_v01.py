"""P16 wave-8 task tables. Run on the cluster: python3 -I -B - P15 < this file. Writes by exclusive create into
the new folder P15/inputs/wave8: WAVE8_G1_TASKS_20261114_v01.tsv and WAVE8_G2_TASKS_20261114_v01.tsv (11 versions x 131
sets each; version-major, so task 0 = set 000 WT and task 131 = set 000 KI_18w_mid in group 1). Sets: 000 (all-1.0 file
popsets_wave8/popset_k000.json) then the 130 kept W2 sets in the order of receipts/WAVE2_KI_TASKS.tsv. Seeds: set 000 -> 20261110,
W2 set k -> 100000000 + 1000000 x k. Every factor file is re-hashed against its P15 pin; every W2 set file against popsets_wave2/SHA256SUMS.
Output folders: results/wave8_syn_<version>/set<NNN>_WT (WT) or set<NNN>_KI_18w_mid (every KI version)."""
import hashlib, os, sys

P15 = sys.argv[1]
P13 = "${PROJECT_ROOT}/p13_zang_pop_20261109"
F = {  # version -> (factor file, P15 pin)
 "WT": (P13 + "/inputs/factors/WT.json", "bdbb0f458069e46aefc13f6dfe0dd488de397fcf1546821dca2c3b804157c204"),
 "KI_18w_mid": (P13 + "/inputs/factors/KI_18w_mid.json", "660f8462aec5b6671a279c5041ee539909b8a8395de3c7eed30994ca283c937a"),
 "KI_1c_restore_Cacna1a": (P15 + "/inputs/factors_wave1_1c1d/KI_1c_restore_Cacna1a.json", "5b11a8fcbf75fb4999bf739a22096b5216f50f240b78efc59af8b4fde6b51fbe"),
 "KI_4_protein_Kcnma1_Cacna1g": (P15 + "/inputs/factors_wave3/KI_4_protein_Kcnma1_Cacna1g.json", "367ed13bdaa0f07fa2067a76df2c53e96df998ed4c8e44f42da0b49187555b22"),
 "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85": (P15 + "/inputs/factors_wave4/KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85.json", "a657ad914d32cfb3d1705def32f16eb3abe817750200958f1f3a24a8db15268f"),
 "KI_9_area_A0p8": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.8.json", "4d83728d4925efb313a2b69f7fe74d64f49007284e5286beea99f3912c100426"),
 "KI_9_area_A0p75": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.75.json", "0d45a06cb93b656289c2281b599252cc9771f93cbce576228f48f4400c039a41"),
 "KI_9_area_A0p9": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.9.json", "25264b8faa433761c2468e04a2947351c3a438d50c7d3670f8801eb1b8d9f5c7"),
 "KI_9_area_A0p85": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.85.json", "7253210c3c233de005824a6374235c6052925c200e045f5c857ca12cf8da5c26"),
 "KI_9_area_A0p7": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.7.json", "a52848b5320cafc7fcd7316ef3f94e98ebbe0aa75a7a0040215e4f1d0800e247"),
 "KI_9_area_A0p6": (P15 + "/inputs/factors_wave9/P8_atrophy_A0.6.json", "f2566bff023b75b0b20ce7ef8096ac7abcb3e3537ff04a4d3c94653f527def4a"),
 "KI_3_only_Kcnc3": (P15 + "/inputs/factors_wave3/KI_3_only_Kcnc3.json", "120bd84515cc864a8efa0c3dd0d97ba37646557bbf1dc5d83cd01166605a9982"),
 "KI_3_only_Kcnma1": (P15 + "/inputs/factors_wave3/KI_3_only_Kcnma1.json", "10610f9327d94f7b5563f743aee0980c6a8b9e13741fdc836f7ffff44f27a4e9"),
 "KI_3_only_Kcnd3": (P15 + "/inputs/factors_wave3/KI_3_only_Kcnd3.json", "04b388690f34cd9ad2520e9ddcddc44c853cb1e4bc1a7984561e49d126ec9581"),
 "KI_3_only_Cacna1a": (P15 + "/inputs/factors_wave3/KI_3_only_Cacna1a.json", "541a8b940179df4c80cf1ed17d5bc25e3abc74ff4e10a110df46055d8c9c83ea"),
 "KI_3_only_Cacna1g": (P15 + "/inputs/factors_wave3/KI_3_only_Cacna1g.json", "65f5be13402bb0077714805dc8c90695ef80d4fb8b0a15bb8dd9cb9a526ccb20"),
 "KI_1c_restore_Kcnc3": (P15 + "/inputs/factors_wave1_1c1d/KI_1c_restore_Kcnc3.json", "df44d4080056959e957a431436c19fc5f7978538259d6fbdc4a987479ca6d7f0"),
 "KI_1c_restore_Kcnma1": (P15 + "/inputs/factors_wave1_1c1d/KI_1c_restore_Kcnma1.json", "f1d8bcb702484db8cfd3087407c467a91e284e3e10c42601548c4513b4e733e7"),
 "KI_1c_restore_Kcnd3": (P15 + "/inputs/factors_wave1_1c1d/KI_1c_restore_Kcnd3.json", "d54bd6e4d307beecd363e16e065c83151702d64dec3e91ceae7fd08f92c881bf"),
 "KI_1c_restore_Cacna1g": (P15 + "/inputs/factors_wave1_1c1d/KI_1c_restore_Cacna1g.json", "1b6b5bd5f3ff27a8ca557033b3568e488e649ac63a18f5213eaec5d75eb9fc1f"),
 "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3": (P15 + "/inputs/factors_wave3/KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3.json", "7227127e499155e4665f212896cbb18c886aa26da5bfa690e8de4bd494e92b2e"),
 "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3": (P15 + "/inputs/factors_wave3/KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3.json", "0fe043257c838d8df8adf3b12081fad2ccd8c85786d0c3828bd6905bda31bf1d"),
}
G1 = ["WT", "KI_18w_mid", "KI_1c_restore_Cacna1a", "KI_4_protein_Kcnma1_Cacna1g", "KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85",
      "KI_9_area_A0p8", "KI_9_area_A0p75", "KI_9_area_A0p9", "KI_9_area_A0p85", "KI_9_area_A0p7", "KI_9_area_A0p6"]
G2 = ["KI_3_only_Kcnc3", "KI_3_only_Kcnma1", "KI_3_only_Kcnd3", "KI_3_only_Cacna1a", "KI_3_only_Cacna1g", "KI_1c_restore_Kcnc3",
      "KI_1c_restore_Kcnma1", "KI_1c_restore_Kcnd3", "KI_1c_restore_Cacna1g", "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3",
      "KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3"]
SET0 = (P15 + "/inputs/popsets_wave8/popset_k000.json", "70d8938fb44715159d90b4a76d515d9810efdb4cb38183036a78156c467fbad3")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


for v, (p, pin) in list(F.items()) + [("set000", SET0)]:
    if sha(p) != pin:
        sys.exit("STOP: hash of %s (%s)" % (v, p))
TK = P15 + "/receipts/WAVE2_KI_TASKS.tsv"
if sha(TK) != open(P15 + "/receipts/WAVE2_KI_TASKS.sha256").read().split()[0]:
    sys.exit("STOP: WAVE2_KI_TASKS hash")
kept = [int(l.split("\t")[1]) for l in open(TK).read().splitlines()[1:]]
psum = {}
for l in open(P15 + "/inputs/popsets_wave2/SHA256SUMS"):
    h, n = l.split()
    psum[n] = h
if len(kept) != 130:
    sys.exit("STOP: %d kept sets" % len(kept))
sets = [(0, SET0[0], SET0[1])]
for k in kept:
    n = "popset_k%03d.json" % k
    p = P15 + "/inputs/popsets_wave2/" + n
    if sha(p) != psum[n]:
        sys.exit("STOP: set file " + n)
    sets.append((k, p, psum[n]))
os.mkdir(P15 + "/inputs/wave8")
for name, group in (("WAVE8_G1_TASKS_20261114_v01.tsv", G1), ("WAVE8_G2_TASKS_20261114_v01.tsv", G2)):
    rows = ["task\tversion\tset_id\tsyn_seed\tfactor_file\tfactor_sha256\tpopset_file\tpopset_sha256\tout_root\tout_folder"]
    for v in group:
        for k, sp, ssh in sets:
            seed = 20261110 if k == 0 else 100000000 + 1000000 * k
            folder = "set%03d_%s" % (k, "WT" if v == "WT" else "KI_18w_mid")
            rows.append("\t".join([str(len(rows) - 1), v, str(k), str(seed), F[v][0], F[v][1], sp, ssh, P15 + "/results/wave8_syn_" + v, folder]))
    data = ("\n".join(rows) + "\n").encode("ascii")
    open(P15 + "/inputs/wave8/" + name, "xb").write(data)
    print("%s: %d tasks (%d versions x %d sets); sha256 %s" % (name, len(rows) - 1, len(group), len(sets), hashlib.sha256(data).hexdigest()))
print("set000 task index in G1: WT 0, KI_18w_mid %d" % len(sets))

#%% debug
