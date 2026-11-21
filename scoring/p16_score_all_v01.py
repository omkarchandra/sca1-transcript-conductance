"""P16 (synaptic-input runs, wave 8): score every finished variant root and print one line per variant (14 Nov 2026).
This is p15_score_all_v05.py (bfa035d6...) with the list of the wave-8 versions, 131 sets per root, the scoring script
p16_variant_score_v02.py, and the last entry of each list line used as the submission group (column "group", in the
place of the column aux). Everything else is unchanged.

v05 (14 Nov 2026) = v04 with two appended columns, n_rate_recorded and n_all_three_recorded: as n_rate_0.85-1.33 and
n_all_three, but with the recorded firing-rate-ratio range of Bushart et al. 2021, Supplementary Table 1 (bounds computed below
from the four firing frequencies), counted from the per-set table <variant>.tsv (columns rate_ratio, cv_ratio, d_trough_mv).

Usage: python3 -I -B p16_score_all_v01.py SCORER_PY WT_ROOT COPY_PARENT OUT_DIR TABLE_OUT

SCORER_PY    p16_variant_score_v02.py (its sha256 is checked against the value below)
WT_ROOT      the copy of the WT root (folders set<NNN>_WT)
COPY_PARENT  the folder that holds the copies of the variant roots; a variant is looked for under the names
             <root>_copy and <root>, where <root> is the name of its result folder on the cluster (third entry of the
             table below)
OUT_DIR      an existing folder for the per-variant files <variant>.tsv and <variant>.summary.txt
TABLE_OUT    the one-line-per-variant table; must not exist

For every variant of the table, in the table's order: the copy is used only if it holds 131 set folders that all
have SHA256SUMS. A variant that already has both files in OUT_DIR is not scored again. Otherwise the scoring script is run once, as
"python3 -I -B SCORER_PY WT_ROOT <copy> <factor sha256> OUT_DIR/<variant>.tsv"; if it ends without error its printed
summary is saved as OUT_DIR/<variant>.summary.txt. If it stops, its message is shown in the status column and the
script goes on with the next variant. Every file is written by exclusive create; nothing is changed or removed.

The entries of the table are read from the saved summaries; nothing is computed again. A summary is used only if
it names the factor file of the variant and the sha256 of the .tsv beside it. Columns (KI relative to the WT run of
the same parameter set; sets with 3 or more window spikes in both runs): the five factors of the variant (1 = WT
level); group (submission group 1, 2 or 3); status; sets; both_fire; rate ratio quartiles; sets with rate ratio 0.85-1.33; median CV ratio, over all
sets and over the sets without pauses; sets with CV ratio 1.5 or more, 1.8 or more; median trough difference (mV);
sets with the KI trough above WT; sets with all three (rate ratio 0.85-1.33, CV ratio 1.5 or more, trough above WT);
KI states; first 12 characters of the sha256 of the per-set table.
"""
import hashlib, os, re, subprocess, sys

SCORER_SHA = "1e0a200d7d6a6e066e6d48796b84fec47f527a5b37a4185274a950fd184d9a5e"
GENES = ("Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g")
# variant, sha256 of its factor file, name of its result folder on the cluster, factors of the five genes, submission group
VARIANTS = [
    ("KI_18w_mid", "660f8462aec5b6671a279c5041ee539909b8a8395de3c7eed30994ca283c937a", "wave8_syn_KI_18w_mid", (0.6666093556941, 0.7190924725561185, 0.5876215270426435, 0.6722697088324129, 0.6491953709552294), "1"),
    ("KI_1c_restore_Cacna1a", "5b11a8fcbf75fb4999bf739a22096b5216f50f240b78efc59af8b4fde6b51fbe", "wave8_syn_KI_1c_restore_Cacna1a", (0.6666093556941, 0.7190924725561185, 0.5876215270426435, 1.0, 0.6491953709552294), "1"),
    ("KI_4_protein_Kcnma1_Cacna1g", "367ed13bdaa0f07fa2067a76df2c53e96df998ed4c8e44f42da0b49187555b22", "wave8_syn_KI_4_protein_Kcnma1_Cacna1g", (1.0, 0.78, 1.0, 1.0, 0.75), "1"),
    ("KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85", "a657ad914d32cfb3d1705def32f16eb3abe817750200958f1f3a24a8db15268f", "wave8_syn_KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p85", (1.0, 0.78, 1.0, 0.85, 0.75), "1"),
    ("P8_atrophy_A0.9", "25264b8faa433761c2468e04a2947351c3a438d50c7d3670f8801eb1b8d9f5c7", "wave8_syn_KI_9_area_A0p9", (0.7406770618823333, 0.798991636173465, 0.6529128078251594, 0.7469663431471254, 0.7213281899502548), "1"),
    ("P8_atrophy_A0.85", "7253210c3c233de005824a6374235c6052925c200e045f5c857ca12cf8da5c26", "wave8_syn_KI_9_area_A0p85", (0.7842463008165882, 0.8459911441836688, 0.6913194435795805, 0.7909055398028387, 0.7637592599473287), "1"),
    ("P8_atrophy_A0.8", "4d83728d4925efb313a2b69f7fe74d64f49007284e5286beea99f3912c100426", "wave8_syn_KI_9_area_A0p8", (0.8332616946176249, 0.8988655906951482, 0.7345269088033043, 0.8403371360405161, 0.8114942136940367), "1"),
    ("P8_atrophy_A0.75", "0d45a06cb93b656289c2281b599252cc9771f93cbce576228f48f4400c039a41", "wave8_syn_KI_9_area_A0p75", (0.8888124742588, 0.958789963408158, 0.7834953693901913, 0.8963596117765505, 0.8655938279403058), "1"),
    ("P8_atrophy_A0.7", "a52848b5320cafc7fcd7316ef3f94e98ebbe0aa75a7a0040215e4f1d0800e247", "wave8_syn_KI_9_area_A0p7", (0.952299079563, 1.0272749607944551, 0.8394593243466336, 0.9603852983320185, 0.9274219585074707), "1"),
    ("P8_atrophy_A0.6", "f2566bff023b75b0b20ce7ef8096ac7abcb3e3537ff04a4d3c94653f527def4a", "wave8_syn_KI_9_area_A0p6", (1.1110155928235, 1.1984874542601975, 0.9793692117377392, 1.120449514720688, 1.0819922849253825), "1"),
    ("KI_3_only_Kcnc3", "120bd84515cc864a8efa0c3dd0d97ba37646557bbf1dc5d83cd01166605a9982", "wave8_syn_KI_3_only_Kcnc3", (0.6666093556941, 1.0, 1.0, 1.0, 1.0), "2"),
    ("KI_3_only_Kcnma1", "10610f9327d94f7b5563f743aee0980c6a8b9e13741fdc836f7ffff44f27a4e9", "wave8_syn_KI_3_only_Kcnma1", (1.0, 0.7190924725561185, 1.0, 1.0, 1.0), "2"),
    ("KI_3_only_Kcnd3", "04b388690f34cd9ad2520e9ddcddc44c853cb1e4bc1a7984561e49d126ec9581", "wave8_syn_KI_3_only_Kcnd3", (1.0, 1.0, 0.5876215270426435, 1.0, 1.0), "2"),
    ("KI_3_only_Cacna1a", "541a8b940179df4c80cf1ed17d5bc25e3abc74ff4e10a110df46055d8c9c83ea", "wave8_syn_KI_3_only_Cacna1a", (1.0, 1.0, 1.0, 0.6722697088324129, 1.0), "2"),
    ("KI_3_only_Cacna1g", "65f5be13402bb0077714805dc8c90695ef80d4fb8b0a15bb8dd9cb9a526ccb20", "wave8_syn_KI_3_only_Cacna1g", (1.0, 1.0, 1.0, 1.0, 0.6491953709552294), "2"),
    ("KI_1c_restore_Kcnc3", "df44d4080056959e957a431436c19fc5f7978538259d6fbdc4a987479ca6d7f0", "wave8_syn_KI_1c_restore_Kcnc3", (1.0, 0.7190924725561185, 0.5876215270426435, 0.6722697088324129, 0.6491953709552294), "2"),
    ("KI_1c_restore_Kcnma1", "f1d8bcb702484db8cfd3087407c467a91e284e3e10c42601548c4513b4e733e7", "wave8_syn_KI_1c_restore_Kcnma1", (0.6666093556941, 1.0, 0.5876215270426435, 0.6722697088324129, 0.6491953709552294), "2"),
    ("KI_1c_restore_Kcnd3", "d54bd6e4d307beecd363e16e065c83151702d64dec3e91ceae7fd08f92c881bf", "wave8_syn_KI_1c_restore_Kcnd3", (0.6666093556941, 0.7190924725561185, 1.0, 0.6722697088324129, 0.6491953709552294), "2"),
    ("KI_1c_restore_Cacna1g", "1b6b5bd5f3ff27a8ca557033b3568e488e649ac63a18f5213eaec5d75eb9fc1f", "wave8_syn_KI_1c_restore_Cacna1g", (0.6666093556941, 0.7190924725561185, 0.5876215270426435, 0.6722697088324129, 1.0), "2"),
    ("KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3", "7227127e499155e4665f212896cbb18c886aa26da5bfa690e8de4bd494e92b2e", "wave8_syn_KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3", (0.6666093556941, 0.78, 1.0, 1.0, 0.75), "2"),
    ("KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3", "0fe043257c838d8df8adf3b12081fad2ccd8c85786d0c3828bd6905bda31bf1d", "wave8_syn_KI_4_protein_Kcnma1_Cacna1g_plus_Kcnd3", (1.0, 0.78, 0.5876215270426435, 1.0, 0.75), "2"),
    ("KI_3_only_Kcnc3_Kcnma1", "b9f105ce7cf168c75ea197c77550e93e85dc33d8635095c5aee6ced6d6834cee", "wave8_syn_KI_3_only_Kcnc3_Kcnma1", (0.6666093556941, 0.7190924725561185, 1.0, 1.0, 1.0), "3"),
    ("KI_3_only_Kcnc3_Kcnd3", "34083ac4e53a384b54a7629fb89c5a9af5effbde840a30770b52f702f439aa6e", "wave8_syn_KI_3_only_Kcnc3_Kcnd3", (0.6666093556941, 1.0, 0.5876215270426435, 1.0, 1.0), "3"),
    ("KI_3_only_Kcnc3_Cacna1g", "f058363e93a0bd0ccda509746210fa3d699604197807a9e5d9567ad1626763e0", "wave8_syn_KI_3_only_Kcnc3_Cacna1g", (0.6666093556941, 1.0, 1.0, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnma1_Kcnd3", "993b670e1fc5c4b145f83744ec40d35a89804aa65bab767fced5d4aab4560648", "wave8_syn_KI_3_only_Kcnma1_Kcnd3", (1.0, 0.7190924725561185, 0.5876215270426435, 1.0, 1.0), "3"),
    ("KI_2b_only_Kcnma1_Cacna1g", "aa182bf2b36f35b758ae9235baa62ef89f560506db69c16d35de725c8ac322ee", "wave8_syn_KI_2b_only_Kcnma1_Cacna1g", (1.0, 0.7190924725561185, 1.0, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnd3_Cacna1g", "70130bc2a12bbf1fb43b5fb13a36394e80b1a65690a9612c00b99bbab4f5c2bf", "wave8_syn_KI_3_only_Kcnd3_Cacna1g", (1.0, 1.0, 0.5876215270426435, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnc3_Kcnma1_Kcnd3", "caff5fb8d0a17ecfdafee8c4e79c3c141431546b548eda1cd71830447f9ae3ba", "wave8_syn_KI_3_only_Kcnc3_Kcnma1_Kcnd3", (0.6666093556941, 0.7190924725561185, 0.5876215270426435, 1.0, 1.0), "3"),
    ("KI_3_only_Kcnc3_Kcnma1_Cacna1g", "89453f3af78955f8c2d9aa48761748d82098b5dd75596c1c846bc7e6c9b0cdaf", "wave8_syn_KI_3_only_Kcnc3_Kcnma1_Cacna1g", (0.6666093556941, 0.7190924725561185, 1.0, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnc3_Kcnd3_Cacna1g", "676f6ba8381a37a8263bc7d12f084e7d42059a73e986a698ac6a760b9de344cc", "wave8_syn_KI_3_only_Kcnc3_Kcnd3_Cacna1g", (0.6666093556941, 1.0, 0.5876215270426435, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnma1_Kcnd3_Cacna1g", "35885d4c2124156a830676a005e0230070a4cc9cf124219e90a863661a3224ad", "wave8_syn_KI_3_only_Kcnma1_Kcnd3_Cacna1g", (1.0, 0.7190924725561185, 0.5876215270426435, 1.0, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnc3_Cacna1a", "af845509b03ace32b3476e4ceeef7b650fc615a149712c2d55faae44838a9f23", "wave8_syn_KI_3_only_Kcnc3_Cacna1a", (0.6666093556941, 1.0, 1.0, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Kcnma1_Cacna1a", "13874e5eab4530cbcbbb1d8e1aae9fac696c5eb1720a6e6b4115a63e9ebdf5f1", "wave8_syn_KI_3_only_Kcnma1_Cacna1a", (1.0, 0.7190924725561185, 1.0, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Kcnd3_Cacna1a", "9563472a18a9a394e2b76e1c813c8bd9da72198746d5f118c75bf2a37de34ae4", "wave8_syn_KI_3_only_Kcnd3_Cacna1a", (1.0, 1.0, 0.5876215270426435, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Cacna1a_Cacna1g", "fef91b4a6ffbd7d824600721ac0573364c107ff3429848334265af1d10db5708", "wave8_syn_KI_3_only_Cacna1a_Cacna1g", (1.0, 1.0, 1.0, 0.6722697088324129, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnc3_Kcnma1_Cacna1a", "0339b5e53d54557b1da973b48a2898f09e93501815a56d1fdccf6e8bb1dddb04", "wave8_syn_KI_3_only_Kcnc3_Kcnma1_Cacna1a", (0.6666093556941, 0.7190924725561185, 1.0, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Kcnc3_Kcnd3_Cacna1a", "ff3cef75f7af83eb9032224d13301e7b0afae539f50992def38d1eee2c3cb6c5", "wave8_syn_KI_3_only_Kcnc3_Kcnd3_Cacna1a", (0.6666093556941, 1.0, 0.5876215270426435, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Kcnc3_Cacna1a_Cacna1g", "bfe8ff885528b299f78a96d3ba9564e408573c8214943618bff3ed35971bafdc", "wave8_syn_KI_3_only_Kcnc3_Cacna1a_Cacna1g", (0.6666093556941, 1.0, 1.0, 0.6722697088324129, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnma1_Kcnd3_Cacna1a", "735afd2bd39f9e35668af2bfe16b21f65d8fb1ff66fbca3eec0b485ac7e7e503", "wave8_syn_KI_3_only_Kcnma1_Kcnd3_Cacna1a", (1.0, 0.7190924725561185, 0.5876215270426435, 0.6722697088324129, 1.0), "3"),
    ("KI_3_only_Kcnma1_Cacna1a_Cacna1g", "99eeed1006ef972048e87268c5955450401b4ebe109ab1911a510e780ec33a0f", "wave8_syn_KI_3_only_Kcnma1_Cacna1a_Cacna1g", (1.0, 0.7190924725561185, 1.0, 0.6722697088324129, 0.6491953709552294), "3"),
    ("KI_3_only_Kcnd3_Cacna1a_Cacna1g", "6702c6624ce7b3a9d60e28a42c29760ce022db6effb9988eb88e9a994c234acd", "wave8_syn_KI_3_only_Kcnd3_Cacna1a_Cacna1g", (1.0, 1.0, 0.5876215270426435, 0.6722697088324129, 0.6491953709552294), "3"),
    ("KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3_Kcnd3", "1940e87b678313438b3fa6db72239c8297698169a3e6229426ae544ca55fe4e6", "wave8_syn_KI_4_protein_Kcnma1_Cacna1g_plus_Kcnc3_Kcnd3", (0.6666093556941, 0.78, 0.5876215270426435, 1.0, 0.75), "3"),
    ("KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p95", "106e921da57e838ee8dbeda719dae66c832e0d1eba256ac63314a6756a66da9f", "wave8_syn_KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p95", (1.0, 0.78, 1.0, 0.95, 0.75), "3"),
    ("KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p90", "34489efdf2075a044b8912ce7911353f4b8a118d0b688dc2cbbf73af344f7fb1", "wave8_syn_KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p90", (1.0, 0.78, 1.0, 0.9, 0.75), "3"),
    ("KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p80", "7d5b1070f759ce90d9660c0c312ce62a5cc0be97beafd960c227025178a46d92", "wave8_syn_KI_5_protein_Kcnma1_Cacna1g_Cacna1a_0p80", (1.0, 0.78, 1.0, 0.8, 0.75), "3"),
]
EXTRA_NAMES = {}
# Bushart et al. 2021, Supplementary Table 1, firing frequency (Hz): anterior WT, anterior KI, posterior WT, posterior KI
ANT_WT, ANT_KI, POST_WT, POST_KI = 31.65, 42.04, 27.81, 23.78
REC_LO, REC_HI = POST_KI / POST_WT, ANT_KI / ANT_WT


def sha_of(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def find_copy(parent, variant, root):
    names = [root + "_copy", root] + EXTRA_NAMES.get(variant, [])
    found = [n for n in names if os.path.isdir(os.path.join(parent, n))]
    if not found:
        return None, "not copied"
    if len(found) > 1:
        return None, "more than one copy: " + " ".join(found)
    path = os.path.join(parent, found[0])
    sets = [n for n in os.listdir(path) if re.fullmatch(r"set\d+_KI_18w_mid", n)]
    done = [n for n in sets if os.path.isfile(os.path.join(path, n, "SHA256SUMS"))]
    if len(sets) != 131 or len(done) != 131:
        return None, "incomplete: %d set folders, %d finished" % (len(sets), len(done))
    return path, "complete"


PATTERNS = [
    ("ident", r"^condition (\S+); factor file sha256 ([0-9a-f]{64}); sets (\d+);"),
    ("states", r"^KI states: fires (\d+), pauses (\d+), block (\d+), rest (\d+)$"),
    ("both", r"^sets with 3 or more spikes in both runs: (\d+);"),
    ("rate", r"^firing rate ratio KI/WT: n \d+; min \S+, quartiles (\S+) / (\S+) / (\S+), max "),
    ("cv", r"^CV ratio KI/WT: n \d+; min \S+, quartiles \S+ / (\S+) / \S+, max "),
    ("cvbands", r"^  CV ratio below 1\.0: (\d+); 1\.0 to below 1\.5: (\d+); 1\.5 to below 1\.8: (\d+); 1\.8 to below 2\.8: (\d+); 2\.8 or more: (\d+)$"),
    ("cvnop", r"^CV ratio KI/WT, sets without pauses: n \d+; min \S+, quartiles \S+ / (\S+) / \S+, max "),
    ("trough", r"^trough difference KI - WT \(mV\): n \d+; min \S+, quartiles \S+ / (\S+) / \S+, max "),
    ("above", r"^  KI trough above WT: (\d+); equal: "),
    ("band", r"^sets with rate ratio 0\.85-1\.33: (\d+); of those with CV ratio 1\.5 or more: \d+; with trough difference above 0: \d+; with both: (\d+)$"),
    ("out", r"^OUT_TSV .*: (\d+) rows, sha256 ([0-9a-f]{64})$"),
]


def measures(text, fsha, tsv_sha):
    """The table entries of one variant, taken from the saved summary of the scoring script (no recomputation)."""
    g = {}
    for key, pat in PATTERNS:
        hits = re.findall(pat, text, flags=re.M)
        if len(hits) != 1:
            return None, "STOP: summary line not found once: " + key
        g[key] = hits[0] if isinstance(hits[0], tuple) else (hits[0],)
    if g["ident"][1] != fsha:
        return None, "STOP: the summary is for another factor file"
    if g["out"][1] != tsv_sha or g["out"][0] != g["ident"][2]:
        return None, "STOP: summary and table in OUT_DIR do not agree"
    c = [int(x) for x in g["cvbands"]]
    return [g["ident"][2], g["both"][0], g["rate"][0], g["rate"][1], g["rate"][2], g["band"][0], g["cv"][0], g["cvnop"][0],
            "%d" % (c[2] + c[3] + c[4]), "%d" % (c[3] + c[4]), g["trough"][0], g["above"][0], g["band"][1],
            g["states"][0], g["states"][1], g["states"][2], g["states"][3]], None


def recorded_counts(tsv):
    """Sets with a rate ratio in [REC_LO, REC_HI]; of those, sets also with CV ratio 1.5 or more and trough difference above 0."""
    n_rate = n_all = 0
    with open(tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_rr, i_cv, i_dt = head.index("rate_ratio"), head.index("cv_ratio"), head.index("d_trough_mv")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rr, cv, dt = float(f[i_rr]), float(f[i_cv]), float(f[i_dt])
            if REC_LO <= rr <= REC_HI:
                n_rate += 1
                n_all += cv >= 1.5 and dt > 0
    return ["%d" % n_rate, "%d" % n_all]


def main(scorer, wt_root, parent, out_dir, table_out):
    if sha_of(scorer) != SCORER_SHA:
        sys.exit("STOP: SCORER_PY does not have the expected sha256")
    for d in (wt_root, parent, out_dir):
        if not os.path.isdir(d):
            sys.exit("STOP: not a folder: " + d)
    if os.path.exists(table_out):
        sys.exit("STOP: TABLE_OUT exists: " + table_out)
    cols = ["variant"] + list(GENES) + ["group", "status", "sets", "both_fire", "rate_q25", "rate_median", "rate_q75",
            "n_rate_0.85-1.33", "cv_ratio_median", "cv_ratio_median_no_pause", "n_cv_ge_1.5", "n_cv_ge_1.8",
            "d_trough_median_mv", "n_trough_above", "n_all_three", "ki_fires", "ki_pauses", "ki_block", "ki_rest",
            "tsv_sha256_12", "n_rate_recorded", "n_all_three_recorded"]
    lines, new, old = ["\t".join(cols)], [], []
    for variant, fsha, root, factors, aux in VARIANTS:
        head = [variant] + [("1" if f == 1.0 else "%.3f" % f) for f in factors] + [aux]
        tsv = os.path.join(out_dir, variant + ".tsv")
        summ = os.path.join(out_dir, variant + ".summary.txt")
        have_t, have_s = os.path.isfile(tsv), os.path.isfile(summ)
        if have_t and have_s:
            status = "scored earlier"
        elif have_t or have_s:
            lines.append("\t".join(head + ["STOP: only one of the two files is in OUT_DIR"]))
            continue
        else:
            path, state = find_copy(parent, variant, root)
            if path is None:
                lines.append("\t".join(head + [state]))
                continue
            p = subprocess.run([sys.executable, "-I", "-B", scorer, wt_root, path, fsha, tsv],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if p.returncode != 0 or not os.path.isfile(tsv):
                msg = (p.stderr.strip().split("\n") or [""])[-1]
                lines.append("\t".join(head + ["scoring script stopped: " + msg[:300].replace("\t", " ")]))
                continue
            with open(summ, "x") as fh:
                fh.write(p.stdout)
            status = "scored now"
        with open(summ) as fh:
            vals, problem = measures(fh.read(), fsha, sha_of(tsv))
        if problem:
            lines.append("\t".join(head + [problem]))
            continue
        (new if status == "scored now" else old).append(variant)
        lines.append("\t".join(head + [status] + vals + [sha_of(tsv)[:12]] + recorded_counts(tsv)))
    data = ("\n".join(lines) + "\n").encode("ascii", "replace")
    with open(table_out, "xb") as fh:
        fh.write(data)
    sys.stdout.write(data.decode("ascii"))
    print("variants in the table: %d; scored now: %d; scored earlier: %d; without a result: %d"
          % (len(VARIANTS), len(new), len(old), len(VARIANTS) - len(new) - len(old)))
    print("TABLE_OUT %s: sha256 %s" % (table_out, hashlib.sha256(data).hexdigest()))
    print("recorded rate-ratio range used for n_rate_recorded: %.6f (%s/%s) to %.6f (%s/%s), both ends included"
          % (REC_LO, POST_KI, POST_WT, REC_HI, ANT_KI, ANT_WT))


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
