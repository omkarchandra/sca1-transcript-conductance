# P7 transcriptional output per Purkinje nucleus

P7: per-mouse Purkinje / granule / all-other median UMI and detected genes from the h5ad raw counts; KI/WT per age.

| age | KI/WT R_gran (UMI) | KI per-mouse R_gran | WT per-mouse R_gran | all KI below lowest WT | KI/WT R_all (UMI) | KI/WT R_gran (genes) | KI/WT R_all (genes) | KI/WT Purkinje share | reading |
|---|---|---|---|---|---|---|---|---|---|
| 5 w | 0.962 | 7.652, 7.000, 7.344 | 7.679, 7.529, 7.663 | no | 0.960 | 0.960 | 0.962 | 1.371 | descriptive (5 weeks) |
| 12 w | 0.855 | 9.893, 7.314, 8.774 | 9.670, 10.796, 9.924 | no | 0.855 | 0.983 | 0.983 | 0.818 | no evidence |
| 18 w | 0.787 | 6.681, 6.072, 5.312 | 6.585, 8.181, 8.197 | no | 0.802 | 0.983 | 0.998 | 0.764 | no evidence |
| 24 w | 0.864 | 8.043, 9.386, 9.989 | 9.832, 10.860, 11.042 | no | 0.861 | 0.963 | 0.962 | 1.207 | no evidence |
| 30 w | 0.831 | 6.395, 6.954, 6.428 | 7.406, 8.153, 8.244 | yes | 0.825 | 0.955 | 0.950 | 0.941 | candidate reduction of transcriptional output per Purkinje nucleus |

Caveats: nuclear UMI reflects nuclear transcript content and capture, not membrane area; the Purkinje share reflects dissociation as well as cell loss.

P7 rule: at 12 weeks or later, if KI/WT R_gran < 0.85 and all three KI mice are below the lowest WT mouse of that age, report “candidate reduction of transcriptional output per Purkinje nucleus”. Otherwise report “no evidence”. 5 weeks is descriptive. Caveats to state: nuclear UMI reflects nuclear transcript content and capture, not membrane area; the Purkinje share reflects dissociation as well as cell loss.
