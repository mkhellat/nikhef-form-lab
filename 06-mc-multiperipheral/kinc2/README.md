# kinc2 (vendored, third-party)

These files are copied verbatim from `kinc2.tar.gz`, the reference
implementation handed out with uam19 `part5.pdf` and described in

    J.A.M. Vermaseren, Nucl. Phys. B229 (1983) 347-371

They are **not our code** — nothing here was written for this repo. They
are vendored so that the `kinc2-driver/` programs are self-contained,
buildable checks with no dependency on paths outside the repo. Every
other file in `notes/` and `README.md` that cites `pickin.c`, `orient.c`,
`pi0.c`, `gamgam.c`, `peripp.c` etc. quotes line-numbered excerpts from
this exact copy.

## Included subset

| Files | Used by |
|---|---|
| `axolib.h`, `orient.c`, `pickin.c`, `mapt.c`, `mapw.c`, `mapla.c`, `ipow.c`, `iipow.c`, `ranf.c` | all drivers |
| `pi0.c` | π⁰ (2→3) drivers (`scan_pi0_scaling`, `check_gram_minors`, `mc_xsec_s_scan`) |
| `gamgam.c`, `peripp.c` | eemm (2→4) driver (`check_eemm_gamgam`) |

The full archive (including `eepi.c`, `eemm.c`, `eee.c`, `mgoto2.c`,
`vegas.c`, and the plotting driver code) is not reproduced here —
`eemm.c` itself is only a thin VEGAS callback around `gamgam`+`peripp`;
our driver calls those two directly.
