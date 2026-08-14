# kinc2 (vendored, third-party)

These files are copied verbatim from `kinc2.tar.gz`, the reference
implementation handed out with uam19 `part5.pdf` and described in

    J.A.M. Vermaseren, Nucl. Phys. B229 (1983) 347-371

They are **not our code** — nothing here was written for this repo. They
are vendored so that `scripts/scan_pi0_scaling.c` is a self-contained,
buildable program with no dependency on paths outside the repo. Every
other file in `notes/` and `README.md` that cites `pickin.c`, `orient.c`,
`pi0.c` etc. quotes line-numbered excerpts from this exact copy.

Only the subset needed to build `scan_pi0_scaling.c` is included:
`axolib.h`, `orient.c`, `pickin.c`, `mapt.c`, `mapw.c`, `mapla.c`,
`ipow.c`, `iipow.c`, `ranf.c`, `pi0.c`. The full archive (including
`eepi.c`, `mgoto2.c`, `gamgam.c`, `vegas.c`, and the plotting/VEGAS
driver code) is not reproduced here.
