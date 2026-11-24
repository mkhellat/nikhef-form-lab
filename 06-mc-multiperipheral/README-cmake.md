# CMake build (Topic 6 C gates)

Real multi-platform CMake for the numeric evaluators and kinc2-driver
gates. This is what `.github/workflows/cmake-multiplatform.yml` builds
on Ubuntu, macOS, and Windows — not an empty CI stub.

## Configure / build / test

```bash
cd 06-mc-multiperipheral
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

On multi-config generators (Visual Studio):

```bash
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Targets

| Target | Role | CTest? |
|---|---|---|
| `check_stable_c` | (3.5) + ΣMij vs fixed refs | yes |
| `digit_loss_mij` | ΣMij double vs LD vs cosθ | yes |
| `check_eemm_gamgam` | 2→4 gamgam/peripp vs ΣMij | yes |
| `digit_loss_eemm_gamgam` | ΣMij vs L⊗L(naive) on gamgam | yes |
| `mc_digit_loss_eemm_gamgam` | gamgam frac_bad rate, soft/mid/hard |t| | yes |
| `digit_loss_eemm_cdc` | CDC-60 naive vs IEEE; ΣMij stays healthy | yes |
| `digit_loss_naive` / `mc_digit_loss` | Sec 6–7 scans | build only |

`digit_loss_eemm_gamgam` and `mc_digit_loss_eemm_gamgam` self-skip the
LD compare (print PASS skip) when `long double` is not wider than
`double` (typical MSVC), so the Windows job stays green without
pretending extended precision exists. `digit_loss_eemm_cdc` does not
skip: it compares the 48-bit CDC word to IEEE (long double when that
is wider, otherwise double) and still requires the soft ΣMᵢⱼ
`frac_bad` to stay low.

Makefiles under `numeric/` and `kinc2-driver/` remain the local
day-to-day path; CMake is the portable / CI path.
