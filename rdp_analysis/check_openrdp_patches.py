#!/usr/bin/env python3
"""
Refuse to run OpenRDP unless the local bug patches are applied.

The four fixes in patches/openrdp-local-fixes.patch live in the *installed*
openrdp package, not in this repo, so a fresh clone or a
`pip install --upgrade openrdp` silently reverts them. Running unpatched is not
a clean failure:

  - MaxChi crashes partway through every threshold (so you get no output at all)
  - SiScan loops forever, leaking ~5 GB/h until the kernel kills it, which took
    23 h and 40 h to happen the first time
  - MaxChi and Chimaera silently discard perfectly-associated windows, so any
    numbers you do get disagree with the committed results for no visible reason

That last one is the reason this check exists: it is the failure that looks like
success. Every runner calls this first.

Checks marker counts per file, not just presence, so a half-applied patch (e.g.
`patch` rejected one hunk) is caught too.

Usage:
    python3 check_openrdp_patches.py          # exit 0 if patched, 1 if not
    from check_openrdp_patches import require_patches; require_patches()
"""
import os
import sys

MARKER = "LOCAL PATCH"

# file -> (expected marker count, what those patches fix)
EXPECTED = {
    "common.py": (1, "calculate_chi2 validity guard (row/column sums, not the "
                     "original wrong condition)"),
    "maxchi.py": (1, "clamp left/right before indexing chi2_values"),
    "chimaera.py": (1, "clamp left/right before indexing chi2_values"),
    "siscan.py": (2, "find_signal index must not move backwards -- the infinite "
                     "loop that caused the OOM kills (one fix per maj1/maj2 loop)"),
}

REMEDY = """
To fix, apply the patch to the installed package (see rdp_analysis/README.md):

    cd "$(python3 -c 'import openrdp, os; print(os.path.dirname(os.path.dirname(openrdp.__file__)))')"
    patch -p1 < {patch}

Then re-run this check. Do not work around it by editing this script.
"""


def check():
    """Return (ok, list of problem strings)."""
    try:
        import openrdp
    except ImportError:
        return False, ["openrdp is not installed or not importable"]

    pkg = os.path.dirname(openrdp.__file__)
    problems = []
    for fname, (want, what) in sorted(EXPECTED.items()):
        path = os.path.join(pkg, fname)
        if not os.path.exists(path):
            problems.append(f"{fname}: missing from {pkg}")
            continue
        with open(path, encoding="utf-8", errors="replace") as handle:
            got = handle.read().count(MARKER)
        if got != want:
            problems.append(
                f"{fname}: found {got} '{MARKER}' marker(s), expected {want} -- {what}")
    return not problems, problems


def require_patches():
    """Exit non-zero with an explanation if the patches are not in place."""
    ok, problems = check()
    if ok:
        return
    patch = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "patches", "openrdp-local-fixes.patch")
    sys.stderr.write("ERROR: the local OpenRDP patches are not applied.\n\n")
    for p in problems:
        sys.stderr.write(f"  - {p}\n")
    sys.stderr.write(REMEDY.format(patch=patch))
    sys.exit(1)


if __name__ == "__main__":
    require_patches()
    import openrdp
    print(f"OpenRDP patches OK ({os.path.dirname(openrdp.__file__)})")
