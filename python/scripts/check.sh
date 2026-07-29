#!/usr/bin/env bash
# check.sh — end-to-end sanity check for polybuilder.
#
# Runs:
#   1. pytest (unit tests, with coverage summary)
#   2. ruff   (lint)
#   3. mypy   (type check)
#   4. CLI smoke tests: three real builds using the CLI, verifying that
#      the expected residue blocks appear in the generated .rtp file.
#
# Prints a green PASS / red FAIL summary at the end and exits non-zero
# if anything failed.

set -u   # unset variables are errors; we DO NOT set -e because we want to
         # keep running after a failing step so the summary is complete.

# -------- resolve script dir + project root --------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT"

# Make src/ importable whether polybuilder is installed or not.
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

# -------- pretty printing --------
if [[ -t 1 ]]; then
    GREEN='\033[0;32m'; RED='\033[0;31m'; YELLOW='\033[0;33m'; BOLD='\033[1m'; RESET='\033[0m'
else
    GREEN=''; RED=''; YELLOW=''; BOLD=''; RESET=''
fi

declare -a STEPS=()
declare -a RESULTS=()
FAILURES=0

run_step() {
    local label="$1"; shift
    echo
    echo -e "${BOLD}${YELLOW}==>${RESET} ${BOLD}${label}${RESET}"
    echo -e "    ${YELLOW}$ $*${RESET}"
    if "$@"; then
        echo -e "    ${GREEN}\xE2\x9C\x94 ${label} passed${RESET}"
        STEPS+=("$label")
        RESULTS+=("PASS")
    else
        echo -e "    ${RED}\xE2\x9C\x98 ${label} FAILED${RESET}"
        STEPS+=("$label")
        RESULTS+=("FAIL")
        FAILURES=$((FAILURES + 1))
    fi
}

# -------- CLI smoke test helper --------
SMOKE_DIR="$(mktemp -d -t polybuilder-check-XXXXXX)"
trap 'rm -rf "$SMOKE_DIR"' EXIT

# Runs the CLI and greps the RTP for expected residue names.
# Usage: cli_smoke <label> <subdir> <expected_pattern> -- <cli args...>
cli_smoke() {
    local label="$1"; shift
    local subdir="$1"; shift
    local expected="$1"; shift
    [[ "$1" == "--" ]] || { echo "cli_smoke bad usage"; return 2; }
    shift

    local out="$SMOKE_DIR/$subdir"
    mkdir -p "$out"

    echo
    echo -e "${BOLD}${YELLOW}==>${RESET} ${BOLD}CLI smoke: ${label}${RESET}"
    echo -e "    ${YELLOW}$ polybuilder $*${RESET}"

    if ! python -m polybuilder.cli "$@" --no-editconf --output-dir "$out" > "$out/stdout.log" 2> "$out/stderr.log"; then
        echo -e "    ${RED}\xE2\x9C\x98 CLI exited non-zero${RESET}"
        echo "    --- stderr ---"; sed 's/^/    /' "$out/stderr.log"
        STEPS+=("CLI smoke: $label"); RESULTS+=("FAIL"); FAILURES=$((FAILURES + 1))
        return
    fi

    local rtp="$out/rearranged_polymer.rtp"
    if [[ ! -s "$rtp" ]]; then
        echo -e "    ${RED}\xE2\x9C\x98 expected RTP not produced${RESET}"
        STEPS+=("CLI smoke: $label"); RESULTS+=("FAIL"); FAILURES=$((FAILURES + 1))
        return
    fi

    local missing=""
    for pat in $expected; do
        if ! grep -q "^\[ $pat \]" "$rtp"; then
            missing="$missing $pat"
        fi
    done
    if [[ -n "$missing" ]]; then
        echo -e "    ${RED}\xE2\x9C\x98 missing residue block(s):${missing}${RESET}"
        echo "    RTP residue lines:"; grep -E '^\[' "$rtp" | sed 's/^/      /'
        STEPS+=("CLI smoke: $label"); RESULTS+=("FAIL"); FAILURES=$((FAILURES + 1))
        return
    fi

    echo -e "    ${GREEN}\xE2\x9C\x94 CLI smoke: ${label} passed${RESET} (RTP: $rtp)"
    STEPS+=("CLI smoke: $label"); RESULTS+=("PASS")
}

# -------- environment sanity --------
echo -e "${BOLD}polybuilder check — $(date '+%Y-%m-%d %H:%M:%S')${RESET}"
echo -e "Working dir : $ROOT"
echo -e "Python      : $(python --version 2>&1)"
python -c "import polybuilder; print('polybuilder :', polybuilder.__version__)" 2>/dev/null \
    || echo -e "${YELLOW}polybuilder is not importable — did you run 'pip install -e .'?${RESET}"
python -c "import rdkit; print('rdkit       :', rdkit.__version__)" 2>/dev/null \
    || { echo -e "${RED}rdkit missing — install first: pip install rdkit${RESET}"; exit 2; }

# -------- 1. pytest --------
run_step "unit tests"         python -m pytest tests -q --tb=short

# -------- 2. ruff --------
if command -v ruff >/dev/null 2>&1; then
    run_step "ruff (lint)"    ruff check src tests
else
    echo -e "${YELLOW}ruff not installed — skipping lint. (pip install ruff)${RESET}"
fi

# -------- 3. mypy --------
if command -v mypy >/dev/null 2>&1; then
    run_step "mypy (types)"   mypy src
else
    echo -e "${YELLOW}mypy not installed — skipping type check. (pip install mypy)${RESET}"
fi

# -------- 4. CLI smoke tests --------
LIB="$ROOT/examples/cpp_monomers.py"

# 4a. MMA homopolymer with KPS on both ends.
cli_smoke "MMA + KPS(both)" "mma_kps" \
    "KPH MCF MCL KPT" -- \
    build --comonomer MMA --n-first 5 --n-comonomer 0 --repeats 1 --cap \
          --initiator KPS --ends both

# 4b. MMA/DMAPS sulfobetaine copolymer (uses --bridge on a sulfobetaine).
cli_smoke "MMA + DMAPS (bridge=3)" "mma_dmaps" \
    "MCF DMR MCR" -- \
    build --comonomer DMAPS --n-first 3 --n-comonomer 1 --repeats 2 --bridge 3

# 4c. NIPAM + DABCO copolymer with KPS ends (the two-arbitrary-residues case).
cli_smoke "NIPAM + DABCO + KPS(both)" "nipam_dabco" \
    "KPH NIF DAR KPT" -- \
    build --library "$LIB" \
          --first-residue NIPAM --n-first 3 \
          --comonomer     DABCO --n-comonomer 1 \
          --repeats 3 --cap \
          --initiator KPS --ends both

# 4d. list-monomers / list-initiators exercise the plugin loader.
run_step "list-monomers with plugin" \
    python -m polybuilder.cli list-monomers --library "$LIB"
run_step "list-initiators with plugin" \
    python -m polybuilder.cli list-initiators --library "$LIB"

# -------- summary --------
echo
echo -e "${BOLD}================ Summary ================${RESET}"
for i in "${!STEPS[@]}"; do
    step="${STEPS[$i]}"; result="${RESULTS[$i]}"
    if [[ "$result" == "PASS" ]]; then
        echo -e "  ${GREEN}\xE2\x9C\x94 PASS${RESET}  $step"
    else
        echo -e "  ${RED}\xE2\x9C\x98 FAIL${RESET}  $step"
    fi
done
echo

if [[ "$FAILURES" -eq 0 ]]; then
    echo -e "${BOLD}${GREEN}All ${#STEPS[@]} checks passed.${RESET}"
    echo "Sample outputs kept in: $SMOKE_DIR"
    trap - EXIT   # keep outputs around for inspection
    exit 0
else
    echo -e "${BOLD}${RED}${FAILURES} of ${#STEPS[@]} checks FAILED.${RESET}"
    echo "Sample outputs kept in: $SMOKE_DIR"
    trap - EXIT
    exit 1
fi
