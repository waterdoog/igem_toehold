import os
import sys
import tempfile
import subprocess
import random
import re
import argparse
import shutil

# ------------------ Design Constraints ------------------ #
NUCLEOTIDES = ['A', 'C', 'G', 'U']
GC_MIN = 0.40
GC_MAX = 0.55
DG_MIN = -2.0
DG_MAX = 0.0
MIN_MFE_FREQ = 0.3
MIN_ENSEMBLE_DIVERSITY = 2.0
RBS_PATTERN = 'AGGAGG'
AUG = 'AUG'

# ------------------ Utilities ------------------ #
def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)

def has_4plus_repeats(seq):
    return bool(re.search(r'(A{4,}|C{4,}|G{4,}|U{4,})', seq))

def generate_candidate_gc_controlled(n):
    for _ in range(1000):
        seq = []
        gc_count = 0

        for i in range(n):
            possible = NUCLEOTIDES[:]

            if len(seq) >= 3 and seq[-1] == seq[-2] == seq[-3]:
                possible.remove(seq[-1])

            if i > 0:
                gc_ratio = gc_count / i
                if gc_ratio > GC_MAX:
                    possible = [b for b in possible if b not in ['G', 'C']]
                elif gc_ratio < GC_MIN:
                    possible = [b for b in possible if b in ['G', 'C']] or possible

            if not possible:
                break

            b = random.choice(possible)
            seq.append(b)
            if b in ['G', 'C']:
                gc_count += 1

        final_seq = ''.join(seq)
        final_gc = gc_count / n
        if len(final_seq) == n and GC_MIN <= final_gc <= GC_MAX and not has_4plus_repeats(final_seq):
            return final_seq
    raise RuntimeError("Unable to generate valid sequence after 1000 attempts.")

# ------------------ RNAfold Integration ------------------ #
def run_rnafold_p(seq, tempdir):
    seq_file = os.path.join(tempdir, "seq.fa")
    with open(seq_file, 'w') as f:
        f.write(">seq\n" + seq + '\n')

    cmd = ['RNAfold', '-p', '--noPS', '--noLP']
    proc = subprocess.run(cmd, input=seq, capture_output=True, text=True, cwd=tempdir)
    if proc.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {proc.stderr}")

    lines = proc.stdout.strip().split('\n')
    structure = None
    mfe = None

    for line in lines:
        m = re.search(r'^([\.()\[\]{}]+)\s+\(\s*(-?[\d\.]+)\s*\)', line)
        if m:
            structure = m.group(1)
            mfe = float(m.group(2))
            break
    if structure is None or mfe is None:
        raise ValueError("Failed to parse structure and MFE")

    mfe_freq = None
    ensemble_div = None
    for line in lines:
        m_freq = re.search(r'frequency of mfe structure in ensemble\s*([\d\.]+)', line)
        ed = re.search(r'ensemble diversity\s*([\d\.]+)', line)
        if m_freq:
            mfe_freq = float(m_freq.group(1))
        if ed:
            ensemble_div = float(ed.group(1))

    return mfe, structure, ensemble_div, mfe_freq

def fold_trigger_plus_switch(trigger, switch, tempdir):
    combined_seq = trigger + switch
    mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(combined_seq, tempdir)
    return combined_seq, structure, mfe

def is_rbs_exposed(structure, trigger_len, switch_seq):
    rbs_match = re.search(RBS_PATTERN, switch_seq)
    aug_match = re.search(AUG, switch_seq)
    if not rbs_match or not aug_match:
        return False
    rbs_start = trigger_len + rbs_match.start()
    rbs_end   = trigger_len + rbs_match.end()
    aug_start = trigger_len + aug_match.start()
    aug_end   = trigger_len + aug_match.end()
    return all(structure[i] == '.' for i in range(rbs_start, rbs_end)) and \
           all(structure[i] == '.' for i in range(aug_start, aug_end))

def score_candidate(mfe, gc, mfe_freq, ensemble_div, no_4_repeats, rbs_exposed, n):
    if mfe > 0 or mfe < -2 or gc < 40 or gc > 55 or not no_4_repeats or mfe_freq is None or mfe_freq < 0.3 or not rbs_exposed:
        return float('-inf')
    norm_ensemble = ensemble_div / (n ** 1.5) if ensemble_div else 0
    mfe_score = 1 - abs(mfe / 2)
    gc_score = 1 - abs(gc - 47.5) / 15
    return 10 * mfe_score + 5 * gc_score + 5 * norm_ensemble + 5 * mfe_freq

# ------------------ Main Program ------------------ #
def main():
    parser = argparse.ArgumentParser(description="Trigger RNA design tool")
    parser.add_argument('n', type=int, help='Length of trigger RNA (e.g. 36)')
    parser.add_argument('--switch', type=str, required=True, help='Switch sequence (must contain RBS and AUG)')
    args = parser.parse_args()

    n = args.n
    switch_seq = args.switch.upper()

    if not shutil.which('RNAfold'):
        print("❌ RNAfold not found in PATH. Install ViennaRNA first.", file=sys.stderr)
        sys.exit(1)

    print(f"🔁 Running until 10 valid trigger RNA candidates of length {n} are found...")

    candidates = []
    attempts = 0
    max_attempts = 100000

    while len(candidates) < 10 and attempts < max_attempts:
        attempts += 1
        try:
            seq = generate_candidate_gc_controlled(n)
        except Exception as e:
            print(f"⚠️ Generation failed: {e}")
            continue

        gc = 100 * gc_content(seq)
        no4 = not has_4plus_repeats(seq)

        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(seq, tmpdir)
                _, combined_structure, _ = fold_trigger_plus_switch(seq, switch_seq, tmpdir)
                rbs_exp = is_rbs_exposed(combined_structure, n, switch_seq)
            except Exception as e:
                print(f"⚠️ RNAfold error: {e}")
                continue

        score = score_candidate(mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp, n)
        if score >= 14:
            candidates.append((score, seq, mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp))
            print(f"✅ {len(candidates)}/10 candidates found (score = {score:.2f})")

    if len(candidates) == 0:
        print("❌ No valid candidates found.")
        return

    candidates.sort(key=lambda x: x[0], reverse=True)
    top10 = candidates[:10]

    print("\n🏆 Top 10 Trigger RNA Candidates:")
    print(f"{'Rank':<4} {'Sequence':<{n}} {'MFE':>6} {'GC%':>5} {'MFEfreq':>7} {'EnsDiv':>8} {'No4+':>7} {'RBSExp':>7} {'Score':>7}")
    for i, c in enumerate(top10, 1):
        score, seq, mfe, gc, mfe_freq, ens_div, no4, rbs = c
        print(f"{i:<4} {seq:<{n}} {mfe:6.2f} {gc:5.1f} {mfe_freq:7.2f} {ens_div:8.2f} {str(no4):7} {str(rbs):7} {score:7.2f}")

if __name__ == '__main__':
    main()
