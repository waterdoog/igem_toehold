import os
import sys
import tempfile
import subprocess
import random
import re
import argparse
import shutil
from typing import Tuple, List

# =============================================================================
# Trigger RNA Design Tool for Toehold Switches
# =============================================================================
# This program designs trigger RNA sequences that activate toehold riboswitches.
# Toehold switches are RNA devices that control gene expression in response to
# specific trigger RNAs. The trigger binds to the switch, exposing the ribosome
# binding site (RBS) and start codon (AUG) to activate translation.
# =============================================================================

# ------------------ Design Constraints ------------------ #
# Biological constraints for functional RNA sequences
NUCLEOTIDES = ['A', 'C', 'G', 'U']  # RNA nucleotides (note: U not T)
GC_MIN = 0.40      # Minimum GC content (40%) for stable folding
GC_MAX = 0.55      # Maximum GC content (55%) to avoid over-stability
DG_MIN = -2.0      # Minimum free energy threshold (kcal/mol)
DG_MAX = 0.0       # Maximum free energy threshold (kcal/mol)
MIN_MFE_FREQ = 0.3 # Minimum frequency of minimum free energy structure
MIN_ENSEMBLE_DIVERSITY = 2.0  # Minimum structural diversity in ensemble
RBS_PATTERN = 'AGGAGG'  # Ribosome binding site sequence (Shine-Dalgarno)
AUG = 'AUG'        # Start codon for translation initiation

# ------------------ Utilities ------------------ #
def gc_content(seq):
    """Calculate GC content as fraction (0-1)"""
    gc = seq.count('G') + seq.count('C')
    return gc / len(seq)

def has_4plus_repeats(seq):
    """Check if sequence has 4 or more consecutive identical nucleotides"""
    return bool(re.search(r'(A{4,}|C{4,}|G{4,}|U{4,})', seq))

def generate_candidate_gc_controlled(n):
    """
    Generate random RNA sequence with controlled GC content and no long repeats

    Args:
        n: Desired sequence length

    Returns:
        Valid RNA sequence meeting constraints

    Raises:
        RuntimeError: If no valid sequence found after 1000 attempts
    """
    for _ in range(1000):
        seq = []
        gc_count = 0

        for i in range(n):
            possible = NUCLEOTIDES[:]  # Start with all nucleotides

            # Avoid 4+ consecutive repeats
            if len(seq) >= 3 and seq[-1] == seq[-2] == seq[-3]:
                possible.remove(seq[-1])

            # Control GC content dynamically during generation
            if i > 0:
                gc_ratio = gc_count / i
                if gc_ratio > GC_MAX:
                    # Too GC-rich: favor A/U
                    possible = [b for b in possible if b not in ['G', 'C']]
                elif gc_ratio < GC_MIN:
                    # Too GC-poor: favor G/C
                    possible = [b for b in possible if b in ['G', 'C']] or possible

            if not possible:
                break

            b = random.choice(possible)
            seq.append(b)
            if b in ['G', 'C']:
                gc_count += 1

        final_seq = ''.join(seq)
        final_gc = gc_count / n
        # Check all constraints: length, GC content, no long repeats
        if len(final_seq) == n and GC_MIN <= final_gc <= GC_MAX and not has_4plus_repeats(final_seq):
            return final_seq
    raise RuntimeError("Unable to generate valid sequence after 1000 attempts.")

# ------------------ RNAfold Integration ------------------ #
def run_rnafold_p(seq, tempdir):
    """
    Run RNAfold to predict RNA secondary structure and thermodynamic properties

    Args:
        seq: RNA sequence to fold
        tempdir: Temporary directory for RNAfold output files

    Returns:
        mfe: Minimum free energy (kcal/mol) - lower = more stable
        structure: Dot-bracket notation of predicted structure
        ensemble_div: Ensemble diversity (structural flexibility)
        mfe_freq: Frequency of MFE structure in thermodynamic ensemble

    Raises:
        RuntimeError: If RNAfold fails to run
        ValueError: If RNAfold output cannot be parsed
    """
    seq_file = os.path.join(tempdir, "seq.fa")
    with open(seq_file, 'w') as f:
        f.write(">seq\n" + seq + '\n')

    # RNAfold with partition function (-p) for ensemble properties
    cmd = ['RNAfold', '-p', '--noPS', '--noLP']
    proc = subprocess.run(cmd, input=seq, capture_output=True, text=True, cwd=tempdir)
    if proc.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {proc.stderr}")

    lines = proc.stdout.strip().split('\n')
    structure = None
    mfe = None

    # Parse MFE structure line: GAAAACCCUUU (((((...)))) -2.30
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
    # Parse ensemble properties from partition function output
    for line in lines:
        m_freq = re.search(r'frequency of mfe structure in ensemble\s*([\d\.]+)', line)
        ed = re.search(r'ensemble diversity\s*([\d\.]+)', line)
        if m_freq:
            mfe_freq = float(m_freq.group(1))
        if ed:
            ensemble_div = float(ed.group(1))

    return mfe, structure, ensemble_div, mfe_freq

def fold_trigger_plus_switch(trigger, switch, tempdir):
    """
    Fold combined trigger + switch sequence to check interaction

    Args:
        trigger: Trigger RNA sequence
        switch: Toehold switch sequence
        tempdir: Temporary directory for RNAfold

    Returns:
        combined_seq: Concatenated sequence
        structure: Predicted structure of trigger-switch complex
        mfe: Minimum free energy of the complex
    """
    combined_seq = trigger + switch
    mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(combined_seq, tempdir)
    return combined_seq, structure, mfe

def is_rbs_exposed(structure, trigger_len, switch_seq):
    """
    Check if RBS and AUG are exposed (single-stranded) in the complex

    In toehold switches, the trigger binds to expose the RBS and start codon
    for ribosome access. Single-stranded regions are marked with '.' in
    dot-bracket notation.

    Args:
        structure: Dot-bracket notation of trigger-switch complex
        trigger_len: Length of trigger sequence
        switch_seq: Switch sequence (must contain RBS and AUG)

    Returns:
        True if both RBS and AUG are in single-stranded regions
    """
    rbs_match = re.search(RBS_PATTERN, switch_seq)
    aug_match = re.search(AUG, switch_seq)
    if not rbs_match or not aug_match:
        return False
    # Calculate positions in the combined sequence
    rbs_start = trigger_len + rbs_match.start()
    rbs_end   = trigger_len + rbs_match.end()
    aug_start = trigger_len + aug_match.start()
    aug_end   = trigger_len + aug_match.end()
    # Check that all positions are single-stranded ('.')
    return all(structure[i] == '.' for i in range(rbs_start, rbs_end)) and \
           all(structure[i] == '.' for i in range(aug_start, aug_end))

def score_candidate(mfe, gc, mfe_freq, ensemble_div, no_4_repeats, rbs_exposed, n):
    """
    Score candidate trigger RNA based on multiple criteria

    A good trigger RNA should:
    - Have stable but not too stable folding (MFE -2 to 0 kcal/mol)
    - Have optimal GC content (~47.5%)
    - Have frequent MFE structure in ensemble (>30%)
    - Have structural diversity
    - Have no long repeats
    - Expose RBS and AUG when bound to switch

    Args:
        mfe: Minimum free energy (kcal/mol)
        gc: GC content percentage (0-100)
        mfe_freq: MFE structure frequency in ensemble (0-1)
        ensemble_div: Ensemble diversity measure
        no_4_repeats: Boolean, no 4+ consecutive nucleotides
        rbs_exposed: Boolean, RBS and AUG are accessible
        n: Sequence length

    Returns:
        Score (0-20+), higher is better. Returns -inf if invalid.
    """
    # Strict filtering: reject if any constraint violated
    if (mfe > 0 or mfe < -2 or gc < 40 or gc > 55 or
        not no_4_repeats or mfe_freq is None or mfe_freq < 0.3 or not rbs_exposed):
        return float('-inf')

    # Normalize ensemble diversity by sequence length
    norm_ensemble = ensemble_div / (n ** 1.5) if ensemble_div else 0

    # Score components (0-1 each, then weighted)
    mfe_score = 1 - abs(mfe / 2)        # Optimal MFE around -1 kcal/mol
    gc_score = 1 - abs(gc - 47.5) / 15  # Optimal GC around 47.5%

    # Weighted combination: MFE(40%), GC(20%), Ensemble(20%), MFE_freq(20%)
    return 10 * mfe_score + 5 * gc_score + 5 * norm_ensemble + 5 * mfe_freq

# ------------------ Main Program ------------------ #
def main():
    """
    Command-line interface for trigger RNA design

    Usage: python trigger_mrna_generator.py 36 --switch AGGAGGUAAAUG...

    This generates 10 high-scoring trigger RNA sequences that should activate
    the provided toehold switch by exposing the RBS and start codon.
    """
    parser = argparse.ArgumentParser(description="Trigger RNA design tool for toehold switches")
    parser.add_argument('n', type=int, help='Length of trigger RNA (e.g. 36)')
    parser.add_argument('--switch', type=str, required=True,
                       help='Toehold switch sequence (must contain RBS "AGGAGG" and AUG start codon)')
    args = parser.parse_args()

    n = args.n
    switch_seq = args.switch.upper()

    # Check RNAfold dependency
    if not shutil.which('RNAfold'):
        print("❌ RNAfold not found in PATH. Install ViennaRNA first.", file=sys.stderr)
        print("   On Ubuntu/Debian: sudo apt install vienna-rna", file=sys.stderr)
        print("   On macOS: brew install vienna-rna", file=sys.stderr)
        sys.exit(1)

    print(f"🔬 Designing trigger RNAs of length {n} for toehold switch...")
    print(f"   Target switch: {switch_seq[:20]}{'...' if len(switch_seq) > 20 else ''}")
    print(f"🔁 Searching until 10 valid candidates are found...")

    candidates = []
    attempts = 0
    max_attempts = 100000  # Prevent infinite loops

    while len(candidates) < 10 and attempts < max_attempts:
        attempts += 1
        try:
            # Generate random sequence with constraints
            seq = generate_candidate_gc_controlled(n)
        except Exception as e:
            print(f"⚠️ Sequence generation failed: {e}")
            continue

        # Evaluate sequence properties
        gc = 100 * gc_content(seq)  # Convert to percentage
        no4 = not has_4plus_repeats(seq)

        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                # Predict RNA folding
                mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(seq, tmpdir)
                # Check interaction with switch
                _, combined_structure, _ = fold_trigger_plus_switch(seq, switch_seq, tmpdir)
                rbs_exp = is_rbs_exposed(combined_structure, n, switch_seq)
            except Exception as e:
                print(f"⚠️ RNA folding analysis failed: {e}")
                continue

        # Score the candidate (strict criteria)
        score = score_candidate(mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp, n)
        if score >= 14:  # High-quality threshold
            candidates.append((score, seq, mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp))
            print(f"✅ Found candidate {len(candidates)}/10 (score: {score:.2f}) | Attempts: {attempts:,}")

    if len(candidates) == 0:
        print("❌ No valid trigger RNAs found after 100,000 attempts.")
        print("   Try adjusting the switch sequence or using a different length.")
        return

    # Rank by score and display results
    candidates.sort(key=lambda x: x[0], reverse=True)
    top10 = candidates[:10]

    print("\n🏆 Top 10 Trigger RNA Candidates:")
    print(f"{'Rank':<4} {'Sequence':<{n}} {'MFE':>6} {'GC%':>5} {'MFEfreq':>7} {'EnsDiv':>8} {'No4+':>7} {'RBSExp':>7} {'Score':>7}")
    for i, c in enumerate(top10, 1):
        score, seq, mfe, gc, mfe_freq, ens_div, no4, rbs = c
        print(f"{i:<4} {seq:<{n}} {mfe:6.2f} {gc:5.1f} {mfe_freq:7.2f} {ens_div:8.2f} {str(no4):7} {str(rbs):7} {score:7.2f}")

    print(f"\n📊 Summary: Found {len(candidates)} high-quality candidates in {attempts:,} attempts")

if __name__ == '__main__':
    main()


def trigger_candidate(n, switch_seq):
    """
    Find the single best trigger RNA candidate for a toehold switch

    This function performs the same search as main() but returns only the
    highest-scoring trigger RNA sequence and its properties.

    Args:
        n: Length of trigger RNA to design
        switch_seq: Toehold switch sequence containing RBS and AUG

    Returns:
        Tuple of (score, sequence, mfe, gc, mfe_freq, ensemble_div, no_repeats, rbs_exposed)
        Returns None if no valid candidates found
    """
    if not shutil.which('RNAfold'):
        print("❌ RNAfold not found in PATH. Install ViennaRNA first.", file=sys.stderr)
        return None

    print(f"🔬 Searching for optimal trigger RNA of length {n}...")

    candidates = []
    attempts = 0
    max_attempts = 100000

    while len(candidates) < 10 and attempts < max_attempts:
        attempts += 1
        try:
            seq = generate_candidate_gc_controlled(n)
        except Exception as e:
            print(f"⚠️ Sequence generation failed: {e}")
            continue

        gc = 100 * gc_content(seq)
        no4 = not has_4plus_repeats(seq)

        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(seq, tmpdir)
                _, combined_structure, _ = fold_trigger_plus_switch(seq, switch_seq, tmpdir)
                rbs_exp = is_rbs_exposed(combined_structure, n, switch_seq)
            except Exception as e:
                print(f"⚠️ RNA structure prediction failed: {e}")
                continue

        score = score_candidate(mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp, n)
        if score >= 14:  # High quality threshold
            candidates.append((score, seq, mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp))
            print(f"✅ Found candidate {len(candidates)}/10 | Score: {score:.2f} | Attempts: {attempts:,}")

    if len(candidates) == 0:
        print("❌ No valid trigger RNAs found.")
        return None

    # Return the highest-scoring candidate
    candidates.sort(key=lambda x: x[0], reverse=True)
    top = candidates[0]
    score, seq, mfe, gc, mfe_freq, ens_div, no4, rbs = top

    print("\n🏆 Best Trigger RNA Found:")
    print(f"   Sequence: {seq}")
    print(f"   Score: {score:.2f}")
    print(f"   Properties: MFE={mfe:.2f} kcal/mol, GC%={gc:.1f}, MFEfreq={mfe_freq:.2f}, EnsDiv={ens_div:.2f}")
    return top


# ------------------ Fast Trigger Generator (returns 1st valid) ------------------ #
def trigger_candidate_fast(n, switch_seq):
    """
    Quickly find the first valid trigger RNA candidate

    This function returns the first trigger RNA that meets the quality criteria,
    rather than searching for multiple candidates or the absolute best one.
    Useful when you just need any working trigger RNA sequence.

    Args:
        n: Length of trigger RNA to design
        switch_seq: Toehold switch sequence containing RBS and AUG

    Returns:
        First valid candidate tuple, or None if search fails
    """
    if not shutil.which('RNAfold'):
        print("❌ RNAfold not found in PATH. Install ViennaRNA first.", file=sys.stderr)
        return None

    attempts = 0
    max_attempts = 100000

    while attempts < max_attempts:
        attempts += 1
        try:
            # Generate random sequence with basic constraints
            seq = generate_candidate_gc_controlled(n)
        except Exception:
            continue

        # Quick evaluation of basic properties
        gc = 100 * gc_content(seq)
        no4 = not has_4plus_repeats(seq)

        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                # Predict RNA structure and check switch interaction
                mfe, structure, ensemble_div, mfe_freq = run_rnafold_p(seq, tmpdir)
                _, combined_structure, _ = fold_trigger_plus_switch(seq, switch_seq, tmpdir)
                rbs_exp = is_rbs_exposed(combined_structure, n, switch_seq)
            except Exception:
                continue

        # Score candidate (returns -inf if invalid)
        score = score_candidate(mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp, n)
        if score >= 14:  # Quality threshold
            print(f"✅ Found valid trigger RNA | Sequence: {seq} | Score: {score:.2f} | Attempts: {attempts:,}")
            return (score, seq, mfe, gc, mfe_freq, ensemble_div, no4, rbs_exp)

    print("❌ No valid trigger RNA found after 100,000 attempts.")
    return None
