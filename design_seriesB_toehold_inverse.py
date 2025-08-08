import subprocess
import random
from Bio.Seq import Seq

# -------------------- ViennaRNA Calls -------------------- #

def rnafold(seq):
    process = subprocess.Popen(['RNAfold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(seq)
    if process.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {stderr}")
    lines = stdout.strip().split('\n')
    struct = lines[1].split()[0]
    energy = lines[1].split('(')[-1].strip(')')
    return struct, float(energy), struct

def rnainverse(target_structure, start_seq=None):
    input_str = target_structure if not start_seq else f"{start_seq}\n{target_structure}"
    process = subprocess.Popen(['RNAinverse'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input_str)
    if process.returncode != 0:
        raise RuntimeError(f"RNAinverse failed: {stderr}")
    lines = stdout.strip().split('\n')
    for line in reversed(lines):
        parts = line.strip().split()
        if len(parts) >= 2:
            seq = parts[0]
            struct = parts[1]
            energy = None
            for p in parts:
                if p.startswith('(') and p.endswith(')'):
                    try:
                        energy = float(p.strip('()'))
                        break
                    except ValueError:
                        continue
            return seq, struct, energy if energy is not None else 0.0
    raise RuntimeError("Could not parse RNAinverse output.")

# -------------------- Sequence Utilities -------------------- #

def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

def random_rna_seq(n):
    return ''.join(random.choices('ACGU', k=n))

# -------------------- Toehold Generator -------------------- #

def generate_toehold_via_inverse(trigger_mrna, toxin_seq, max_attempts=100):
    trigger_len = len(trigger_mrna)
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    linker_seq = "AACCUGGCGGCAGCGCAAAAG"

    def build_target_structure_exact(trigger_len):
        return (
            '.' * trigger_len +
            '(((((((((((.........................))))))))))).' +
            '.' * len(linker_seq)
        )

    target_structure = build_target_structure_exact(trigger_len)
    rc_trigger = reverse_complement_rna(trigger_mrna)

    for attempt in range(max_attempts):
        try:
            seq, inverse_struct, _ = rnainverse(target_structure)
        except Exception as e:
            print(f"⚠️ RNAinverse failed at attempt {attempt+1}: {e}")
            continue

        if seq[:trigger_len] != rc_trigger:
            seq = rc_trigger + seq[trigger_len:]

        struct, energy, folded_struct = rnafold(seq)
        if folded_struct != target_structure:
            print(f"❌ Attempt {attempt+1}: RNAfold structure mismatch.")
            continue

        conserved_start = trigger_len + 11
        linker_start = trigger_len + 25 + 11 + 1

        actual_conserved = seq[conserved_start:conserved_start + len(conserved_seq)]
        actual_linker = seq[linker_start:linker_start + len(linker_seq)]

        if actual_conserved != conserved_seq:
            print(f"❌ Attempt {attempt+1}: Conserved seq mismatch: {actual_conserved}")
            continue

        if actual_linker != linker_seq:
            print(f"❌ Attempt {attempt+1}: Linker seq mismatch: {actual_linker}")
            continue

        print(f"✅ Found valid toehold after {attempt+1} attempts.")
        return seq + toxin_seq

    print("❌ Failed to find valid toehold matching all constraints after max attempts.")
    return None

# -------------------- Design Driver -------------------- #

def design_seriesB_toehold(toxin_sequence, num_trials=10, trigger_length=25):
    results = []
    triggers = generate_trigger_candidates(toxin_sequence, num_trials=num_trials, trigger_len=trigger_length)

    for idx, trig in enumerate(triggers):
        print(f"\n🔍 Trying trigger #{idx + 1}")
        print(f"Trigger mRNA (5'->3'): {trig}")
        rc_trigger = reverse_complement_rna(trig)
        print(f"Reverse complement for toehold 5' region: {rc_trigger}")

        toehold_seq = generate_toehold_via_inverse(trig, toxin_sequence)
        if toehold_seq is None:
            print("❌ RNAinverse+RNAfold constraints failed. Skipping this trigger.\n")
            continue

        toehold_only = toehold_seq[:-len(toxin_sequence)]
        struct, energy, dotbracket = rnafold(toehold_seq)

        print(f"🧬 Toehold (no toxin): {toehold_only}")
        print(f"🧬 Full sequence: {toehold_seq}")
        print(f"📐 Structure:\n{dotbracket}")
        print(f"ΔG = {energy:.2f} kcal/mol")

        results.append({
            'trigger_mRNA': trig,
            'toehold_sequence': toehold_seq,
            'structure': dotbracket,
            'energy': energy
        })

    return results

# -------------------- Helpers -------------------- #

def generate_trigger_candidates(toxin_seq, num_trials=100, trigger_len=25):
    valid = []
    attempts = 0
    while len(valid) < num_trials and attempts < num_trials * 10:
        candidate = random_rna_seq(trigger_len)
        if validate_trigger_mrna(candidate):
            valid.append(candidate)
        attempts += 1
    return valid

def validate_trigger_mrna(seq, min_gc=0.4, max_hairpin_energy=-5.0):
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    if gc_content < min_gc:
        return False
    struct, energy, _ = rnafold(seq)
    return energy > max_hairpin_energy

# -------------------- CLI -------------------- #

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Design Series B toehold switches using RNAinverse.")
    parser.add_argument("--toxin", type=str, required=True, help="Toxin gene RNA sequence (5' to 3').")
    parser.add_argument("--trials", type=int, default=10, help="Number of trigger mRNA candidates to generate.")
    parser.add_argument("--trigger_length", type=int, default=25, help="Length of trigger mRNA.")
    args = parser.parse_args()

    toxin_sequence = args.toxin.upper().replace('T', 'U').replace(' ', '').replace('\n','')
    designs = design_seriesB_toehold(toxin_sequence, num_trials=args.trials, trigger_length=args.trigger_length)

    if not designs:
        print("❌ No valid toehold switches found.")
    else:
        for i, res in enumerate(designs):
            print(f"\n=== ✅ Final Toehold design {i+1} ===")
            print(f"Trigger mRNA (5'->3'): {res['trigger_mRNA']}")
            print(f"Complementary Toehold Trigger: {reverse_complement_rna(res['trigger_mRNA'])}")
            print(f"Toehold Switch + Toxin:\n{res['toehold_sequence']}")
            print(f"Structure:\n{res['structure']}")
            print(f"ΔG = {res['energy']:.2f} kcal/mol")
