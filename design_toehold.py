import subprocess
import random
from Bio.Seq import Seq
from trigger_mrna_generator import trigger_candidate_fast
import RNA  

def rnafold(seq):
    """
    Fold RNA sequence using ViennaRNA Python API.
    Returns: (structure, mfe_energy, structure)
    """
    structure, mfe = RNA.fold(seq)
    return structure, mfe, structure

def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

def is_region_unpaired(dot_bracket, start, length):
    region = dot_bracket[start:start+length]
    return all(c == '.' for c in region)

def no_interaction_rnacofold(seq1, seq2):
    """
    Returns True if no base-pair interaction between seq1 and seq2 is predicted.
    Uses ViennaRNA Python bindings (no subprocess).
    """
    input_str = f"{seq1}&{seq2}"
    
    # Cofolding using RNA.cofold()
    fc = RNA.cofold(input_str)
    struct, mfe = fc

    # Find index of '&' separator
    amp_index = len(seq1)

    # Map base pairs using stack
    pair_stack = []
    pairs = {}

    for idx, c in enumerate(struct):
        if c == '(':
            pair_stack.append(idx)
        elif c == ')':
            if not pair_stack:
                raise ValueError("Unbalanced dot-bracket structure")
            start_idx = pair_stack.pop()
            pairs[start_idx] = idx
            pairs[idx] = start_idx

    # Check for cross-strand base pairing (across `&`)
    for i in range(amp_index):
        paired_idx = pairs.get(i)
        if paired_idx is not None and paired_idx > amp_index:
            return False  # cross-strand interaction detected

    return True  # no interaction

def validate_toxin(toxin_seq, toehold_seq):
    struct, energy, _ = rnafold(toxin_seq)
    paired_nt = sum(1 for c in struct if c != '.')
    max_paired = int(0.3 * len(toxin_seq))
    if paired_nt > max_paired:
        return False
    if not no_interaction_rnacofold(toehold_seq, toxin_seq):
        return False
    return True

def generate_toehold(trigger_mrna, toxin_seq):
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    toehold_trigger = reverse_complement_rna(trigger_mrna)
    full_seq = toehold_trigger + conserved_seq + universal_linker + toxin_seq
    return full_seq

def validate_toehold_structure(full_seq, toxin_len):
    struct, energy, dotbracket = rnafold(full_seq)
    toehold_len = len(full_seq) - len("GGACUUUAGAACAGAGGAGAUAAAGAUG") - len("AACCUGGCGGCAGCGCAAAAG") - toxin_len
    conserved_start = toehold_len
    conserved_end = conserved_start + 25
    universal_linker_start = conserved_end + 11
    universal_linker_end = universal_linker_start + 22
    if not is_region_unpaired(dotbracket, conserved_start, 25):
        return False
    if not is_region_unpaired(dotbracket, universal_linker_start, 22):
        return False
    return True

def design_seriesB_toehold(toxin_seq, trigger_length=25):
    attempts = 0
    output_path = "/mnt/c/Users/tsdn/Desktop/Igem/igem_toehold/toehold_results.txt"
    
    while True:
        if attempts % 1000 == 0 and attempts > 0:
            print(f"📣 Attempted {attempts} trigger+switch designs...")

        best = trigger_candidate_fast(trigger_length, "GGAAGGAGGUAACAAUG")
        attempts += 1
        if best is None:
            continue

        best_seq, best_score = best[1], best[0]
        print(f"[Attempt {attempts}] Best sequence: {best_seq}, Score: {best_score:.2f}")

        rev_comp = reverse_complement_rna(best_seq)
        struct, energy, _ = rnafold(rev_comp)
        if sum(1 for c in struct if c != '.') > int(len(rev_comp)):
            print("❌ Reverse complement is too structured (not linear)")
            continue

        toehold_seq = generate_toehold(best_seq, toxin_seq)
        if not validate_toehold_structure(toehold_seq, len(toxin_seq)):
            print("❌ Toehold structure invalid")
            continue
        if not validate_toxin(toxin_seq, toehold_seq):
            print("❌ Toxin structure or interaction invalid")
            continue

        print("✅ Valid toehold switch found!")
        with open(output_path, "w") as f:
            f.write(f"Trigger mRNA (5'->3'): {best_seq}\n")
            f.write(f"Complementary Toehold Trigger region (5'->3'): {rev_comp}\n")
            f.write(f"Score: {best_score:.2f}\n")
            f.write(f"Full Toehold Switch + Toxin sequence:\n{toehold_seq}\n")

        return [{
            'trigger_mRNA': best_seq,
            'toehold_sequence': toehold_seq,
            'score': best_score
        }]

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Design Series B toehold switches with toxin gene input.")
    parser.add_argument("--toxin", type=str, required=True, help="Toxin gene RNA sequence (5' to 3').")
    parser.add_argument("--trigger_length", type=int, default=25, help="Length of trigger mRNA sequence.")
    args = parser.parse_args()

    toxin_sequence = args.toxin.upper().replace('T', 'U').replace(' ', '').replace('\n','')
    trigger_len = args.trigger_length

    print(f"Starting design with toxin gene length {len(toxin_sequence)}, trigger length {trigger_len}")

    designs = design_seriesB_toehold(toxin_sequence, trigger_length=trigger_len)

    if not designs:
        print("No valid toehold switches found with given constraints.")
    else:
        for i, res in enumerate(designs):
            print(f"\n--- Toehold design {i+1} ---")
            print(f"Trigger mRNA (5'->3'): {res['trigger_mRNA']}")
            comp_trigger = reverse_complement_rna(res['trigger_mRNA'])
            print(f"Complementary Toehold Trigger region (5'->3'): {comp_trigger}")
            print(f"Full Toehold Switch + Toxin seq:\n{res['toehold_sequence']}")
