import subprocess
import random
from Bio.Seq import Seq
from trigger_mrna_generator import trigger_candidate_fast
import RNA  

def rnafold(seq):
    structure, mfe = RNA.fold(seq)
    return structure, mfe, structure

def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

def is_region_unpaired(dot_bracket, start, length):
    region = dot_bracket[start:start+length]
    return all(c == '.' for c in region)

def no_interaction_rnacofold(seq1, seq2):
    input_str = f"{seq1}&{seq2}"
    fc = RNA.cofold(input_str)
    struct, mfe = fc
    amp_index = len(seq1)
    pair_stack = []
    pairs = {}
    for idx, c in enumerate(struct):
        if c == '(': pair_stack.append(idx)
        elif c == ')':
            if not pair_stack:
                raise ValueError("Unbalanced dot-bracket structure")
            start_idx = pair_stack.pop()
            pairs[start_idx] = idx
            pairs[idx] = start_idx
    for i in range(amp_index):
        paired_idx = pairs.get(i)
        if paired_idx is not None and paired_idx > amp_index:
            return False
    return True

def validate_toxin(toxin_seq, toehold_seq):
    struct, energy, _ = rnafold(toxin_seq)
    paired_nt = sum(1 for c in struct if c != '.')
    max_paired = int(0.3 * len(toxin_seq))
    if paired_nt > max_paired:
        print("📛 Failure reason: Toxin structure too folded (>30% paired).\n")
        return False
    if not no_interaction_rnacofold(toehold_seq, toxin_seq):
        print("📛 Failure reason: Toehold and toxin sequence interact.\n")
        return False
    return True

def generate_toehold(trigger_mrna, toxin_seq):
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    toehold_trigger = reverse_complement_rna(trigger_mrna)
    print (f"🧪 toehold_trigger: {toehold_trigger}", flush=True)
    full_seq = toehold_trigger + conserved_seq + universal_linker + toxin_seq
    return full_seq

def validate_toehold_structure(full_seq, toxin_len):
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"


    design_region = full_seq[:len(full_seq) - toxin_len]

    struct, energy, dotbracket = rnafold(design_region)
    print(f"🧬 Switch without toxin:{design_region}\n🧪 Dot-bracket: {dotbracket}\n", flush=True)

    toehold_len = len(design_region) - len(conserved_seq) - len(universal_linker)
    conserved_start = toehold_len
    conserved_end = conserved_start + len(conserved_seq)
    linker_start = conserved_end
    linker_end = linker_start + len(universal_linker)

    conserved_seq_actual = design_region[conserved_start:conserved_end]
    linker_seq_actual = design_region[linker_start:linker_end]

    print(f"🧬 Conserved region:{conserved_seq_actual}\n🧪 Dot-bracket: {dotbracket[conserved_start:conserved_end]}", flush=True)
    print(f"🧬 Universal linker region:{linker_seq_actual}\n🧪 Dot-bracket: {dotbracket[linker_start:linker_end]}", flush=True)

    if not is_region_unpaired(dotbracket, conserved_start, len(conserved_seq)):
        print("📛 Failure reason: Conserved region is not unpaired.\n")
        return False
    if not is_region_unpaired(dotbracket, linker_start, len(universal_linker)):
        print("📛 Failure reason: Universal linker region is not unpaired.\n")
        return False
    return True


def design_seriesB_toehold(toxin_seq, trigger_length=25):
    attempts = 0
    output_path = "/mnt/c/Users/tsdn/Desktop/Igem/igem_toehold/toehold_results.txt"

    while True:
        if attempts % 1000 == 0 and attempts > 0:
            print(f"📣 Attempted {attempts} trigger+switch designs...", flush=True)

        best = trigger_candidate_fast(trigger_length, "GGAAGGAGGUAACAAUG")
        attempts += 1
        if best is None:
            continue

        best_seq, best_score = best[1], best[0]
        print(f"[Attempt {attempts}] Best sequence: {best_seq}, Score: {best_score:.2f}", flush=True)

        struct_trigger, _, dot_trigger = rnafold(best_seq)
        print(f"🧬 Trigger: {best_seq}\n🧪 Dot-bracket: {dot_trigger}", flush=True)

        rev_comp = reverse_complement_rna(best_seq)
        struct_rc, _, dot_rc = rnafold(rev_comp)
        print(f"🧬 Reverse Complement: {rev_comp}\n🧪 Dot-bracket: {dot_rc}", flush=True)

        if sum(1 for c in struct_rc if c != '.') > int(len(rev_comp)):
            print("❌ Reverse complement is too structured (not linear)", flush=True)
            print("📛 Failure reason: Reverse complement structure not linear.\n")
            continue

        toehold_seq = generate_toehold(best_seq, toxin_seq)
        if not validate_toehold_structure(toehold_seq, len(toxin_seq)):
            print("❌ Toehold structure invalid", flush=True)
            continue
        if not validate_toxin(toxin_seq, toehold_seq):
            print("❌ Toxin structure or interaction invalid", flush=True)
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