import subprocess
import random
from Bio.Seq import Seq
from trigger_mrna_generator import get_top_trigger_candidates

# Run RNAfold to compute structure; if name_prefix is given, output .ps file for visualization
def rnafold(seq, name_prefix=None):
    if name_prefix:
        # Run RNAfold and capture both stdout and .ps file output
        cmd = ['RNAfold']
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(seq)
        with open(f"{name_prefix}.fa", "w") as f:
            f.write(f">{name_prefix}\n{seq}\n")
    else:
        # Run without plotting
        cmd = ['RNAfold', '--noPS']
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate(seq)

    if process.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {stderr}")

    lines = stdout.strip().split('\n')
    structure_line = lines[1]
    struct = structure_line.split()[0]
    energy = structure_line.split('(')[-1].strip(')')
    return struct, float(energy), struct

# Convert RNA sequence to its reverse complement (A<->U, C<->G)
def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

# Check if a region is sufficiently unpaired in dot-bracket notation
def is_region_unpaired(dot_bracket, start, length, threshold=0.8):
    region = dot_bracket[start:start+length]
    return region.count('.') / length >= threshold

# Ensure no base pairing occurs between two RNA sequences using RNAcofold
def no_interaction_rnacofold(seq1, seq2):
    input_str = f"{seq1}&{seq2}"
    process = subprocess.Popen(['RNAcofold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input_str)
    if process.returncode != 0:
        raise RuntimeError(f"RNAcofold failed: {stderr}")
    lines = stdout.strip().split('\n')
    struct_line = lines[1]
    struct = struct_line.split()[0]
    amp_index = struct.find('&')

    # Detect if any base pair crosses the boundary between seq1 & seq2
    pair_stack = []
    pairs = {}
    for idx, c in enumerate(struct):
        if c == '(':
            pair_stack.append(idx)
        elif c == ')':
            if not pair_stack:
                raise ValueError("Unbalanced RNAcofold structure")
            start_idx = pair_stack.pop()
            pairs[start_idx] = idx
            pairs[idx] = start_idx

    for i in range(amp_index):
        paired_idx = pairs.get(i)
        if paired_idx is not None and paired_idx > amp_index:
            return False
    return True

# Generate a random RNA sequence of length n
def random_rna_seq(n):
    return ''.join(random.choices('ACGU', k=n))

# Validate folding of toxin and absence of interaction with toehold
def validate_toxin(toxin_seq, toehold_seq):
    struct, energy, _ = rnafold(toxin_seq)
    paired_nt = sum(1 for c in struct if c != '.')
    max_paired = int(0 * len(toxin_seq))  # Allow up to 0 paired
    if paired_nt > max_paired:
        return False
    if not no_interaction_rnacofold(toehold_seq, toxin_seq):
        return False
    return True

# Assemble full toehold switch sequence
def generate_toehold(trigger_mrna, toxin_seq):
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    toehold_trigger = reverse_complement_rna(trigger_mrna)
    return toehold_trigger + conserved_seq + universal_linker + toxin_seq

# Check if conserved_seq and linker are sufficiently unpaired
def validate_toehold_structure(full_seq, toxin_sequence):
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"
    struct, energy, dotbracket = rnafold(full_seq)
    toehold_len = len(full_seq) - len(conserved_seq) - len(universal_linker) - len(toxin_sequence)
    conserved_start = toehold_len
    conserved_end = conserved_start + 25
    universal_linker_start = conserved_end + 11
    if not is_region_unpaired(dotbracket, conserved_start, 25, threshold=0.8):
        return False
    if not is_region_unpaired(dotbracket, universal_linker_start, len(universal_linker), threshold=0.8):
        return False
    return True

# Validate trigger: GC content and folding energy
def validate_trigger_mrna(seq, min_gc=0.3, max_hairpin_energy=-8.0):
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    if gc_content < min_gc:
        return False
    struct, energy, _ = rnafold(seq)
    if energy < max_hairpin_energy:
        return False
    return True

# Repeatedly try to design valid constructs until `required` designs are found
def design_using_trigger_candidates(toxin_sequence, trigger_len=36, num_required=1):
    # Construct full switch template without trigger
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
    universal_linker = "AACCUGGCGGCAGCGCAAAAG"
    switch_backbone = conserved_seq + universal_linker + toxin_sequence

    print(f"🔍 Calling trigger_mrna_generator to get trigger candidates...")
    trigger_candidates = get_top_trigger_candidates(trigger_len, switch_backbone, topk=5)

    if not trigger_candidates:
        print("❌ No trigger candidates found.")
        return []

    results = []
    rejected = 0

    for i, trigger in enumerate(trigger_candidates, 1):
        if i % 500 == 0:
            print(f"📊 Tried {i} candidates so far, {len(results)} valid designs found, {rejected} rejected.")

        toehold_seq = generate_toehold(trigger, toxin_sequence)

        # Step 1: Validate toehold structure
        if not validate_toehold_structure(toehold_seq, toxin_sequence):
            rejected += 1
            continue

        # Step 2: Validate toxin folding & interaction
        if not validate_toxin(toxin_sequence, toehold_seq):
            rejected += 1
            continue

        # ✅ Success
        design_id = f"design_{len(results)+1}"
        rnafold(toehold_seq, name_prefix=design_id)

        print(f"\n✅ --- Toehold Design {len(results)+1} ---")
        print(f"Trigger mRNA (5'→3'): {trigger}")
        print(f"Toehold Trigger (5'→3'): {reverse_complement_rna(trigger)}")
        print(f"Saved: {design_id}.fa + .ps\n")

        results.append({
            'trigger_mRNA': trigger,
            'toehold_sequence': toehold_seq,
        })

        if len(results) >= num_required:
            break

    print(f"\n🎯 Finished: {len(results)} designs found from {i} candidates, {rejected} rejected.\n")
    return results


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Design Series B toehold switches with toxin gene input.")
    parser.add_argument("--toxin", type=str, required=True, help="Toxin gene RNA sequence (5' to 3').")
    parser.add_argument("--trigger_length", type=int, default=25, help="Length of trigger mRNA sequence.")
    args = parser.parse_args()

    toxin_sequence = args.toxin.upper().replace('T', 'U').replace(' ', '').replace('\n','')
    trigger_len = args.trigger_length
    switch_backbone = "GGACUUUAGAACAGAGGAGAUAAAGAUG" + "AACCUGGCGGCAGCGCAAAAG" + toxin_sequence

    print(f"Designing toehold switches for toxin length {len(toxin_sequence)} with trigger length {trigger_len}...")

    designs = design_using_trigger_candidates(toxin_sequence, trigger_len=trigger_len, num_required=2)

    if not designs:
        print("❌ No valid toehold switches found.")
    else:
        for i, res in enumerate(designs):
            print(f"\n✅ --- Toehold Design {i+1} ---")
            print(f"Trigger mRNA (5'→3'): {res['trigger_mRNA']}")
            print(f"Toehold Trigger (5'→3'): {reverse_complement_rna(res['trigger_mRNA'])}")
            print(f"Full Toehold + Toxin:\n{res['toehold_sequence']}")
