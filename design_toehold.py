import subprocess
import random
from Bio.Seq import Seq

# External command call to RNAfold, assumes RNAfold is installed and accessible
def rnafold(seq):
    """
    Runs RNAfold on the given RNA sequence.
    Returns tuple: (mfe_structure, free_energy, dot_bracket)
    """
    process = subprocess.Popen(['RNAfold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(seq)
   
    if process.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {stderr}")
   
    lines = stdout.strip().split('\n')
    # RNAfold output: line1=seq, line2=structure (energy)
    structure_line = lines[1]
    struct = structure_line.split()[0]
    energy = structure_line.split('(')[-1].strip(')')
    return struct, float(energy), struct

# Reverse complement for RNA using Biopython, output RNA bases
def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

# Check that all positions in a region are unpaired (dot) in dot-bracket notation
def is_region_unpaired(dot_bracket, start, length):
    region = dot_bracket[start:start+length]
    return all(c == '.' for c in region)

# Check no base pairing between two sequences using RNAcofold (ViennaRNA)
def no_interaction_rnacofold(seq1, seq2):
    """
    Returns True if no base pairing is predicted between seq1 and seq2 when cofolded.
    Uses RNAcofold from ViennaRNA.
    RNAcofold input: "seq1&seq2"
    Output structure notation uses '&' marking boundary.
    """
    input_str = f"{seq1}&{seq2}"
    process = subprocess.Popen(['RNAcofold', '--noPS'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(input_str)
    if process.returncode != 0:
        raise RuntimeError(f"RNAcofold failed: {stderr}")
    lines = stdout.strip().split('\n')
    struct_line = lines[1]
    struct = struct_line.split()[0]

    # Check for pairs crossing the '&' delimiter (position of '&' is len(seq1))
    # Paired bases are '(' or ')', unpaired '.'
    # If any '(' or ')' after seq1 length corresponds to pairing between two strands, it's interaction.
    amp_index = struct.find('&')
    # Left and right parts structure
    left_struct = struct[:amp_index]
    right_struct = struct[amp_index+1:]

    # RNAcofold uses parentheses to represent base pairs; pairs between strands will form
    # pairs crossing the '&' symbol. To detect interactions, check for any base pair
    # spanning the boundary.
    # We'll parse the dot bracket and see if any pair crosses between the sequences.
    # We'll use stack method to determine pairs and if pairs cross boundary.

    pair_stack = []
    pairs = {}

    for idx, c in enumerate(struct):
        if c == '(':
            pair_stack.append(idx)
        elif c == ')':
            if not pair_stack:
                raise ValueError("Unbalanced RNAcofold pairing notation")
            start_idx = pair_stack.pop()
            pairs[start_idx] = idx
            pairs[idx] = start_idx

    # Check if any pairing crosses the & position boundary
    for i in range(amp_index):
        paired_idx = pairs.get(i)
        if paired_idx is not None and paired_idx > amp_index:
            # base pair crosses strands => interaction
            return False

    # No cross strand base pairing found
    return True

# Generate random RNA sequence of length n
def random_rna_seq(n):
    return ''.join(random.choices('ACGU', k=n))

# Validate toxin gene - must be mostly unpaired (linear) in folding and no binding to toehold
def validate_toxin(toxin_seq, toehold_seq, check_boundaries=True):
    # 1. Check toxin folding: tolerate some pairing inside toxin gene,
    # but ideally as linear as possible - define threshold of max paired positions
    struct, energy, _ = rnafold(toxin_seq)
    # Count paired nucleotides within toxin
    paired_nt = sum(1 for c in struct if c != '.')
    # Threshold: tolerate up to 30% paired nucleotides - can be parameterized
    max_paired = int(0.3 * len(toxin_seq))

    if paired_nt > max_paired:
        return False

    # 2. Check no interaction between toxin and toehold
    if not no_interaction_rnacofold(toehold_seq, toxin_seq):
        return False

    return True

# Generate toehold switch sequence based on inputs
def generate_toehold(trigger_mrna, toxin_seq):
    """
    Builds Series B toehold switch:
    Format (from user input):
    25.11(GGACUUUAGAACAGAGGAGAUAAAGAUG11).AACCUGGCGGCAGCGCAAAAG + downstream toxin gene
   
    Steps:
    - toehold region (N positions) + fixed parts
    - toehold 5' region complementary (reverse comp) to trigger mRNA
    """
    # Universal linker from prompt
    universal_linker = "AACCUGGCGGCAGCGCAAAAG" # length = 22
   
    # Fixed hairpin conserved sequence from user:
    # "GGACUUUAGAACAGAGGAGAUAAAGAUG" (length 25)
    conserved_seq = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
   
    # Convert trigger mRNA to 5'-3' sequence if not already
    # The toehold trigger domain is complementary of trigger mRNA 5'-3'
    # So toehold 5' end of the switch trigger region:
    toehold_trigger = reverse_complement_rna(trigger_mrna)
   
    # Construct the full toehold sequence:
    # Structure notation says: 25.11( conserved_seq 11 ).universal_linker
    # The 25 corresponds to conserved_seq length
    # The 11 after 25 probably refers to loop length (we keep fixed)
    # The toehold trigger (positions N) goes before conserved sequence as 'N' or complementary region
    # We assume toehold_trigger length = length of trigger_mrna (user defined)
   
    # Place toehold_trigger at 5' end
    # Then conserved sequence (the hairpin forming region)
    # Then universal linker
    # Finally, toxin gene
   
    full_seq = toehold_trigger + conserved_seq + universal_linker + toxin_seq
    return full_seq

# Check that first 25 + universal linker region (defined positions) are unpaired on structure
def validate_toehold_structure(full_seq):
    struct, energy, dotbracket = rnafold(full_seq)
    # First 25 nt = conserved_seq positions (should be unpaired)
    # Universal linker occupies next 22nt after conserved
    # Positions 0-24 and (25+11+toehold length) to (25+11+toehold length + 21)
    # but toehold_trigger length is variable
   
    # We actually want the USER defined toehold trigger length at 5' end not paired
    # The prompt says the first 25 nucleotides of the toehold (5' to 3') AND the universal linker must be unpaired
   
    # The "first 25 nucleotides on the toehold" probably means the first 25 nt of the entire toehold switch RNA sequence (excluding trigger)
    # We do have to clarify positions:
    # Positions in full_seq:
    # 0 to len(toehold_trigger)-1 = toehold trigger region (complementary to trigger mRNA)
    # Then conserved_seq (positions from len(toehold_trigger) to len(toehold_trigger)+24)
    # Then 11nt loop?
    # Then universal linker (positions next 22 nt)
   
    # But the conserved_seq length is 25 nt, and prompt says 25.11(GGAC....11), so after 25 + 11 loop then universal linker.
    # We do not have explicit loop nt sequence here, but structure implies loop of length 11
    # For simplicity, assume no sequence change in loop, loop part is part of conserved sequence or counted separately
    # To approximate, we'll ensure:
    #  - first 25 nt starting at position len(toehold_trigger) are unpaired (conserved sequence)
    #  - universal linker right after that 25+11 nt are unpaired (22 nt universal linker)
   
    # Let's calculate indices:
    toehold_len = len(full_seq) - len(conserved_seq) - len(universal_linker) - len(toxin_sequence)
    conserved_start = toehold_len
    conserved_end = conserved_start + 25
    # Loop of 11 nt after conserved, so universal linker starts at:
    universal_linker_start = conserved_end + 11
    universal_linker_end = universal_linker_start + len(universal_linker)
   
    # Check conserved sequence 25nt unpaired
    if not is_region_unpaired(dotbracket, conserved_start, 25):
        #print(f"Conserved region paired: {dotbracket[conserved_start:conserved_start+25]}")
        return False
   
    # Check universal linker 22nt unpaired
    if not is_region_unpaired(dotbracket, universal_linker_start, len(universal_linker)):
        #print(f"Universal linker region paired: {dotbracket[universal_linker_start:universal_linker_start+len(universal_linker)]}")
        return False
   
    return True

# Validate trigger mRNA (reuse previous criteria):
def validate_trigger_mrna(seq, min_gc=0.4, max_hairpin_energy=-5.0):
    """
    Validates the candidate trigger mRNA sequence.
    - Check GC content (>= min_gc)
    - Folding energy should meet threshold (previous trigger mRNA code)
    We fold with RNAfold and check MFE.
    """
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)
    if gc_content < min_gc:
        return False
   
    struct, energy, _ = rnafold(seq)
    # energy threshold - at least > max_hairpin_energy (meaning less negative energy)
    # to avoid too stable hairpins
    if energy < max_hairpin_energy:
        return False
   
    return True

# Generate trigger candidates from random RNA sequences fulfilling constraints
def generate_trigger_candidates(toxin_seq, num_trials=100, trigger_len=25):
    """
    Generates trigger mRNA candidates following previous design rules.
    Returns list of valid triggers.
    """
    valid_triggers = []
    attempts = 0
    while len(valid_triggers) < num_trials and attempts < num_trials * 10:
        candidate = random_rna_seq(trigger_len)
        if validate_trigger_mrna(candidate):
            valid_triggers.append(candidate)
        attempts += 1
    return valid_triggers

# Main function integrating all checks and generation
def design_seriesB_toehold(toxin_sequence, num_trials=10, trigger_length=25):
    results = []
    triggers = generate_trigger_candidates(toxin_sequence, num_trials=num_trials, trigger_len=trigger_length)

    for trig in triggers:
        # Generate toehold switch sequence
        toehold_seq = generate_toehold(trig, toxin_sequence)

        # Validate toehold secondary structure for unpaired conserved + universal linker regions
        if not validate_toehold_structure(toehold_seq):
            #print(f"Toehold failed structure check for trigger {trig}")
            continue

        # Validate toxin folding and interaction with toehold
        if not validate_toxin(toxin_sequence, toehold_seq):
            #print(f"Toxin interaction/folding invalid for trigger {trig}")
            continue
       
        # All validations passed
        results.append({
            'trigger_mRNA': trig,
            'toehold_sequence': toehold_seq,
        })
    return results


if __name__ == "__main__":
    import argparse
   
    parser = argparse.ArgumentParser(description="Design Series B toehold switches with toxin gene input.")
    parser.add_argument("--toxin", type=str, required=True, help="Toxin gene RNA sequence (5' to 3').")
    parser.add_argument("--trials", type=int, default=10, help="Number of trigger mRNA candidates to generate and test.")
    parser.add_argument("--trigger_length", type=int, default=25, help="Length of trigger mRNA sequence.")
    args = parser.parse_args()

    toxin_sequence = args.toxin.upper().replace('T', 'U').replace(' ', '').replace('\n','')
    num_trials = args.trials
    trigger_len = args.trigger_length

    print(f"Starting design with toxin gene length {len(toxin_sequence)}, {num_trials} trials, trigger length {trigger_len}")

    designs = design_seriesB_toehold(toxin_sequence, num_trials=num_trials, trigger_length=trigger_len)

    if not designs:
        print("No valid toehold switches found with given constraints.")
    else:
        for i, res in enumerate(designs):
            print(f"\n--- Toehold design {i+1} ---")
            print(f"Trigger mRNA (5'->3'): {res['trigger_mRNA']}")
            # The complementary to trigger_mrna 5'->3' for validation display
            comp_trigger = reverse_complement_rna(res['trigger_mRNA'])
            print(f"Complementary Toehold Trigger region (5'->3'): {comp_trigger}")
            print(f"Full Toehold Switch + Toxin seq:\n{res['toehold_sequence']}")