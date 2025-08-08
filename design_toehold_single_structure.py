import subprocess
from Bio.Seq import Seq

# === Fixed Constants === #
TRIGGER_LEN = 25
CONSERVE_SEQ = "GGACUUUAGAACAGAGGAGAUAAAGAUG"
LINKER_SEQ = "AACCUGGCGGCAGCGCAAAAG"
CONSERVE_LEN = len(CONSERVE_SEQ)
LINKER_LEN = len(LINKER_SEQ)
LOOP_LEN = 11

# === Utils === #
def reverse_complement_rna(seq):
    return str(Seq(seq).reverse_complement()).replace('T', 'U')

def rnafold(seq):
    result = subprocess.run(['RNAfold', '--noPS'],
                            input=seq,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
    if result.returncode != 0:
        raise RuntimeError(f"RNAfold failed: {result.stderr}")
    lines = result.stdout.strip().split('\n')
    structure = lines[1].split()[0]
    energy = float(lines[1].split('(')[-1].strip(')'))
    return structure, energy

def rnainverse(target_structure, max_attempts=100):
    """Keep retrying RNAinverse until structure matches"""
    for attempt in range(max_attempts):
        result = subprocess.run(['RNAinverse'],
                                input=target_structure,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True)
        if result.returncode != 0:
            raise RuntimeError(f"RNAinverse failed: {result.stderr}")
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if not line or line.startswith('%'): continue
            parts = line.strip().split()
            if len(parts) < 2: continue
            seq, struct = parts[0], parts[1]
            struct_check, energy = rnafold(seq)
            if struct_check == target_structure:
                print(f"[INFO] Found valid sequence at attempt {attempt+1}")
                return seq, struct_check, energy
    return None, None, None

def build_target_structure():
    """Builds fixed target structure: trigger . + hairpin + linker"""
    return (
        '.' * TRIGGER_LEN +
        '(((((((((((.................................))))))))))).' +
        '.' * LINKER_LEN
    )

# === Main Design === #
def design_toehold_fixed(trigger_seq, toxin_seq):
    print(f"[INFO] Designing toehold for trigger: {trigger_seq}")
    rc_trigger = reverse_complement_rna(trigger_seq)
    print(f"[INFO] Reverse complement: {rc_trigger}")

    target_struct = build_target_structure()
    print(f"[INFO] Target structure:\n{target_struct}")

    candidate, struct, energy = rnainverse(target_struct)
    if not candidate:
        print("[FAIL] Could not generate valid sequence after retries.")
        return

    print(f"[PASS] Candidate sequence:\n{candidate}")
    print(f"[PASS] Structure: {struct}")
    print(f"[PASS] ΔG = {energy:.2f} kcal/mol")

    # Validate subsequences
    conserve_start = TRIGGER_LEN + LOOP_LEN
    linker_start = conserve_start + CONSERVE_LEN + 1

    actual_conserve = candidate[conserve_start:conserve_start + CONSERVE_LEN]
    actual_linker = candidate[linker_start:linker_start + LINKER_LEN]

    if actual_conserve != CONSERVE_SEQ:
        print(f"[FAIL] Conserved sequence mismatch:\nExpected: {CONSERVE_SEQ}\nActual:   {actual_conserve}")
        return
    if actual_linker != LINKER_SEQ:
        print(f"[FAIL] Linker sequence mismatch:\nExpected: {LINKER_SEQ}\nActual:   {actual_linker}")
        return

    full_seq = candidate + toxin_seq
    print(f"[✅ FINAL] Full switch + toxin:\n{full_seq}")

# === CLI Entry === #
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Design a toehold switch using RNAinverse and fixed structure.")
    parser.add_argument("--trigger", required=True, help="Trigger RNA sequence (5'->3')")
    parser.add_argument("--toxin", required=True, help="Toxin gene RNA sequence (5'->3')")
    args = parser.parse_args()

    trigger = args.trigger.upper().replace('T', 'U').replace(' ', '')
    toxin = args.toxin.upper().replace('T', 'U').replace(' ', '')

    design_toehold_fixed(trigger, toxin)
