import os
import multiprocessing
from design_toehold import design_seriesB_toehold, reverse_complement_rna

def run_parallel_design(toxin_rna: str, trigger_len: int = 25, output_dir: str = "./batch_results") -> None:
    toxin_clean = toxin_rna.replace("T", "U").replace(" ", "").replace("\n", "").upper()
    os.makedirs(output_dir, exist_ok=True)

    result = design_seriesB_toehold(toxin_clean, trigger_length=trigger_len)
    if not result:
        print(f"❌ No valid toehold found for toxin: {toxin_clean[:30]}...")
        return

    res = result[0]
    trigger = res['trigger_mRNA']
    comp = reverse_complement_rna(trigger)
    switch = res['toehold_sequence']
    score = res['score']

    # ✅ 打印 switch 详细信息
    print("\n🧪 Toehold Switch Found:")
    print(f"Trigger mRNA (5'->3'): {trigger}")
    print(f"Complementary Trigger (5'->3'): {comp}")
    print(f"Score: {score:.2f}")
    print(f"Full Switch:\n{switch[:100]}... [length={len(switch)} nt]\n")

    # 输出文件名：trigger + pid，避免重复
    filename = os.path.join(output_dir, f"toxin_{trigger}_{os.getpid()}.txt")
    with open(filename, "w") as f:
        f.write(f"Trigger mRNA (5'->3'): {trigger}\n")
        f.write(f"Complementary Toehold Trigger (5'->3'): {comp}\n")
        f.write(f"Score: {score:.2f}\n")
        f.write(f"Full Toehold Switch + Toxin sequence:\n{switch}\n")

    print(f"✅ Saved result for toxin starting with {toxin_clean[:20]}...\n")

if __name__ == "__main__":
    toxin_list = [
        "ATGGTAAGCCGATACGTACCCGATATGGGCGATCTGATTTGGGTTGATTTTGACCCGACAAAAGGTAGCGAGCAAGCTGGACATCGTCCAGCTGTTGTCCTGAGTCCTTTCATGTACAACAACAAAACAGGTATGTGTCTGTGTGTTCCTTGTACAACGCAATCAAAAGGATATCCGTTCGAAGTTGTTTTATCCGGTCAGGAACGTGATGGCGTAGCGTTAGCTGATCAGGTAAAAAGTATCGCCTGGCGGGCAAGAGGAGCAACGAAGAAAGGAACAGTTGCCCCAGAGGAATTACAACTCATTAAAGCCAAAATTAACGTACTGATTGGGTAG",
    ]

    os.makedirs("./batch_results", exist_ok=True)

    with multiprocessing.Pool(processes=min(2, len(toxin_list))) as pool:
        pool.starmap(run_parallel_design, [(seq, 25, "./batch_results") for seq in toxin_list])
