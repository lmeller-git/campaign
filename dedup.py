import glob
import os

OUTPUT_DIR = "./sequences/binders"
INPUT_DIRS = ["./sequences4/binders", "./sequences2/binders2", "./sequences3/binders"]

PARENT_SEQUENCE = "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAGGGGSGGGGSGGGGSDILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPRLLIKYASESISGIPSRFSGSGSGTDFTLSINSVESEDIADYYCQQNNNWPTTFGAGTKLELK"


# def apply_mutations(base_seq, mutations_list):
#     seq_arr = list(base_seq)
#     for mut in mutations_list:
#         pos = int(mut[1:-1]) - 1
#         seq_arr[pos] = mut[-1]
#     return "".join(seq_arr)

# CRADLE_MUTATIONS = ["K5Q", "N70S", "K71R", "N73T", "S87T", "N88D", "R179K", "K183R", "E213D", "S214P", ""]
# BIOBETTER_CRADLE = apply_mutations(PARENT_SEQUENCE, CRADLE_MUTATIONS)
# KNOWN_BIOBETTERS = [BIOBETTER_CRADLE, PARENT_SEQUENCE]

BANNED_REFERENCES = {
    "Parent": PARENT_SEQUENCE,

    # Cradle R2 (10 framework mutations)
    "Cradle_R2": (
        "QVQLQQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSISRDTSKSQVFFKMNSLQTDDTAIYYCARALTYYDYEFAYWGQGTLVTVSAGGGGSGGGGSGGGGSDILLTQSPVILSVSPGERVSFSCRASQSIGTNIHWYQQRTNGSPKLLIRYASESISGIPSRFSGSGSGTDFTLSINSVDPEDIADYYCQQNNNWPTTFGAGTKLELK"
    ),
    # Syndra DSM
    "Syndra_DSM":(
        "QVQLQQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSISRDTSKSQVFFKMNSLQTDDTAVYYCARALTYYDYEFAYWGQGTLVTVSAGGGGSGGGGSGGGGSDILLTQSPVILSVSPGERVSFSCRASQSIGSNIHWYQQRTNGSPKLLIRYASESISGIPSRFSGSGSGTDFTLSINSVDPEDIADYYCQQNNNWPTTFGAGTKLEIK"
        ),
    # Converge Biobetter (VH T61A/S87A/N88D + VL V9A/N32D/N93A)
    "Converge_6_Edit": (
        "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNAPFTSRLSINKDNSKSQVFFKMNSLQADDTAIYYCARALTYYDYEFAYWGQGTLVTVSAGGGGSGGGGSGGGGSDILLTQSPAILSVSPGERVSFSCRASQSIGTDIHWYQQRTNGSPRLLIKYASESISGIPSRFSGSGSGTDFTLSINSVESEDIADYYCQQNNAWPTTFGAGTKLELK"
    )
}

def get_mutation_list(mutated_seq, parent_seq=PARENT_SEQUENCE):
    """Identifies the explicit amino acid substitutions."""
    mutations = []
    for idx, (p, m) in enumerate(zip(parent_seq, mutated_seq)):
        if p != m:
            mutations.append(f"{p}{idx+1}{m}")
    return mutations

def compute_hamming(seq1, seq2):
    """Calculates edit distance between two equal-length strings."""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def is_compliant_sequence(candidate_seq, parent_seq, banned_references, min_distance=5):
    """
    Enforces that a candidate sequence is at least 'min_distance' mutations
    away from both the wildtype parent sequence and all existing patent/banned biobetters.
    """
    # 1. Check distance from the baseline parent sequence
    if compute_hamming(candidate_seq, parent_seq) < min_distance:
        return False

    # 2. Check distance from all patent/banned references
    for ref_name, ref_seq in banned_references.items():
        if compute_hamming(candidate_seq, ref_seq) < min_distance:
            return False

    return True

def main():
    seqs = set()

    # with open("./output_a.txt", "r") as f:
    #     for line in f.readlines():
    #         seqs.add(line)

    for f in INPUT_DIRS:
        fastas = glob.glob(os.path.join(f,"*.fasta"))
        for fasta in fastas:
            with open(fasta, "r") as f:
                _name = f.readline()
                seq = f.readline()
                seqs.add(seq)

    os.mkdir(OUTPUT_DIR)
    for i, seq in enumerate(seqs):
        if not is_compliant_sequence(seq, PARENT_SEQUENCE, BANNED_REFERENCES):
            print(f"seq_{i} is banned")
            continue
        with open(f"{OUTPUT_DIR}/sequence_{i}.fasta", "w") as f:
            _ = f.write(f"> sequence_{i}\n")
            _ = f.write(f"{seq}\n")



if __name__ == "__main__":
    main()
