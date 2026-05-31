import glob
import os

OUTPUT_DIR = "./sequences/binders"
INPUT_DIRS = ["./sequences4/binders", "./sequences2/binders2", "./sequences3/binders"]


def main():
    seqs = set()

    for f in INPUT_DIRS:
        fastas = glob.glob(os.path.join(f,"*.fasta"))
        for fasta in fastas:
            with open(fasta, "r") as f:
                _name = f.readline()
                seq = f.readline()
                seqs.add(seq)

    os.mkdir(OUTPUT_DIR)
    for i, seq in enumerate(seqs):
        with open(f"{OUTPUT_DIR}/sequence_{i}.fasta", "w") as f:
            _ =f.write(f"> sequence_{i}\n")
            _ = f.write(f"{seq}\n")



if __name__ == "__main__":
    main()
