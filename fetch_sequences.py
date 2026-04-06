from Bio import Entrez, SeqIO
import os, ssl

ssl._create_default_https_context = ssl._create_unverified_context

# Аккессии белковых последовательностей (RefSeq):
# PAH  человек : NP_000268.1   — фенилаланин-4-гидроксилаза, 452 aa
# PAH  мышь   : NP_032803.2   — phenylalanine-4-hydroxylase, 453 aa
# QDPR человек : NP_000311.2   — дигидроптеридинредуктаза изоформа 1, 244 aa
# QDPR мышь   : NP_077198.1   — dihydropteridine reductase, 241 aa

SEQUENCES = {
    "PAH_human":  "NP_000268.1",
    "PAH_mouse":  "NP_032803.2",
    "QDPR_human": "NP_000311.2",
    "QDPR_mouse": "NP_077198.1",
}

OUT_DIR = "sequences"
os.makedirs(OUT_DIR, exist_ok=True)

for name, acc in SEQUENCES.items():
    out_path = os.path.join(OUT_DIR, f"{name}.fasta")
    print(f"Скачиваю {name} ({acc}) ...", end=" ", flush=True)
    handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    record.id = name
    record.name = name
    with open(out_path, "w") as f:
        SeqIO.write(record, f, "fasta")
    print(f"OK  ({len(record.seq)} aa)  →  {out_path}")

print("\nВсе последовательности сохранены.")
