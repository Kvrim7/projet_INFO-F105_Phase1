import sys
from typing import Tuple, Optional, TextIO

# Définition de l'alphabet autorisé
ALPHABET = {'A', 'C', 'G', 'N', 'T'}
MAX_SEQ_LEN = 1000

def fasta_open_and_read_header(path: str) -> Tuple[Optional[TextIO], Optional[str]]:
    try:
        fp = open(path, "r")
    except OSError:
        return None, None
    header = fp.readline()
    if header == "" or not header.startswith(">"):
        fp.close()
        return None, None
    return fp, header.strip()[1:]

def fasta_read_bases(fp: Optional[TextIO], n: int) -> Optional[str]:
    if fp is None: return None
    bases = []
    while len(bases) < n:
        c = fp.read(1)
        if c == "": break
        if c in ("\n", "\r"): continue
        if c == ">": return None
        bases.append(c)
    return "".join(bases)

def csv_open(path: str) -> Optional[TextIO]:
    try: return open(path, "w")
    except OSError: return None

def csv_write_mutation(fp: TextIO, name: str, pos: int, ref: str, alt: str) -> None:
    fp.write(f"{name},{pos},{ref},{alt}\n")

def csv_close(fp: Optional[TextIO]) -> None:
    if fp is not None: fp.close()

def main():
    # 1. Vérification des arguments
    if len(sys.argv) != 4:
        print("Usage: python compare_strings.py <ref> <seq> <out>")
        return 1

    # 2. Ouverture des fichiers
    fp_ref, name_ref = fasta_open_and_read_header(sys.argv[1])
    fp_seq, name_seq = fasta_open_and_read_header(sys.argv[2])
    
    if not fp_ref or not fp_seq:
        print("Erreur: Impossible d'ouvrir ou lire les fichiers FASTA.")
        return 1

    # 3. Vérifications (Noms et longueurs)
    if name_ref != name_seq:
        print("Erreur: Les noms des séquences diffèrent.")
        return 1

    seq_ref = fasta_read_bases(fp_ref, MAX_SEQ_LEN)
    seq_comp = fasta_read_bases(fp_seq, MAX_SEQ_LEN)

    if len(seq_ref) != len(seq_comp):
        print("Erreur: Les séquences n'ont pas la même longueur.")
        return 1

    # 4. Ouverture du CSV
    fp_out = csv_open(sys.argv[3])
    if not fp_out:
        print("Erreur: Impossible d'ouvrir le fichier de sortie.")
        return 1

    # 5. Comparaison et écriture
    for i in range(len(seq_ref)):
        if seq_ref[i] not in ALPHABET or seq_comp[i] not in ALPHABET:
            print("Erreur: Caractère invalide dans la séquence.")
            return 1
        if seq_ref[i] != seq_comp[i]:
            csv_write_mutation(fp_out, name_ref, i, seq_ref[i], seq_comp[i])

    # 6. Nettoyage
    csv_close(fp_out)
    fp_ref.close()
    fp_seq.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())