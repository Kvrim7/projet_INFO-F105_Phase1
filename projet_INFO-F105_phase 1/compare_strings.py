import sys
from typing import Tuple, Optional, TextIO

ALPHABET = {'A', 'C', 'G', 'N', 'T'}

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
    bases =[]
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

def print_error(msg: str):
    """Affiche un message d'erreur sur la sortie standard d'erreur."""
    print(f"Erreur: {msg}", file=sys.stderr)

def compare_blocks(buf_ref: str, buf_seq: str, global_pos: int, fp_out: TextIO, seq_name: str) -> bool:
    """Compare deux blocs et écrit les mutations avec leur position globale."""
    for i in range(len(buf_ref)):
        if buf_ref[i] not in ALPHABET or buf_seq[i] not in ALPHABET:
            print_error("Caractère non autorisé détecté.")
            return False
        if buf_ref[i] != buf_seq[i]:
            # La position réelle = position de départ du bloc + index dans le bloc
            csv_write_mutation(fp_out, seq_name, global_pos + i, buf_ref[i], buf_seq[i])
    return True

def process_sequences(path_ref: str, path_seq: str, path_out: str, buffer_size: int) -> bool:
    """Gère l'ouverture des fichiers et la boucle de lecture par blocs."""
    fp_ref, name_ref = fasta_open_and_read_header(path_ref)
    fp_seq, name_seq = fasta_open_and_read_header(path_seq)

    if not fp_ref or not fp_seq:
        print_error("Impossible d'ouvrir ou de lire les fichiers FASTA.")
        if fp_ref: fp_ref.close()
        if fp_seq: fp_seq.close()
        return False

    if name_ref != name_seq:
        print_error("Les noms des séquences diffèrent.")
        fp_ref.close(); fp_seq.close()
        return False

    fp_out = csv_open(path_out)
    if not fp_out:
        print_error("Impossible de créer le fichier CSV.")
        fp_ref.close(); fp_seq.close()
        return False

    global_pos = 0

    # Boucle infinie qui lit bloc par bloc
    while True:
        buf_ref = fasta_read_bases(fp_ref, buffer_size)
        buf_seq = fasta_read_bases(fp_seq, buffer_size)

        # Si les deux buffers sont vides, on a fini de lire les fichiers avec succès
        if not buf_ref and not buf_seq:
            break

        # Si l'un renvoie None (Erreur de format ou second > trouvé)
        if buf_ref is None or buf_seq is None:
            print_error("Erreur de format FASTA durant la lecture.")
            csv_close(fp_out); fp_ref.close(); fp_seq.close()
            return False

        # Vérification si un fichier est plus long que l'autre
        if len(buf_ref) != len(buf_seq):
            print_error("Les séquences n'ont pas la même longueur.")
            csv_close(fp_out); fp_ref.close(); fp_seq.close()
            return False

        # Comparaison du bloc courant
        if not compare_blocks(buf_ref, buf_seq, global_pos, fp_out, name_ref):
            csv_close(fp_out); fp_ref.close(); fp_seq.close()
            return False

        # On avance la position globale pour le prochain bloc
        global_pos += len(buf_ref)

    # Nettoyage propre
    csv_close(fp_out)
    fp_ref.close()
    fp_seq.close()
    return True

def main():
    # 4 arguments attendus au lieu de 3 (ajout de buffer_size)
    if len(sys.argv) != 5:
        print_error("Usage: python compare_strings.py <ref> <seq> <out> <buffer_size>")
        return 1

    try:
        buffer_size = int(sys.argv[4])
        if buffer_size <= 0: raise ValueError
    except ValueError:
        print_error("Le paramètre buffer_size doit être un entier positif.")
        return 1

    # Lancement du processus
    if not process_sequences(sys.argv[1], sys.argv[2], sys.argv[3], buffer_size):
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())