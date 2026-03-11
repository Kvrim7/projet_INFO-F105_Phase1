#include <iostream>
#include <cstring>
#include "utils_fasta.h"
#include "utils_csv.h"

// Définitions selon l'énoncé
constexpr size_t MAX_NAME_LEN = 100;
constexpr size_t MAX_SEQ_LEN = 1000;

bool is_valid_alphabet(char c) {
    return (c == 'A' || c == 'C' || c == 'G' || c == 'N' || c == 'T');
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: ./compare_strings <ref_fasta> <seq_fasta> <output_csv>" << std::endl;
        return 1;
    }

    FILE *fp_ref = nullptr, *fp_seq = nullptr, *fp_out = nullptr;
    char name_ref[MAX_NAME_LEN + 1], name_seq[MAX_NAME_LEN + 1];
    char seq_ref[MAX_SEQ_LEN + 1], seq_comp[MAX_SEQ_LEN + 1];
    size_t n_ref = 0, n_comp = 0;

    // --- 1. Ouverture et lecture FASTA ---
    if (fasta_open_and_read_header(argv[1], &fp_ref, name_ref, MAX_NAME_LEN + 1) != FASTA_OK ||
        fasta_open_and_read_header(argv[2], &fp_seq, name_seq, MAX_NAME_LEN + 1) != FASTA_OK) {
        std::cerr << "Erreur lors de l'ouverture des fichiers FASTA." << std::endl;
        return 1;
    }

    if (fasta_read_bases(fp_ref, seq_ref, MAX_SEQ_LEN, &n_ref) != FASTA_OK ||
        fasta_read_bases(fp_seq, seq_comp, MAX_SEQ_LEN, &n_comp) != FASTA_OK) {
        std::cerr << "Erreur lors de la lecture des séquences." << std::endl;
        return 1;
    }
    seq_ref[n_ref] = '\0'; seq_comp[n_comp] = '\0';

    // --- 2. Vérifications ---
    if (std::strcmp(name_ref, name_seq) != 0) {
        std::cerr << "Erreur: Les noms des séquences diffèrent." << std::endl;
        return 1;
    }
    if (n_ref != n_comp) {
        std::cerr << "Erreur: Les longueurs des séquences diffèrent." << std::endl;
        return 1;
    }

    // --- 3. Comparaison et écriture CSV ---
    if (csv_open(argv[3], &fp_out) != CSV_OK) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier de sortie." << std::endl;
        return 1;
    }

    for (size_t i = 0; i < n_ref; ++i) {
        if (!is_valid_alphabet(seq_ref[i]) || !is_valid_alphabet(seq_comp[i])) {
            std::cerr << "Erreur: Caractère non autorisé détecté." << std::endl;
            return 1;
        }
        if (seq_ref[i] != seq_comp[i]) {
            if (csv_write_mutation(fp_out, name_ref, i, seq_ref[i], seq_comp[i]) != CSV_OK) {
                std::cerr << "Erreur lors de l'écriture CSV." << std::endl;
                return 1;
            }
        }
    }

    // --- 4. Fermeture ---
    fasta_close(fp_ref); fasta_close(fp_seq);
    csv_close(fp_out);
    return 0;
}