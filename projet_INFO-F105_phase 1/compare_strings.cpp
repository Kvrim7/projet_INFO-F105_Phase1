#include <iostream>
#include <cstring>
#include <cstdlib>
#include "utils_fasta.h"
#include "utils_csv.h"

constexpr size_t MAX_NAME_LEN = 256;

// 1. Fonction pour afficher les erreurs
void print_error(const char* msg) {
    std::cerr << "Erreur : " << msg << std::endl;
}

// 2. Fonction pour nettoyer la mémoire et fermer les fichiers
void cleanup_resources(FILE* fp_ref, FILE* fp_seq, FILE* fp_out, char* buf_ref, char* buf_seq) {
    fasta_close(fp_ref);
    fasta_close(fp_seq);
    csv_close(fp_out);
    delete[] buf_ref;
    delete[] buf_seq;
}

// Fonction de validation de l'alphabet
bool is_valid_alphabet(char c) {
    return (c == 'A' || c == 'C' || c == 'G' || c == 'N' || c == 'T');
}

// 3. Fonction pour comparer deux blocs et écrire les mutations
bool compare_blocks(const char* buf_ref, const char* buf_seq, size_t n_read, size_t global_pos, FILE* fp_out, const char* seq_name) {
    for (size_t i = 0; i < n_read; ++i) {
        if (!is_valid_alphabet(buf_ref[i]) || !is_valid_alphabet(buf_seq[i])) {
            print_error("Caractère non autorisé détecté dans la séquence.");
            return false;
        }
        if (buf_ref[i] != buf_seq[i]) {
            // Attention : on utilise (global_pos + i) pour la position réelle dans le génome complet
            if (csv_write_mutation(fp_out, seq_name, global_pos + i, buf_ref[i], buf_seq[i]) != CSV_OK) {
                print_error("Problème lors de l'écriture dans le fichier CSV.");
                return false;
            }
        }
    }
    return true;
}

// 4. Fonction contenant la boucle principale de lecture par blocs
bool process_sequences(const char* ref_path, const char* seq_path, const char* out_path, size_t buffer_size) {
    FILE *fp_ref = nullptr, *fp_seq = nullptr, *fp_out = nullptr;
    char *buf_ref = nullptr, *buf_seq = nullptr;
    char name_ref[MAX_NAME_LEN], name_seq[MAX_NAME_LEN];

    // Allocation dynamique sur le TAS (Heap)
    buf_ref = new char[buffer_size];
    buf_seq = new char[buffer_size];

    // Ouverture des fichiers FASTA
    if (fasta_open_and_read_header(ref_path, &fp_ref, name_ref, MAX_NAME_LEN) != FASTA_OK ||
        fasta_open_and_read_header(seq_path, &fp_seq, name_seq, MAX_NAME_LEN) != FASTA_OK) {
        print_error("Impossible d'ouvrir les fichiers FASTA.");
        cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
        return false;
    }

    // Vérification des noms
    if (std::strcmp(name_ref, name_seq) != 0) {
        print_error("Les noms des séquences diffèrent.");
        cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
        return false;
    }

    // Ouverture du fichier CSV de sortie
    if (csv_open(out_path, &fp_out) != CSV_OK) {
        print_error("Impossible de créer le fichier de sortie CSV.");
        cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
        return false;
    }

    size_t global_pos = 0;
    size_t n_ref = 0, n_seq = 0;

    // Boucle principale de lecture par blocs
    while (true) {
        if (fasta_read_bases(fp_ref, buf_ref, buffer_size, &n_ref) != FASTA_OK ||
            fasta_read_bases(fp_seq, buf_seq, buffer_size, &n_seq) != FASTA_OK) {
            print_error("Erreur de lecture dans les fichiers FASTA.");
            cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
            return false;
        }

        // Si on a atteint la fin des deux fichiers en même temps
        if (n_ref == 0 && n_seq == 0) {
            break; 
        }

        // Vérification que les blocs ont la même taille
        if (n_ref != n_seq) {
            print_error("Les deux séquences n'ont pas la même longueur.");
            cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
            return false;
        }

        // Comparaison du bloc actuel
        if (!compare_blocks(buf_ref, buf_seq, n_ref, global_pos, fp_out, name_ref)) {
            cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
            return false;
        }

        // Mise à jour de la position globale
        global_pos += n_ref;
    }

    // Libération de la mémoire et fermeture
    cleanup_resources(fp_ref, fp_seq, fp_out, buf_ref, buf_seq);
    return true;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: ./compare_strings <ref> <seq> <out> <buffer_size>\n";
        return 1;
    }

    // Conversion du 4ème argument en entier
    size_t buffer_size = std::strtoull(argv[4], nullptr, 10);
    if (buffer_size == 0) {
        print_error("La taille du buffer doit être supérieure à 0.");
        return 1;
    }

    // Lancement du traitement
    if (!process_sequences(argv[1], argv[2], argv[3], buffer_size)) {
        return 1;
    }

    return 0;
}