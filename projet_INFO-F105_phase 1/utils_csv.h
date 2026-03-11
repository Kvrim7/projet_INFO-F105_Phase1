#ifndef UTILS_CSV_H
#define UTILS_CSV_H

#include <cstddef> // size_t
#include <cstdio>  // FILE

/*
 * Codes de retour des fonctions utilitaires pour l’écriture du fichier CSV.
 *
 * Comme pour les fonctions FASTA :
 *  - ces fonctions n’affichent aucun message d’erreur,
 *  - elles ne terminent jamais le programme.
 *
 * Il est de la responsabilité de votre programme :
 *  - de vérifier les valeurs de retour,
 *  - d’afficher un message d’erreur explicite en cas de problème,
 *  - de terminer proprement le programme si nécessaire.
 */
enum CsvStatus {
    CSV_OK = 0,       // Écriture réussie
    CSV_ERR_OPEN = 1, // Impossible d’ouvrir le fichier de sortie
    CSV_ERR_IO = 2    // Erreur d’entrée/sortie lors de l’écriture
};

/*
 * Ouvre un fichier CSV en écriture.
 *
 * Paramètres :
 *  - path : chemin vers le fichier CSV à créer ou écraser.
 *  - out_fp : adresse d’une variable FILE*.
 *             En cas de succès, *out_fp contient un FILE* valide ouvert en écriture.
 *
 * Valeur de retour :
 *  - CSV_OK en cas de succès.
 *  - CSV_ERR_OPEN si le fichier ne peut pas être ouvert.
 */
CsvStatus csv_open(const char* path, FILE** out_fp);

/*
 * Écrit une ligne dans le fichier CSV correspondant à une substitution détectée.
 *
 * Format de la ligne écrite :
 *   nom_sequence,position,caractere_reference,caractere_sequence
 *
 * Paramètres :
 *  - out_fp : FILE* précédemment obtenu via csv_open.
 *  - name : nom de la séquence (issu de l’en-tête FASTA).
 *  - pos : position de la substitution dans la séquence (index commençant à 0).
 *  - ref : caractère dans la séquence de référence.
 *  - alt : caractère dans la séquence à comparer.
 *
 * Exemple de ligne écrite :
 *   chr1,42,A,T
 *
 * Valeur de retour :
 *  - CSV_OK si l’écriture s’est déroulée correctement.
 *  - CSV_ERR_IO en cas d’erreur d’entrée/sortie.
 */
CsvStatus csv_write_mutation(FILE* out_fp,
                             const char* name,
                             size_t pos,
                             char ref,
                             char alt);

/*
 * Ferme un fichier CSV précédemment ouvert.
 *
 * Paramètres :
 *  - out_fp : FILE* à fermer (peut être NULL).
 *
 * Remarque :
 *  - Appeler cette fonction avec out_fp == NULL n’a aucun effet.
 */
void csv_close(FILE* out_fp);

#endif
