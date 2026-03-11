#ifndef UTILS_FASTA_H
#define UTILS_FASTA_H

#include <cstddef> // size_t
#include <cstdio>  // FILE

/*
 * Codes de retour des fonctions utilitaires FASTA.
 *
 * Les fonctions utilitaires :
 *  - n’affichent aucun message,
 *  - ne terminent jamais le programme.
 *
 * Il est de la responsabilité de votre programme :
 *  - de vérifier la valeur de retour,
 *  - d’afficher un message d’erreur explicite si nécessaire,
 *  - de terminer proprement le programme avec un code de retour non nul.
 */
enum FastaStatus {
    FASTA_OK = 0,                // Opération réussie
    FASTA_ERR_OPEN = 1,          // Impossible d’ouvrir le fichier
    FASTA_ERR_EMPTY = 2,         // Fichier vide
    FASTA_ERR_FORMAT = 3,        // Format FASTA invalide (en-tête manquant, plusieurs séquences, etc.)
    FASTA_ERR_NAME_TOO_LONG = 4, // Nom de séquence trop long
    FASTA_ERR_IO = 5             // Erreur d’entrée/sortie lors de la lecture
};

/*
 * Ouvre un fichier FASTA et lit l’en-tête de la séquence.
 *
 * Paramètres :
 *  - path : chemin vers le fichier FASTA à ouvrir.
 *  - out_fp : adresse d’une variable FILE*.
 *             En cas de succès, *out_fp contient un FILE* valide ouvert en lecture.
 *             Le curseur du fichier est positionné au début de la séquence
 *             (juste après la ligne d’en-tête).
 *  - name : tableau de caractères dans lequel le nom de la séquence sera stocké
 *           (le caractère '>' n’est PAS inclus).
 *  - name_cap : capacité du tableau name, incluant l’espace pour le caractère '\0'.
 *
 * Valeur de retour :
 *  - FASTA_OK en cas de succès.
 *  - Un code d’erreur sinon (fichier inexistant, format invalide, nom trop long, etc.).
 *
 * Remarques :
 *  - Pour ce projet, chaque fichier FASTA doit contenir exactement UNE seule séquence.
 *  - Si le fichier ne commence pas par une ligne valide de type '>nom', cela constitue une erreur.
 */
FastaStatus fasta_open_and_read_header(const char* path,
                                       FILE** out_fp,
                                       char* name,
                                       size_t name_cap);

/*
 * Lit des bases de la séquence FASTA dans un buffer.
 *
 * Paramètres :
 *  - fp : FILE* obtenu précédemment via fasta_open_and_read_header.
 *  - buf : buffer dans lequel les bases lues seront stockées.
 *  - cap : capacité du buffer (nombre maximal de bases à lire).
 *  - out_n : pointeur vers une variable de type size_t dans laquelle sera stocké
 *            le nombre de bases effectivement lues.
 *
 * Comportement :
 *  - La fonction lit au maximum 'cap' caractères de la séquence.
 *  - Les caractères '\n' et '\r' sont ignorés silencieusement et ne constituent PAS une erreur.
 *  - Si la séquence contient moins de 'cap' bases restantes, toutes les bases restantes sont lues.
 *  - Lorsque la fin de la séquence est atteinte (fin de fichier), *out_n vaut 0.
 *
 * Valeur de retour :
 *  - FASTA_OK si la lecture s’est déroulée correctement
 *    (y compris lorsque la fin de la séquence est atteinte).
 *  - FASTA_ERR_FORMAT si un second en-tête FASTA ('>') est rencontré
 *    (plusieurs séquences ne sont pas autorisées dans ce projet).
 *  - FASTA_ERR_IO en cas d’erreur d’entrée/sortie.
 */
FastaStatus fasta_read_bases(FILE* fp,
                             char* buf,
                             size_t cap,
                             size_t* out_n);

/*
 * Ferme un fichier FASTA précédemment ouvert.
 *
 * Paramètres :
 *  - fp : FILE* à fermer (peut être NULL).
 *
 * Remarque :
 *  - Appeler cette fonction avec fp == NULL n’a aucun effet.
 */
void fasta_close(FILE* fp);

#endif

