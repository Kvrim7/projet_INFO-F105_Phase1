#include "utils_fasta.h"
#include <cstring> // strlen, memcpy

static void strip_newlines(char* s) {
    if (!s) return;
    size_t n = std::strlen(s);
    while (n > 0 && (s[n - 1] == '\n' || s[n - 1] == '\r')) {
        s[n - 1] = '\0';
        --n;
    }
}

FastaStatus fasta_open_and_read_header(const char* path, FILE** out_fp,
                                       char* name, size_t name_cap) {
    if (!path || !out_fp || !name || name_cap < 2) return FASTA_ERR_FORMAT;

    *out_fp = std::fopen(path, "r");
    if (!*out_fp) return FASTA_ERR_OPEN;

    char line[2048];
    if (!std::fgets(line, (int)sizeof(line), *out_fp)) {
        fasta_close(*out_fp);
        *out_fp = nullptr;
        return FASTA_ERR_EMPTY;
    }

    strip_newlines(line);
    if (line[0] != '>' || line[1] == '\0') {
        fasta_close(*out_fp);
        *out_fp = nullptr;
        return FASTA_ERR_FORMAT;
    }

    const char* src = line + 1;
    size_t nlen = std::strlen(src);
    if (nlen + 1 > name_cap || nlen > 100) { 
        fasta_close(*out_fp);
        *out_fp = nullptr;
        return FASTA_ERR_NAME_TOO_LONG;
    }

    std::memcpy(name, src, nlen + 1);
    return FASTA_OK;
}

FastaStatus fasta_read_bases(FILE* fp, char* buf, size_t cap, size_t* out_n) {
    if (!fp || !buf || cap == 0 || !out_n) return FASTA_ERR_FORMAT;

    size_t n = 0;
    while (n < cap) {
        int c = std::fgetc(fp);
        if (c == EOF) break;

        if (c == '\n' || c == '\r') continue; 
        if (c == '>') return FASTA_ERR_FORMAT; 

        buf[n++] = (char)c;
    }

    if (std::ferror(fp)) return FASTA_ERR_IO;

    *out_n = n;            // 0 => EOF (pas une erreur)
    return FASTA_OK;
}

void fasta_close(FILE* fp) {
    if (fp) std::fclose(fp);
}
