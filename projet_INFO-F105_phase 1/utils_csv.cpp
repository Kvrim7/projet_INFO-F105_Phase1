#include "utils_csv.h"

CsvStatus csv_open(const char* path, FILE** out_fp) {
    if (!path || !out_fp) return CSV_ERR_OPEN;
    *out_fp = std::fopen(path, "w");
    return (*out_fp) ? CSV_OK : CSV_ERR_OPEN;
}

CsvStatus csv_write_mutation(FILE* out_fp, const char* name, size_t pos, char ref, char alt) {
    if (!out_fp || !name) return CSV_ERR_IO;
    if (std::fprintf(out_fp, "%s,%zu,%c,%c\n", name, pos, ref, alt) < 0) return CSV_ERR_IO;
    return std::ferror(out_fp) ? CSV_ERR_IO : CSV_OK;
}

void csv_close(FILE* out_fp) {
    if (out_fp) std::fclose(out_fp);
}
