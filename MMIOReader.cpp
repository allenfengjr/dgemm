/// ---------------------------------------------------------------------------
///   HPC Course Fall 2022
///   Indiana University Bloomington
///
///   Luke D'Alessandro
///
///   NOTES: Reader Taken from: https://math.nist.gov/MatrixMarket/mmio-c.html
///
/// ---------------------------------------------------------------------------
#include "MMIOReader.hpp"
#include "mmio.h"
#include <fmt/format.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

MMIOReader::MMIOReader(std::string path)
{
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        fmt::print(stderr, "{} open failed, {}: {}\n", path, errno, strerror(errno));
        std::exit(EXIT_FAILURE);
    }

    FILE* f = fdopen(fd, "r");
    if (f == nullptr) {
        fmt::print(stderr, "{} fdopen failed, {}: {}\n", path, errno, strerror(errno));
        close(fd);
        std::exit(EXIT_FAILURE);
    }

    MM_typecode type;
    switch (mm_read_banner(f, &type)) {
        case MM_PREMATURE_EOF:    // if all items are not present on first line of file.
        case MM_NO_HEADER:        // if the file does not begin with "%%MatrixMarket".
        case MM_UNSUPPORTED_TYPE: // if not recongizable description.
            fmt::print(stderr, "Could not process Matrix Market banner for {}.\n", path);
            close(fd);
            fclose(f);
            std::exit(1);
    }

    if (!mm_is_coordinate(type)) {
        fmt::print(stderr, "{} does not appear to be a coordinate matrix", path);
        close(fd);
        fclose(f);
        std::exit(1);
    }

    if (mm_is_complex(type) && mm_is_matrix(type) && mm_is_sparse(type)) {
        fmt::print(stderr, "Sorry, this application does not support ");
        fmt::print(stderr, "Market Market type: [{}]\n", mm_typecode_to_str(type));
        close(fd);
        fclose(f);
        exit(1);
    }

    // Check if the matrix is symmetric.
    is_symmetric = mm_is_symmetric(type);

    // Check if the matrix has edge weights or just geometry
    has_weights = not mm_is_pattern(type);


    switch (mm_read_mtx_crd_size(f, &n_rows, &n_columns, &n_non_zeros)) {
        case MM_PREMATURE_EOF:    // if an end-of-file is encountered before processing these three values.
            close(fd);
            fclose(f);
            std::exit(1);
    }

    long n = ftell(f);
    fseek(f, 0L, SEEK_END);
    long bytes = ftell(f);

    _data = static_cast<char*>(mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0));
    if (_data == MAP_FAILED) {
        fmt::print(stderr, "{} mmap failed, {}: {}\n", path, errno, strerror(errno));
        close(fd);
        fclose(f);
        std::exit(1);
    }

    _edges = _data + n;
    _end = _data + bytes;

    close(fd);
    fclose(f);
}

MMIOReader::~MMIOReader()
{
    if (_data && munmap((char*)_data, _end - _data)) {
        fmt::print(stderr, "munmap failed, {}: {}\n", errno, strerror(errno));
    }
}