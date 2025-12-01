#ifndef GLOBALMUTEX_H
#define GLOBALMUTEX_H

// USE_LIBDEFLATE, CONSUMER_USE_LIBDEFLATE, USE_CC_GZ are always enabled (code is written directly)
// WRITER_USE_LIBDEFLATE is always disabled

#include <mutex>
#define Q_lim_se 16
#define Q_lim_pe 8

#define BLOCK_SIZE (1 << 20)
#define SWAP1_SIZE (1 << 18)
#define SWAP2_SIZE (1 << 14)
#define USE_CC_GZ
extern std::mutex globalMutex;
struct Para {
    char* in_buffer;
    char* out_buffer;
    size_t *out_size;
    size_t in_size;
    int level;
};

//struct QChunkItem {
//    char* buffer[64];
//    int buffer_len[64];
//    long long file_offset[64];
//};

#endif
