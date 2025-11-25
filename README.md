# RabbitQC-X

A unified, cross-platform quality control tool for sequencing data that supports both x86 and Sunway architectures.

## Overview

RabbitQC-X combines two high-performance implementations:
- **x86 Platform**: RabbitQCPlus - optimized for x86 architecture with AVX2/AVX512 support
- **Sunway Platform**: SWQC - optimized for Sunway supercomputer architecture with heterogeneous computing

If you encounter issues with RabbitQC-X in your environment, you can always fall back to the original standalone projects:
- **Original RabbitQCPlus (x86)**: [RabbitBio/RabbitQCPlus](https://github.com/RabbitBio/RabbitQCPlus)
- **Original SWQC (Sunway)**: [RabbitBio/SWQC](https://github.com/RabbitBio/SWQC)

These upstream projects are kept as close as possible to their original behavior; RabbitQC-X focuses on unifying the build and code layout.

## Features

### x86 Platform (RabbitQCPlus)
- Single-threaded performance improved by 2x
- 4x+ speedup for gzip-compressed files compared to SOAPnuke
- Integrated and optimized CARE error correction engine (1.3x speedup)
- Automatic instruction set detection (AVX2/AVX512)
- Support for both NGS and TGS data

### Sunway Platform (SWQC)
- First distributed quality control software for Sunway platform
- 1.5-3x faster than fastest x86 software on single node
- 65%-80% scalability on 16 nodes
- Optimized multi-threaded (de)compression framework

## Installation

### System Requirements

#### For x86 Platform:
- 64-bit Linux system
- GCC 7.5.0 or newer (GCC 9.0+ recommended for best performance)
- zlib library
- OpenMP support

#### For Sunway Platform:
- Sunway next-generation platform
- sw9gcc 7.1.0 or newer
- zlib library
- MPI support

### Building from Source

#### 1. Clone the repository

```bash
git clone <repository-url>
cd RabbitQC-X
```

#### 2. Build for x86 platform

```bash
mkdir build
cd build
cmake -DPLATFORM=x86 ..
make -j4
```

#### 3. Build for Sunway platform

```bash
mkdir build-sunway
cd build-sunway
cmake -DPLATFORM=sunway ..
make -j8
```

## Project Structure

```
RabbitQC-X/
├── src/                 # Unified entry point
│   ├── main.cpp        # Unified main program (platform-specific via #ifdef)
│   └── CLI11.hpp       # Command-line parsing library
├── common/              # Shared code between platforms
│   ├── lib/            # Unified library files (merged from both platforms)
│   ├── include/        # Shared header files
│   ├── common.hpp
│   ├── exceptions.hpp
│   └── system.h
├── x86/                 # x86 platform specific code
│   ├── src/            # x86 source files
│   ├── include/        # x86 header files
│   ├── dependencies/   # Third-party dependencies (Thrust)
│   └── CMakeLists.txt
├── sunway/             # Sunway platform specific code
│   ├── src/            # Host source files
│   ├── slave/          # Slave core source files
│   └── CMakeLists.txt
├── CMakeLists.txt      # Top-level CMake configuration
└── README.md           # This file
```

## Usage

### For Next Generation Sequencing (NGS) Data

#### Single-end (SE) reads

```bash
# Plain FASTQ
./rabbitqc-x -w 8 -i in1.fastq -o p1.fastq

# Gzip compressed
./rabbitqc-x -w 8 -i in1.fastq.gz -o p1.fastq.gz
```

#### Paired-end (PE) reads

```bash
# Plain FASTQ
./rabbitqc-x -w 8 -i in1.fastq -I in2.fastq -o p1.fastq -O p2.fastq

# Gzip compressed
./rabbitqc-x -w 16 -i in1.fastq.gz -I in2.fastq.gz -o p1.fastq.gz -O p2.fastq.gz
```

#### With error correction (x86 only)

```bash
# SE mode
./rabbitqc-x -w 32 -i in1.fastq -o p1.fastq --correctWithCare --coverage 30 --pairmode SE

# PE mode
./rabbitqc-x -w 32 -i in1.fastq -I in2.fastq -o p1.fastq -O p2.fastq --correctWithCare --coverage 30 --pairmode PE
```

### For Third Generation Sequencing (TGS) Data (x86 only)

```bash
# Plain FASTQ
./rabbitqc-x -w 4 -i in.fastq --TGS

# Gzip compressed
./rabbitqc-x -w 6 -i in.fastq.gz --TGS
```

## Command-line Options

For complete help information:

```bash
./rabbitqc-x -h
```

Common options:
- `-w, --thread`: Number of threads to use
- `-i, --in1`: Input file 1 (or single-end file)
- `-I, --in2`: Input file 2 (for paired-end)
- `-o, --out1`: Output file 1 (or single-end file)
- `-O, --out2`: Output file 2 (for paired-end)
- `--TGS`: Enable TGS mode (third-generation sequencing)
- `--correctWithCare`: Enable CARE error correction (x86 only)

## Performance

The x86 version has been benchmarked against RabbitQC (v0.0.1), fastp (v0.23.2), SOAPnuke (v2.1.7), Trimmomatic (v0.40), and CARE (v2.0.0) using 370 million Illumina sequencing reads.

Key results:
- Significant speedup on both plain and gzip-compressed FASTQ files
- Improved performance for over-representation analysis (5x faster)
- Efficient error correction with CARE engine integration

## Visual Output

The tool generates comprehensive HTML reports with before/after filtering statistics. See example reports for detailed visualization of quality control metrics.

## License

See LICENSE files in respective platform directories.

## Citation

### For x86 Platform (RabbitQCPlus):

Lifeng Yan, Zekun Yin, Hao Zhang, Zhan Zhao, Mingkai Wang, André Müller, Felix Kallenborn et al. "RabbitQCPlus 2.0: More efficient and versatile quality control for sequencing data." Methods 216 (2023): 39-50.

Lifeng Yan, Zekun Yin, Hao Zhang, Zhan Zhao, Mingkai Wang, André Müller, Robin Kobus, Yanjie Wei, Beifang Niu, Bertil Schmidt, Weiguo Liu. "RabbitQCPlus: More Efficient Quality Control for Sequencing Data," *2022 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*, Las Vegas, NV, USA, 2022, pp. 619-626, doi: 10.1109/BIBM55620.2022.9995332.

Zekun Yin, Hao Zhang, Meiyang Liu, Wen Zhang, Honglei Song, Haidong Lan, Yanjie Wei, Beifang Niu, Bertil Schmidt, Weiguo Liu, RabbitQC: High-speed scalable quality control for sequencing data, Bioinformatics, , btaa719, https://doi.org/10.1093/bioinformatics/btaa719

### For Sunway Platform (SWQC):

Yan L, Yin Z, Zhang T, et al. SWQC: Efficient sequencing data quality control on the next-generation sunway platform[J]. Future Generation Computer Systems, 2025, 164: 107577.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Contact

For questions and support, please open an issue in the repository.

