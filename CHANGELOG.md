# CHANGELOG

## 2025-11-25

- Initial integration of **RabbitQCPlus** (x86) and **SWQC** (Sunway) into a unified project `RabbitQC-X`.
- Added CMake-based build system with `PLATFORM` switch (`x86` / `sunway`).
- Unified the main entry point into `src/main.cpp` with `#ifdef PLATFORM_X86` / `#ifdef PLATFORM_SUNWAY`.
- Moved shared headers and libraries into `common/` and `common/lib/`, keeping platform-specific differences behind `#ifdef` where needed.
- Fixed GCC compatibility on x86 (supporting GCC 4.8.5+), cleaned up compiler warnings, and improved include paths.
- Added clear platform-specific include tags in `src/main.cpp` so headers are resolved from `x86/` or `sunway/` explicitly.


