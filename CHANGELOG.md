# CHANGELOG

## 2025-11-25

### Header Consolidation
- **Merged `Globals.h`** from both platforms into `src/Globals.h`:
  - Platform-specific MPI headers and system includes are conditionally included via `#ifdef PLATFORM_SUNWAY`
  - `slave_num` macro is defined for Sunway platform only
  - Original platform-specific files backed up as `.bak`
  
- **Merged `th_ass.h`** from both platforms into `src/th_ass.h`:
  - Thread assignment arrays differ between platforms (x86 has larger arrays)
  - Platform-specific arrays are conditionally compiled using `#ifdef PLATFORM_X86` / `#ifdef PLATFORM_SUNWAY`
  - Original platform-specific files backed up as `.bak`

- **Merged `cmdinfo.h`** from both platforms into `src/cmdinfo.h`:
  - x86-specific: `notKeepOrder_`, `use_igzip_`, and CARE-related parameters
  - Sunway-specific: `CmdInfoCopy()` method and `splitWrite_` flag
  - Common parameters shared between platforms
  - Original platform-specific files backed up as `.bak`

### Build System Improvements
- Fixed `PLATFORM_SUNWAY` macro definition in `sunway/CMakeLists.txt`:
  - Added `PLATFORM_SUNWAY` to `host_objects` and `slave_objects` targets
  - Added include directories to `rabbitqc-x` target to resolve `src/Globals.h` and `src/cmdinfo.h`
  
- Updated all include paths across the codebase:
  - Changed `#include "Globals.h"` to `#include "src/Globals.h"` in platform-specific files
  - Changed `#include "cmdinfo.h"` to `#include "src/cmdinfo.h"` in platform-specific files
  - Unified `src/main.cpp` to include common headers directly from `src/`

### Code Organization
- Improved include order consistency:
  - System headers (`<...>`) first
  - Unified headers (`src/Globals.h`, `src/cmdinfo.h`) second
  - Platform-specific headers last
- Removed duplicate MPI header includes (now handled by `src/Globals.h`)

## 2025-11-25 (Earlier)

- Initial integration of **RabbitQCPlus** (x86) and **SWQC** (Sunway) into a unified project `RabbitQC-X`.
- Added CMake-based build system with `PLATFORM` switch (`x86` / `sunway`).
- Unified the main entry point into `src/main.cpp` with `#ifdef PLATFORM_X86` / `#ifdef PLATFORM_SUNWAY`.
- Moved shared headers and libraries into `common/` and `common/lib/`, keeping platform-specific differences behind `#ifdef` where needed.
- Fixed GCC compatibility on x86 (supporting GCC 4.8.5+), cleaned up compiler warnings, and improved include paths.
- Added clear platform-specific include tags in `src/main.cpp` so headers are resolved from `x86/` or `sunway/` explicitly.


