# Star Tracker - Embedded-Friendly Star Tracker in C

A lightweight, embedded-friendly star tracker implementation in C. Performs centroid extraction, Hipparcos star catalog management, and attitude determination using TRIAD and QUEST algorithms.

## Features

- **Centroid Processing** — Centroid list management with flux-based sorting, position and error tracking
- **Hipparcos Star Catalog** — Load/filter stars by magnitude, RA/Dec to unit vector conversion, CSV import support
- **Camera Model** — Pinhole camera model with pixel-to-unit-vector projection and FOV computation
- **Attitude Determination**
  - **TRIAD** — Fast 2-vector attitude solution, ideal for real-time embedded use
  - **QUEST** — Optimal N-vector quaternion estimator using Newton-Raphson on Davenport's eigenvalue problem
- **Coordinate Conversions** — DCM ↔ Quaternion ↔ Euler angles, angular distance computation

## Architecture

```
star_tracker/
├── include/
│   ├── centroid.h        # Centroid data structures and API
│   ├── catalog.h         # Hipparcos catalog management
│   ├── star_tracker.h    # Camera model and projections
│   └── attitude.h        # TRIAD, QUEST, DCM, Quaternion
├── src/
│   ├── centroid.c         # Centroid list: add, sort, print
│   ├── catalog.c          # Catalog: CSV load, magnitude filter, unit vectors
│   ├── star_tracker.c     # Pixel→unit vector, angular distance
│   └── attitude.c         # TRIAD & QUEST algorithms, conversions
├── test/
│   └── main.c             # Demo: simulated centroids + attitude recovery
├── data/                   # Catalog CSV files (user-provided)
├── Makefile
└── README.md
```

## Embedded Design Principles

- **No dynamic memory allocation** — All structures use fixed-size arrays (`MAX_CENTROIDS_PER_FRAME = 256`, `CATALOG_MAX_STARS = 9110`)
- **No external dependencies** — Only standard C library (`math.h`, `stdio.h`, `string.h`, `stdint.h`)
- **C99 standard** — Compatible with most embedded compilers (GCC, ARM-GCC, MSVC)
- **Float arithmetic** — Uses `float` (not `double`) for FPU-constrained targets
- **Insertion sort** — O(n²) but no recursion, no malloc, stack-safe

## Building

### Prerequisites

- GCC (or any C99 compiler)
- Make

### Compile & Run

```bash
make clean && make
./build/star_tracker
```

### Cross-compilation (example: ARM Cortex-M)

```bash
make CC=arm-none-eabi-gcc CFLAGS="-mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard -std=c99 -O2 -Iinclude"
```

## Usage Example

### Centroid List

```c
CentroidList frame;
centroid_list_init(&frame, 1);  // frame_id = 1

// Add detected centroids: x, y, flux, x_error, y_error
centroid_list_add(&frame, 512.3f, 510.7f, 45000.0f, 0.12f, 0.11f);
centroid_list_add(&frame, 234.1f, 678.9f, 32000.0f, 0.15f, 0.14f);

// Sort by brightness (brightest first)
centroid_list_sort_by_flux(&frame);
centroid_list_print(&frame);
```

### Camera & Projection

```c
CameraParams cam;
camera_init(&cam, 50.0f, 0.0055f, 1024, 1024);  // 50mm focal, 5.5µm pixel

float ux, uy, uz;
pixel_to_unit_vector(&cam, 512.0f, 512.0f, &ux, &uy, &uz);  // center pixel → boresight
```

### Attitude Determination (TRIAD)

```c
MatchList matches;
match_list_init(&matches);

// Add matched star pairs: body vector + reference vector + weight
match_list_add(&matches, bx1, by1, bz1, rx1, ry1, rz1, 1.0f);
match_list_add(&matches, bx2, by2, bz2, rx2, ry2, rz2, 1.0f);

DCM dcm;
attitude_triad(&matches, &dcm);

float roll, pitch, yaw;
dcm_to_euler(&dcm, &roll, &pitch, &yaw);
```

### Attitude Determination (QUEST)

```c
// Uses N>=2 matched pairs for optimal quaternion estimate
Quaternion q;
attitude_quest(&matches, &q);
quaternion_print(&q);

// Convert to DCM or Euler
DCM dcm;
quaternion_to_dcm(&q, &dcm);
dcm_to_euler(&dcm, &roll, &pitch, &yaw);
```

## Demo Output

```
TRUE attitude (Euler angles):
  Roll:  +15.00 deg
  Pitch: -25.00 deg
  Yaw:   +40.00 deg

TRIAD Euler: Roll=+15.0000  Pitch=-25.0000  Yaw=+40.0000
TRIAD error vs truth: 0.000000 deg

QUEST Euler: Roll=+15.0000  Pitch=-25.0000  Yaw=+40.0000
QUEST error vs truth: 0.000000 deg
```

## Algorithms

### TRIAD

Uses two body-reference vector pairs to construct orthonormal triads in both frames, then computes the Direction Cosine Matrix (DCM) as:

```
C = M_body × M_ref^T
```

Fast and deterministic. Best when exactly 2 high-confidence matches are available.

### QUEST (QUaternion ESTimator)

Solves Wahba's problem optimally by finding the quaternion that maximizes:

```
J(C) = Σ wᵢ · bᵢᵀ · C · rᵢ
```

Uses Newton-Raphson iteration on the characteristic equation of Davenport's K-matrix to find the maximum eigenvalue, then extracts the quaternion via the adjugate method. No matrix decomposition required — fully embedded-friendly.

## Hipparcos Catalog

The implementation supports loading the Hipparcos catalog from CSV format:

```csv
hip_id,ra,dec,mag
32349,101.287,-16.716,-1.46
30438,95.988,-52.696,-0.74
...
```

A sample of the 10 brightest stars is hardcoded for testing. Full catalog CSV can be placed in the `data/` directory and loaded at runtime via `catalog_load_csv()`.

## License

MIT
