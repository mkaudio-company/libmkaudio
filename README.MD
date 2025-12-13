# MKAudio

A C++ library for high-performance audio signal processing, focusing on analog modeling techniques optimized for real-time applications.

## Overview

MKAudio explores two complementary approaches to analog audio modeling:

1. **Numeric Modeling** - Lightweight saturation approximation using mathematical functions
2. **Circuit Simulation** - Real-time transient analysis of actual circuit schematics

Both approaches are designed for low memory footprint and minimal computational overhead, making them suitable for real-time audio processing.

## Modules

### Saturation (Numeric Modeling)

The saturation module implements a novel approach to modeling analog saturation characteristics without circuit simulation or impulse responses. The core idea is to use logarithmic functions (`log2`) with carefully tuned parameters to approximate the nonlinear transfer functions found in analog equipment.

The model supports:
- **Asymmetric saturation** - Independent control over positive and negative signal portions via separate alpha (drive/knee) parameters
- **Compression control** - Beta parameters for adjusting gain/compression characteristics on each polarity
- **DC bias** - Delta parameter shifts the saturation curve's center point
- **Polarity flip** - Optional signal inversion

This approach can be chained: saturation stage → filter → saturation stage → filter, building complex analog-style processing chains while maintaining efficiency.

**Key files:** `saturation.h`, `saturation.cpp`

### Circuit Simulation (Real-Time Transient Analysis)

The simulation module provides authentic analog filtering by performing real-time transient analysis of circuit schematics. Audio signals are treated as source voltages, and the circuit response is computed sample-by-sample using Modified Nodal Analysis (MNA).

Features:
- **Component library** - Resistors, capacitors, and inductors with companion model discretization
- **Templated circuit engine** - Compile-time fixed sizes (`NumDevices`, `NumNodes`) for maximum performance
- **Gaussian elimination solver** - Optimized linear system solver with partial pivoting
- **Preprocessing step** - Static admittance matrix (Y) is built once before audio processing
- **Per-sample processing** - Dynamic current sources updated each sample for reactive components

This enables modeling of actual analog filter topologies (RC filters, LC networks, etc.) with physically accurate frequency response and transient behavior.

**Key files:** `simulation.h`, `simulation.cpp`

## Building

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

Requires CMake 3.22+ and a C++17 compatible compiler.

## Usage

### Saturation

```cpp
#include "saturation.h"

// Create saturation with asymmetric parameters
MKAudio::Saturation sat(
    10.0,   // alphaPlus (positive drive)
    8.0,    // alphaMinus (negative drive)
    1.0,    // betaPlus (positive compression)
    1.0,    // betaMinus (negative compression)
    0.0,    // deltaBias
    false   // flipPolarity
);

double output = sat.process(inputSample);
```

### Circuit Simulation

```cpp
#include "simulation.h"

// Create circuit: 4 devices max, 3 nodes max
MKAudio::Circuit<4, 3> circuit(44100.0);

// Build an RC lowpass filter
circuit.addComponent(std::make_shared<MKAudio::Resistor>(1, 2, 1000.0));   // 1k ohm
circuit.addComponent(std::make_shared<MKAudio::Capacitor>(2, 0, 1e-7));    // 100nF

// Preprocess (call once before audio starts)
circuit.preprocess(0.1);  // source impedance

// Process audio
double output = circuit.process(inputVoltage, 2);  // probe node 2
```

## License

This library is dual-licensed:

**GPL-3.0** for open source projects

**Commercial license** available for closed source usage

For commercial licensing inquiries, contact: minjaekim@mkaudio.company
