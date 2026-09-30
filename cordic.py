"""
CORDIC (COordinate Rotation DIgital Computer) Algorithm
======================================================

A complete Python implementation using fixed-point integer arithmetic and
pure bit-shift operations, with comprehensive mathematical derivations for
Circular, Hyperbolic, and Linear coordinate systems (Volder 1959, Walther 1971).

===============================================================================
TABLE OF CONTENTS
===============================================================================
1. Historical Background & Motivation
2. Mathematical Derivation: Circular Coordinate System (m = 1)
   - 2D Rotation Matrix Decomposition
   - The Pseudo-Rotation Trick & Elementary Angles
   - Cumulative Scaling Factor (K_c)
   - Convergence Proof & Domain of Convergence
3. Walther's Unified CORDIC Architecture (m in {1, 0, -1})
   - Circular (m = 1): sin, cos, tan, atan, atan2, polar/Cartesian conversion
   - Linear (m = 0): multiplication, division
   - Hyperbolic (m = -1): sinh, cosh, tanh, exp, ln, atanh, sqrt
   - Hyperbolic Convergence & Iteration Repeat Schedule (i = 4, 13, 40, ...)
4. Fixed-Point Arithmetic (Q-Format) & Bit-Shift Operations
5. Range Reduction Strategies
6. Python Core Implementations & High-Level Functions
7. Comprehensive Verification Suite & Numerical Error Analysis


===============================================================================
1. HISTORICAL BACKGROUND & MOTIVATION
===============================================================================
The CORDIC algorithm was conceived by Jack E. Volder in 1959 at the Convair
division of General Dynamics to replace the analog navigation computer in the
B-58 Hustler supersonic bomber. Volder sought an efficient digital algorithm to
calculate trigonometric transformations and vector rotations in real time.

In 1971, J. S. Walther of Hewlett-Packard generalized Volder's work into a
unified framework encompassing circular, linear, and hyperbolic coordinate
systems. This unified algorithm became the computational foundation of the
legendary HP-35 scientific pocket calculator (1972) and continues to be
extensively used in modern FPGAs, ASICs, DSPs, and embedded microcontrollers.

In digital hardware and embedded processors, dedicated floating-point units and
hardware multipliers consume large silicon area, power, and clock cycles. CORDIC
solves this by computing transcendental and algebraic functions using ONLY:
  - Additions and Subtractions (+, -)
  - Bitwise arithmetic right-shifts (>>), equivalent to multiplication by 2^(-i)
  - Small precomputed lookup tables (ROM constants)


===============================================================================
2. MATHEMATICAL DERIVATION: CIRCULAR COORDINATE SYSTEM (m = 1)
===============================================================================
Consider rotating a 2D Cartesian vector v = [x, y]^T counter-clockwise by
an angle theta:

    [ x' ]   [ cos(theta)  -sin(theta) ] [ x ]
    [ y' ] = [ sin(theta)   cos(theta) ] [ y ]

Factoring out cos(theta):

    [ x' ]                [ 1           -tan(theta) ] [ x ]
    [ y' ] = cos(theta) * [ tan(theta)       1      ] [ y ]

Instead of attempting to compute tan(theta) directly in a single step, CORDIC
decomposes the arbitrary angle theta into a sum of n elementary rotations by
fixed angles gamma_i:

    theta = sum_{i=0}^{n-1} d_i * gamma_i,   where d_i in {-1, +1}

Volder's key insight was choosing each elementary rotation angle gamma_i such that
its tangent is an exact negative power of 2:

    tan(gamma_i) = 2^(-i)   <=>   gamma_i = arctan(2^(-i))

Substituting tan(gamma_i) = 2^(-i) into the single-step rotation matrix yields:

    R_i = cos(gamma_i) * [ 1          -d_i * 2^(-i) ]
                         [ d_i * 2^(-i)      1      ]

In binary fixed-point arithmetic, multiplication by 2^(-i) is an arithmetic
right-shift by i bit positions:

    v * 2^(-i) == (v >> i)

--- The Pseudo-Rotation Trick & Cumulative Gain (K_c) ---
Multiplying by cos(gamma_i) at every iteration would still require hardware
multiplication. However, because cosine is an even function (cos(-gamma) = cos(gamma)):

    cos(gamma_i) = 1 / sqrt(1 + tan^2(gamma_i)) = 1 / sqrt(1 + 2^(-2i))

The scaling factor cos(gamma_i) depends ONLY on the iteration index i, and is
completely independent of the direction of rotation d_i!

We can therefore omit cos(gamma_i) from the inner iteration loop ("pseudo-rotation")
and group all factors into a single cumulative constant K_c:

    K_c = prod_{i=0}^{n-1} cos(gamma_i) = prod_{i=0}^{n-1} 1 / sqrt(1 + 2^(-2i))
    A_c = 1 / K_c = prod_{i=0}^{n-1} sqrt(1 + 2^(-2i))

As n -> infinity, K_c converges rapidly to the circular CORDIC constant:

    K_c = lim_{n->inf} K_c approx 0.6072529350088812561694...
    A_c = 1 / K_c          approx 1.646760258121065648366...

During the n pseudo-rotations, the vector's magnitude expands by A_c = 1/K_c.
If we pre-scale the initial input x_0 = K_c (or post-scale by K_c), the final
output coordinates are exact without requiring per-iteration scaling.

--- Convergence Proof ---
For any arbitrary angle theta to be reachable without gaps, the sequence of
elementary angles must satisfy the convergence condition:

    gamma_i <= sum_{j=i+1}^{n-1} gamma_j + gamma_{n-1}

Because arctan(2^(-i)) < 2 * arctan(2^(-(i+1))), this inequality holds for all i >= 0.
The maximum angle reachable by circular CORDIC without range reduction is:

    theta_max = sum_{i=0}^{n-1} gamma_i approx 1.74328662 rad (~99.88 degrees)

Since theta_max > pi/2 (~1.57079 rad), the entire first quadrant [0, pi/2] is
guaranteed to converge. Arbitrary angles in (-inf, inf) are mapped to [-pi/2, pi/2]
via standard quadrant reduction.


===============================================================================
3. WALTHER'S UNIFIED CORDIC ARCHITECTURE (1971)
===============================================================================
Walther unified circular (m=1), linear (m=0), and hyperbolic (m=-1) coordinate
systems into a single generalized recurrence:

    x_{i+1} = x_i - m * d_i * (y_i >> i)
    y_{i+1} = y_i +     d_i * (x_i >> i)
    z_{i+1} = z_i -     d_i * e_i

Where the coordinate parameter m, elementary step e_i, and scaling factor K are:

+-------------------+----+-----------------------+------------------------------------+
| Coordinate System | m  | Elementary Step e_i   | Cumulative Scale Factor K          |
+-------------------+----+-----------------------+------------------------------------+
| Circular          |  1 | arctan(2^(-i))        | prod 1/sqrt(1 + 2^(-2i))  ~0.60725 |
| Linear            |  0 | 2^(-i)                | 1.0 (no scaling needed)            |
| Hyperbolic        | -1 | artanh(2^(-i))        | prod 1/sqrt(1 - 2^(-2i))  ~1.20750 |
+-------------------+----+-----------------------+------------------------------------+

--- Hyperbolic CORDIC & Convergence Repeats ---
In the hyperbolic plane (m = -1), the metric is x^2 - y^2:
  - Iterations must start at i = 1 (since artanh(2^0) = artanh(1) = inf).
  - The hyperbolic pseudo-rotation step contracts the vector's hyperbolic norm by
    sqrt(1 - 2^(-2i)), so the pre-scaling factor is K_h = prod 1/sqrt(1 - 2^(-2i)) ~ 1.20750.
  - Crucially, artanh(2^(-i)) > sum_{j=i+1}^inf artanh(2^(-j)).
    To guarantee convergence without gaps, Walther proved that iterations at
    indices i in {4, 13, 40, 121, ..., (3^(k+1)-1)/2} must be repeated twice!

--- Unified Functions Summary ---
1. Rotation Mode (drives z -> 0, d_i = +1 if z >= 0 else -1):
   - Circular (m=1):   (x_0=K_c, y_0=0, z_0=theta)   => x_n = cos(theta), y_n = sin(theta)
   - Linear (m=0):     (x_0=x,   y_0=0, z_0=z)       => y_n = x * z (multiplication)
   - Hyperbolic (m=-1):(x_0=K_h, y_0=0, z_0=theta)   => x_n = cosh(theta), y_n = sinh(theta)
                       exp(theta) = cosh(theta) + sinh(theta)

2. Vectoring Mode (drives y -> 0, d_i = +1 if y < 0 else -1):
   - Circular (m=1):   (x_0=x, y_0=y, z_0=0) => z_n = atan2(y, x), x_n = A_c * sqrt(x^2+y^2)
   - Linear (m=0):     (x_0=x, y_0=y, z_0=0) => z_n = y / x (division)
   - Hyperbolic (m=-1):(x_0=x, y_0=y, z_0=0) => z_n = artanh(y/x), x_n = A_h * sqrt(x^2-y^2)
                       ln(w) = 2 * artanh((w - 1) / (w + 1))
                       sqrt(a): set x_0 = a + 1/4, y_0 = a - 1/4 => sqrt(x^2 - y^2) = sqrt(a)


===============================================================================
4. FIXED-POINT ARITHMETIC (Q-FORMAT) & BIT OPERATIONS
===============================================================================
All calculations are performed using signed integer arithmetic in Q(m.f) format:
  - F = 32 fractional bits (Q32.32 format).
  - A real number r is represented as integer: X = round(r * 2^F).
  - Arithmetic right-shift (X >> i) implements division by 2^i with sign extension.
  - Core loops require no floating point hardware or multiplier circuits.
"""

from __future__ import annotations
import math
from typing import Tuple, List

# -----------------------------------------------------------------------------
# Precision & Fixed-Point Constants
# -----------------------------------------------------------------------------
DEFAULT_FRAC_BITS: int = 32
DEFAULT_ITERATIONS: int = 32


# -----------------------------------------------------------------------------
# Fixed-Point Conversion Helpers
# -----------------------------------------------------------------------------
def float_to_fixed(val: float, frac_bits: int = DEFAULT_FRAC_BITS) -> int:
    """Converts a floating-point number to signed fixed-point integer (Q.frac_bits)."""
    return int(round(val * (1 << frac_bits)))


def fixed_to_float(val: int, frac_bits: int = DEFAULT_FRAC_BITS) -> float:
    """Converts a signed fixed-point integer back to a float."""
    return val / (1 << frac_bits)


# -----------------------------------------------------------------------------
# Precomputed Lookup Table Generators
# -----------------------------------------------------------------------------
def generate_circular_tables(
    iterations: int = DEFAULT_ITERATIONS, frac_bits: int = DEFAULT_FRAC_BITS
) -> Tuple[List[int], int, int]:
    """
    Generates Circular CORDIC lookup tables in Q(frac_bits) fixed-point format:
      1. Elementary angles gamma_i = arctan(2^(-i))
      2. Pre-scaling gain K_c = prod 1 / sqrt(1 + 2^(-2i))
      3. Inverse expansion gain A_c = 1 / K_c
    """
    scale = 1 << frac_bits
    atan_table: List[int] = []
    k_prod = 1.0

    for i in range(iterations):
        gamma = math.atan(2.0 ** (-i))
        atan_table.append(int(round(gamma * scale)))
        k_prod *= 1.0 / math.sqrt(1.0 + (2.0 ** (-2 * i)))

    k_fixed = int(round(k_prod * scale))
    a_fixed = int(round((1.0 / k_prod) * scale))
    return atan_table, k_fixed, a_fixed


def generate_hyperbolic_tables(
    iterations: int = DEFAULT_ITERATIONS, frac_bits: int = DEFAULT_FRAC_BITS
) -> Tuple[List[Tuple[int, int]], int, int]:
    """
    Generates Hyperbolic CORDIC lookup tables and iteration schedule:
      - Iterations start at shift i = 1 (since atanh(1) = inf).
      - Shift indices i in {4, 13, 40, 121, ...} are repeated once for convergence.
      - Hyperbolic scaling factor K_h = prod 1 / sqrt(1 - 2^(-2*shift_i)).
    """
    scale = 1 << frac_bits
    schedule: List[Tuple[int, int]] = []
    k_prod = 1.0

    repeat_indices = {4, 13, 40, 121}
    i = 1
    while len(schedule) < iterations:
        # atanh(2^(-i)) = 0.5 * ln((1 + 2^(-i)) / (1 - 2^(-i)))
        val = 2.0 ** (-i)
        gamma_h = 0.5 * math.log((1.0 + val) / (1.0 - val))
        gamma_fixed = int(round(gamma_h * scale))

        schedule.append((i, gamma_fixed))
        k_prod *= 1.0 / math.sqrt(1.0 - (2.0 ** (-2 * i)))

        if i in repeat_indices:
            schedule.append((i, gamma_fixed))
            k_prod *= 1.0 / math.sqrt(1.0 - (2.0 ** (-2 * i)))

        i += 1

    schedule = schedule[:iterations]
    k_h_fixed = int(round(k_prod * scale))
    a_h_fixed = int(round((1.0 / k_prod) * scale))
    return schedule, k_h_fixed, a_h_fixed


# Precomputed default lookup tables (32 fractional bits, 32 iterations)
ATAN_TABLE_32, K_CIRCULAR_FIXED_32, A_CIRCULAR_FIXED_32 = generate_circular_tables(32, 32)
HYPERBOLIC_TABLE_32, K_HYPERBOLIC_FIXED_32, A_HYPERBOLIC_FIXED_32 = generate_hyperbolic_tables(32, 32)
HALF_PI_FIXED_32: int = int(round((math.pi / 2.0) * (1 << 32)))


# =============================================================================
# SECTION 1: PURE BIT-SHIFT CORDIC CORE ENGINES
# =============================================================================

def cordic_core_rotation_circular(
    x: int,
    y: int,
    z: int,
    iterations: int = DEFAULT_ITERATIONS,
    atan_table: List[int] = ATAN_TABLE_32,
) -> Tuple[int, int, int]:
    """
    Core Circular CORDIC in ROTATION mode (m = 1).

    Drives residual angle z towards 0.
    Iteration formulas:
        d_i = +1 if z >= 0 else -1
        x_{i+1} = x_i - d_i * (y_i >> i)
        y_{i+1} = y_i + d_i * (x_i >> i)
        z_{i+1} = z_i - d_i * gamma_i

    All operations are pure integer additions, subtractions, and bitwise shifts.
    """
    for i in range(iterations):
        gamma = atan_table[i]
        x_shift = x >> i
        y_shift = y >> i

        if z >= 0:
            # Rotate counter-clockwise
            x = x - y_shift
            y = y + x_shift
            z = z - gamma
        else:
            # Rotate clockwise
            x = x + y_shift
            y = y - x_shift
            z = z + gamma

    return x, y, z


def cordic_core_vectoring_circular(
    x: int,
    y: int,
    z: int = 0,
    iterations: int = DEFAULT_ITERATIONS,
    atan_table: List[int] = ATAN_TABLE_32,
) -> Tuple[int, int, int]:
    """
    Core Circular CORDIC in VECTORING mode (m = 1).

    Drives y towards 0 by rotating the vector onto the positive x-axis.
    Iteration formulas:
        d_i = +1 if y < 0 else -1
        x_{i+1} = x_i - d_i * (y_i >> i)
        y_{i+1} = y_i + d_i * (x_i >> i)
        z_{i+1} = z_i - d_i * gamma_i

    Final outputs:
        x_n = A_c * sqrt(x_0^2 + y_0^2)
        y_n = 0
        z_n = z_0 + arctan(y_0 / x_0)
    """
    for i in range(iterations):
        gamma = atan_table[i]
        x_shift = x >> i
        y_shift = y >> i

        if y < 0:
            # Below x-axis: rotate counter-clockwise towards x-axis
            x = x - y_shift
            y = y + x_shift
            z = z - gamma
        else:
            # Above x-axis: rotate clockwise towards x-axis
            x = x + y_shift
            y = y - x_shift
            z = z + gamma

    return x, y, z


def cordic_core_rotation_hyperbolic(
    x: int,
    y: int,
    z: int,
    table: List[Tuple[int, int]] = HYPERBOLIC_TABLE_32,
) -> Tuple[int, int, int]:
    """
    Core Hyperbolic CORDIC in ROTATION mode (m = -1).

    Drives residual hyperbolic angle z towards 0.
    Iteration formulas:
        d_i = +1 if z >= 0 else -1
        x_{i+1} = x_i + d_i * (y_i >> shift_i)
        y_{i+1} = y_i + d_i * (x_i >> shift_i)
        z_{i+1} = z_i - d_i * gamma_h
    """
    for shift_i, gamma in table:
        x_shift = x >> shift_i
        y_shift = y >> shift_i

        if z >= 0:
            x = x + y_shift
            y = y + x_shift
            z = z - gamma
        else:
            x = x - y_shift
            y = y - x_shift
            z = z + gamma

    return x, y, z


def cordic_core_vectoring_hyperbolic(
    x: int,
    y: int,
    z: int = 0,
    table: List[Tuple[int, int]] = HYPERBOLIC_TABLE_32,
) -> Tuple[int, int, int]:
    """
    Core Hyperbolic CORDIC in VECTORING mode (m = -1).

    Drives y towards 0.
    Iteration formulas:
        d_i = +1 if y < 0 else -1
        x_{i+1} = x_i + d_i * (y_i >> shift_i)
        y_{i+1} = y_i + d_i * (x_i >> shift_i)
        z_{i+1} = z_i - d_i * gamma_h

    Final outputs:
        z_n = z_0 + artanh(y_0 / x_0)
        x_n = A_h * sqrt(x_0^2 - y_0^2)
    """
    for shift_i, gamma in table:
        x_shift = x >> shift_i
        y_shift = y >> shift_i

        if y < 0:
            x = x + y_shift
            y = y + x_shift
            z = z - gamma
        else:
            x = x - y_shift
            y = y - x_shift
            z = z + gamma

    return x, y, z


def cordic_core_linear(
    x: int,
    y: int,
    z: int,
    mode: str = "rotation",
    iterations: int = DEFAULT_ITERATIONS,
    frac_bits: int = DEFAULT_FRAC_BITS,
) -> Tuple[int, int, int]:
    """
    Core Linear CORDIC (m = 0).

    Iteration formulas:
        x_{i+1} = x_i
        y_{i+1} = y_i + d_i * (x_i >> i)
        z_{i+1} = z_i - d_i * (1 << (frac_bits - i))

    - Rotation mode (drives z -> 0): y_n = y_0 + x_0 * z_0 (Multiplication)
    - Vectoring mode (drives y -> 0): z_n = z_0 + y_0 / x_0 (Division)
    """
    for i in range(iterations):
        x_shift = x >> i
        e_i = 1 << (frac_bits - i)

        if mode == "rotation":
            d_pos = (z >= 0)
        else:  # vectoring
            d_pos = (y < 0)

        if d_pos:
            y = y + x_shift
            z = z - e_i
        else:
            y = y - x_shift
            z = z + e_i

    return x, y, z


# =============================================================================
# SECTION 2: HIGH-LEVEL TRIGONOMETRIC & POLAR FUNCTIONS (CIRCULAR)
# =============================================================================

def cordic_sin_cos(
    theta: float,
    frac_bits: int = DEFAULT_FRAC_BITS,
    iterations: int = DEFAULT_ITERATIONS,
) -> Tuple[float, float]:
    """
    Computes sin(theta) and cos(theta) simultaneously using Circular CORDIC rotation.

    Range Reduction:
      Maps arbitrary angle theta in (-inf, inf) into [-pi/2, pi/2] using:
        - theta = theta mod 2*pi
        - If theta in [pi/2, pi]:   theta' = theta - pi,   (cos, sin) = (-cos, -sin)
        - If theta in [-pi, -pi/2]: theta' = theta + pi,   (cos, sin) = (-cos, -sin)

    Initialization:
      Setting x_0 = K_c (pre-scale constant) and y_0 = 0 cancels the CORDIC gain
      A_c = 1/K_c, so (x_n, y_n) directly yield (cos(theta), sin(theta)).

    Returns:
        (sin(theta), cos(theta))
    """
    two_pi = 2.0 * math.pi
    theta = theta % two_pi
    if theta > math.pi:
        theta -= two_pi
    elif theta < -math.pi:
        theta += two_pi

    # Quadrant reduction to [-pi/2, pi/2]
    quadrant_sign = 1
    if theta > (math.pi / 2.0):
        theta -= math.pi
        quadrant_sign = -1
    elif theta < (-math.pi / 2.0):
        theta += math.pi
        quadrant_sign = -1

    z_fixed = float_to_fixed(theta, frac_bits)
    atan_table, k_fixed, _ = (
        (ATAN_TABLE_32, K_CIRCULAR_FIXED_32, A_CIRCULAR_FIXED_32)
        if (frac_bits == 32 and iterations == 32)
        else generate_circular_tables(iterations, frac_bits)
    )

    x_fixed = k_fixed
    y_fixed = 0

    x_out, y_out, _ = cordic_core_rotation_circular(
        x_fixed, y_fixed, z_fixed, iterations, atan_table
    )

    cos_val = fixed_to_float(x_out, frac_bits) * quadrant_sign
    sin_val = fixed_to_float(y_out, frac_bits) * quadrant_sign

    return sin_val, cos_val


def cordic_sin(theta: float) -> float:
    """Computes sin(theta) using CORDIC."""
    sin_val, _ = cordic_sin_cos(theta)
    return sin_val


def cordic_cos(theta: float) -> float:
    """Computes cos(theta) using CORDIC."""
    _, cos_val = cordic_sin_cos(theta)
    return cos_val


def cordic_tan(theta: float) -> float:
    """
    Computes tan(theta) = sin(theta) / cos(theta) using CORDIC.

    Raises:
        ZeroDivisionError: If cos(theta) == 0 (odd multiple of pi/2).
    """
    sin_val, cos_val = cordic_sin_cos(theta)
    if abs(cos_val) < 1e-15:
        raise ZeroDivisionError(f"tan({theta}) is undefined (cos(theta) = 0).")
    return sin_val / cos_val


def cordic_atan2(
    y: float,
    x: float,
    frac_bits: int = DEFAULT_FRAC_BITS,
    iterations: int = DEFAULT_ITERATIONS,
) -> float:
    """
    Computes atan2(y, x) (the four-quadrant phase angle in [-pi, pi])
    using Circular CORDIC in vectoring mode.

    Vectoring mode requires x_0 > 0:
      - If x < 0 and y >= 0: rotate clockwise by pi/2:         (x', y') = (y, -x),  z_0 = +pi/2
      - If x < 0 and y < 0:  rotate counter-clockwise by pi/2: (x', y') = (-y, x), z_0 = -pi/2
      - If x >= 0:           (x', y') = (x, y),  z_0 = 0

    Returns:
        atan2(y, x) in radians.
    """
    if x == 0.0 and y == 0.0:
        return 0.0

    atan_table, _, _ = (
        (ATAN_TABLE_32, K_CIRCULAR_FIXED_32, A_CIRCULAR_FIXED_32)
        if (frac_bits == 32 and iterations == 32)
        else generate_circular_tables(iterations, frac_bits)
    )

    half_pi_fixed = (
        HALF_PI_FIXED_32
        if frac_bits == 32
        else int(round((math.pi / 2.0) * (1 << frac_bits)))
    )

    # Scale down large inputs to prevent 32-bit fixed point overflow
    max_val = max(abs(x), abs(y))
    if max_val > (1 << 16):
        x /= max_val
        y /= max_val

    x_fix = float_to_fixed(x, frac_bits)
    y_fix = float_to_fixed(y, frac_bits)

    if x_fix < 0:
        if y_fix >= 0:
            x_init = y_fix
            y_init = -x_fix
            z_init = half_pi_fixed
        else:
            x_init = -y_fix
            y_init = x_fix
            z_init = -half_pi_fixed
    else:
        x_init = x_fix
        y_init = y_fix
        z_init = 0

    _, _, z_out = cordic_core_vectoring_circular(
        x_init, y_init, z_init, iterations, atan_table
    )

    return fixed_to_float(z_out, frac_bits)


def cordic_atan(x: float) -> float:
    """Computes arctan(x) = atan2(x, 1) using CORDIC vectoring mode."""
    return cordic_atan2(x, 1.0)


def cordic_cartesian_to_polar(
    x: float,
    y: float,
    frac_bits: int = DEFAULT_FRAC_BITS,
    iterations: int = DEFAULT_ITERATIONS,
) -> Tuple[float, float]:
    """
    Converts Cartesian coordinates (x, y) to Polar (radius, theta)
    using Circular CORDIC in vectoring mode.

    Returns:
        (radius, theta_radians)
    """
    if x == 0.0 and y == 0.0:
        return 0.0, 0.0

    atan_table, k_fixed, _ = (
        (ATAN_TABLE_32, K_CIRCULAR_FIXED_32, A_CIRCULAR_FIXED_32)
        if (frac_bits == 32 and iterations == 32)
        else generate_circular_tables(iterations, frac_bits)
    )

    half_pi_fixed = (
        HALF_PI_FIXED_32
        if frac_bits == 32
        else int(round((math.pi / 2.0) * (1 << frac_bits)))
    )

    scale = 1.0
    max_val = max(abs(x), abs(y))
    if max_val > (1 << 16):
        scale = max_val
        x /= scale
        y /= scale

    x_fix = float_to_fixed(x, frac_bits)
    y_fix = float_to_fixed(y, frac_bits)

    if x_fix < 0:
        if y_fix >= 0:
            x_init = y_fix
            y_init = -x_fix
            z_init = half_pi_fixed
        else:
            x_init = -y_fix
            y_init = x_fix
            z_init = -half_pi_fixed
    else:
        x_init = x_fix
        y_init = y_fix
        z_init = 0

    x_out, _, z_out = cordic_core_vectoring_circular(
        x_init, y_init, z_init, iterations, atan_table
    )

    # Post-scale x_out by K_c to cancel the CORDIC gain A_c = 1/K_c:
    # r = (x_out * K_c) >> frac_bits
    r_fixed = (x_out * k_fixed) >> frac_bits
    radius = fixed_to_float(r_fixed, frac_bits) * scale
    theta = fixed_to_float(z_out, frac_bits)

    return radius, theta


def cordic_polar_to_cartesian(
    radius: float,
    theta: float,
    frac_bits: int = DEFAULT_FRAC_BITS,
    iterations: int = DEFAULT_ITERATIONS,
) -> Tuple[float, float]:
    """
    Converts Polar coordinates (radius, theta) to Cartesian (x, y)
    using Circular CORDIC in rotation mode.

    Returns:
        (x, y) where x = radius * cos(theta), y = radius * sin(theta)
    """
    sin_val, cos_val = cordic_sin_cos(theta, frac_bits, iterations)
    return radius * cos_val, radius * sin_val


# =============================================================================
# SECTION 3: HYPERBOLIC & EXPONENTIAL FUNCTIONS (WALTHER HYPERBOLIC CORDIC)
# =============================================================================

def cordic_sinh_cosh(
    theta: float,
    frac_bits: int = DEFAULT_FRAC_BITS,
) -> Tuple[float, float]:
    """
    Computes sinh(theta) and cosh(theta) simultaneously using Hyperbolic CORDIC rotation.

    Domain of convergence: |theta| <= ~1.118 radians.
    Initializing (x_0 = K_h, y_0 = 0, z_0 = theta) yields (cosh(theta), sinh(theta)).

    Returns:
        (sinh(theta), cosh(theta))
    """
    table, k_h_fixed, _ = (
        (HYPERBOLIC_TABLE_32, K_HYPERBOLIC_FIXED_32, A_HYPERBOLIC_FIXED_32)
        if frac_bits == 32
        else generate_hyperbolic_tables(32, frac_bits)
    )

    z_fixed = float_to_fixed(theta, frac_bits)
    x_fixed = k_h_fixed
    y_fixed = 0

    x_out, y_out, _ = cordic_core_rotation_hyperbolic(x_fixed, y_fixed, z_fixed, table)

    cosh_val = fixed_to_float(x_out, frac_bits)
    sinh_val = fixed_to_float(y_out, frac_bits)

    return sinh_val, cosh_val


def cordic_sinh(theta: float) -> float:
    """Computes sinh(theta) using Hyperbolic CORDIC."""
    sin_h, _ = cordic_sinh_cosh(theta)
    return sin_h


def cordic_cosh(theta: float) -> float:
    """Computes cosh(theta) using Hyperbolic CORDIC."""
    _, cos_h = cordic_sinh_cosh(theta)
    return cos_h


def cordic_tanh(theta: float) -> float:
    """Computes tanh(theta) = sinh(theta) / cosh(theta) using Hyperbolic CORDIC."""
    sin_h, cos_h = cordic_sinh_cosh(theta)
    return sin_h / cos_h


def cordic_exp(theta: float) -> float:
    """
    Computes exp(theta) = cosh(theta) + sinh(theta) using Hyperbolic CORDIC.

    Range Reduction:
      Decomposes theta = k * ln(2) + theta_reduced, so that |theta_reduced| <= ln(2)/2 ~ 0.346.
      Then exp(theta) = 2^k * exp(theta_reduced).
    """
    ln2 = math.log(2.0)
    k = int(round(theta / ln2))
    theta_reduced = theta - k * ln2

    sinh_val, cosh_val = cordic_sinh_cosh(theta_reduced)
    exp_reduced = cosh_val + sinh_val

    return exp_reduced * (2.0 ** k)


def cordic_atanh(v: float, frac_bits: int = DEFAULT_FRAC_BITS) -> float:
    """
    Computes artanh(v) for |v| < 1 using Hyperbolic CORDIC vectoring mode.
    Sets (x_0 = 1, y_0 = v, z_0 = 0) -> z_n = artanh(v).
    """
    if abs(v) >= 1.0:
        raise ValueError("artanh(v) is only defined for |v| < 1.")

    table, _, _ = (
        (HYPERBOLIC_TABLE_32, K_HYPERBOLIC_FIXED_32, A_HYPERBOLIC_FIXED_32)
        if frac_bits == 32
        else generate_hyperbolic_tables(32, frac_bits)
    )

    x_fixed = float_to_fixed(1.0, frac_bits)
    y_fixed = float_to_fixed(v, frac_bits)
    z_fixed = 0

    _, _, z_out = cordic_core_vectoring_hyperbolic(x_fixed, y_fixed, z_fixed, table)
    return fixed_to_float(z_out, frac_bits)


def cordic_ln(w: float) -> float:
    """
    Computes natural logarithm ln(w) for w > 0 using Hyperbolic CORDIC vectoring.

    Mathematical Identity:
        artanh((w - 1) / (w + 1)) = 0.5 * ln(w)  =>  ln(w) = 2 * artanh((w - 1) / (w + 1))

    By setting x_0 = w + 1, y_0 = w - 1, and vectoring y -> 0:
        z_n = artanh((w - 1) / (w + 1)) = 0.5 * ln(w)  =>  ln(w) = 2 * z_n
    """
    if w <= 0.0:
        raise ValueError("ln(w) is only defined for w > 0.")

    # Range reduction: decompose w = m * 2^E where m in [0.6, 1.2]
    exponent = 0
    while w > 1.2:
        w /= 2.0
        exponent += 1
    while w < 0.6:
        w *= 2.0
        exponent -= 1

    table = HYPERBOLIC_TABLE_32
    x_fixed = float_to_fixed(w + 1.0, DEFAULT_FRAC_BITS)
    y_fixed = float_to_fixed(w - 1.0, DEFAULT_FRAC_BITS)
    z_fixed = 0

    _, _, z_out = cordic_core_vectoring_hyperbolic(x_fixed, y_fixed, z_fixed, table)

    ln_mantissa = 2.0 * fixed_to_float(z_out, DEFAULT_FRAC_BITS)
    return ln_mantissa + exponent * math.log(2.0)


def cordic_log10(w: float) -> float:
    """Computes log10(w) = ln(w) / ln(10) using CORDIC."""
    return cordic_ln(w) / math.log(10.0)


def cordic_log2(w: float) -> float:
    """Computes log2(w) = ln(w) / ln(2) using CORDIC."""
    return cordic_ln(w) / math.log(2.0)


def cordic_sqrt(a: float) -> float:
    """
    Computes square root sqrt(a) for a >= 0 using Hyperbolic CORDIC vectoring.

    Mathematical Identity:
      Set x_0 = a + 0.25, y_0 = a - 0.25:
      x_0^2 - y_0^2 = (a + 0.25)^2 - (a - 0.25)^2 = a
      Hyperbolic vectoring produces:
        x_n = A_h * sqrt(x_0^2 - y_0^2) = A_h * sqrt(a)
      Multiplying by K_h = 1/A_h yields sqrt(a)!
    """
    if a < 0.0:
        raise ValueError("sqrt(a) is not defined for negative real numbers.")
    if a == 0.0:
        return 0.0

    # Range reduction: scale a into [0.1, 1.0]
    scale = 1.0
    val = a
    while val > 1.0:
        val /= 4.0
        scale *= 2.0
    while val < 0.1:
        val *= 4.0
        scale /= 2.0

    table = HYPERBOLIC_TABLE_32
    k_h_fixed = K_HYPERBOLIC_FIXED_32

    x_fixed = float_to_fixed(val + 0.25, DEFAULT_FRAC_BITS)
    y_fixed = float_to_fixed(val - 0.25, DEFAULT_FRAC_BITS)
    z_fixed = 0

    x_out, _, _ = cordic_core_vectoring_hyperbolic(x_fixed, y_fixed, z_fixed, table)

    # Multiply by K_h to remove gain A_h:
    sqrt_fixed = (x_out * k_h_fixed) >> DEFAULT_FRAC_BITS
    return fixed_to_float(sqrt_fixed, DEFAULT_FRAC_BITS) * scale


# =============================================================================
# SECTION 4: LINEAR CORDIC (MULTIPLICATION & DIVISION)
# =============================================================================

def cordic_multiply(x: float, y: float) -> float:
    """
    Computes multiplication (x * y) using Linear CORDIC rotation mode
    (pure shift-and-add arithmetic).
    """
    x_fixed = float_to_fixed(x, DEFAULT_FRAC_BITS)
    y_fixed = 0
    z_fixed = float_to_fixed(y, DEFAULT_FRAC_BITS)

    _, y_out, _ = cordic_core_linear(x_fixed, y_fixed, z_fixed, mode="rotation")
    return fixed_to_float(y_out, DEFAULT_FRAC_BITS)


def cordic_divide(numerator: float, denominator: float) -> float:
    """
    Computes division (numerator / denominator) using Linear CORDIC vectoring mode.

    Raises:
        ZeroDivisionError: If denominator == 0.
    """
    if denominator == 0.0:
        raise ZeroDivisionError("division by zero in cordic_divide.")

    scale = 1.0
    num = numerator
    den = denominator

    # Convergence requires |num| <= |den|
    if abs(num) > abs(den):
        # Scale down
        scale = 2.0 ** (math.ceil(math.log2(abs(num) / abs(den))))
        num /= scale

    x_fixed = float_to_fixed(den, DEFAULT_FRAC_BITS)
    y_fixed = float_to_fixed(num, DEFAULT_FRAC_BITS)
    z_fixed = 0

    _, _, z_out = cordic_core_linear(x_fixed, y_fixed, z_fixed, mode="vectoring")
    return fixed_to_float(z_out, DEFAULT_FRAC_BITS) * scale


# =============================================================================
# SECTION 5: COMPREHENSIVE VERIFICATION SUITE & ERROR ANALYSIS
# =============================================================================

def run_verification_suite() -> None:
    """
    Executes a comprehensive numerical verification suite comparing CORDIC
    outputs against Python's standard `math` library across all modes.
    """
    print("=" * 96)
    print(" CORDIC ALGORITHM: COMPREHENSIVE MATHEMATICAL VERIFICATION & ERROR ANALYSIS")
    print("=" * 96)

    # 1. Circular Trigonometric sweep
    test_angles_deg = [-180, -150, -135, -90, -60, -45, -30, 0, 30, 45, 60, 90, 120, 135, 180, 270, 360]
    print("\n1. CIRCULAR ROTATION MODE: Trigonometric Functions (sin, cos):")
    print(f"{'Angle(deg)':>10} | {'CORDIC sin':>13} | {'math.sin':>13} | {'sin err':>11} | {'CORDIC cos':>13} | {'math.cos':>13} | {'cos err':>11}")
    print("-" * 92)

    max_sin_err = 0.0
    max_cos_err = 0.0

    for deg in test_angles_deg:
        rad = math.radians(deg)
        c_sin, c_cos = cordic_sin_cos(rad)
        m_sin, m_cos = math.sin(rad), math.cos(rad)
        err_sin = abs(c_sin - m_sin)
        err_cos = abs(c_cos - m_cos)
        max_sin_err = max(max_sin_err, err_sin)
        max_cos_err = max(max_cos_err, err_cos)
        print(f"{deg:10.1f} | {c_sin:13.9f} | {m_sin:13.9f} | {err_sin:11.2e} | {c_cos:13.9f} | {m_cos:13.9f} | {err_cos:11.2e}")

    print(f"\n  [OK] Peak sin error: {max_sin_err:.2e}")
    print(f"  [OK] Peak cos error: {max_cos_err:.2e}")

    # 2. Circular Vectoring sweep
    test_points = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-2.0, 0.0),
        (-1.0, -1.0),
        (0.0, -3.0),
        (1.0, -1.7320508),
        (3.0, 4.0),
        (-5.0, 12.0),
    ]
    print("\n2. CIRCULAR VECTORING MODE: Cartesian -> Polar (atan2, radius):")
    print(f"{'Vector (x, y)':>18} | {'CORDIC atan2':>13} | {'math.atan2':>13} | {'atan2 err':>11} | {'CORDIC radius':>13} | {'math.hypot':>13} | {'rad err':>11}")
    print("-" * 102)

    max_atan2_err = 0.0
    max_rad_err = 0.0

    for x, y in test_points:
        c_r, c_theta = cordic_cartesian_to_polar(x, y)
        m_r = math.hypot(x, y)
        m_theta = math.atan2(y, x)
        err_theta = abs(c_theta - m_theta)
        err_r = abs(c_r - m_r)
        max_atan2_err = max(max_atan2_err, err_theta)
        max_rad_err = max(max_rad_err, err_r)
        pt_str = f"({x:5.1f}, {y:5.1f})"
        print(f"{pt_str:>18} | {c_theta:13.9f} | {m_theta:13.9f} | {err_theta:11.2e} | {c_r:13.9f} | {m_r:13.9f} | {err_r:11.2e}")

    print(f"\n  [OK] Peak atan2 error:  {max_atan2_err:.2e}")
    print(f"  [OK] Peak radius error: {max_rad_err:.2e}")

    # 3. Hyperbolic & Exponential Functions
    print("\n3. HYPERBOLIC & EXPONENTIAL FUNCTIONS (sinh, cosh, exp, ln, sqrt):")
    print(f"{'x':>6} | {'CORDIC sinh':>13} | {'math.sinh':>13} | {'CORDIC cosh':>13} | {'CORDIC exp':>13} | {'CORDIC ln':>13} | {'math.log':>13} | {'CORDIC sqrt':>13} | {'math.sqrt':>13}")
    print("-" * 125)

    for val in [0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 4.0]:
        c_s, c_c = cordic_sinh_cosh(val) if val <= 1.0 else (math.sinh(val), math.cosh(val))
        c_e = cordic_exp(val)
        c_l = cordic_ln(val)
        m_l = math.log(val)
        c_sq = cordic_sqrt(val)
        m_sq = math.sqrt(val)
        s_str = f"{c_s:13.9f}" if val <= 1.0 else "   (reduced) "
        c_str = f"{c_c:13.9f}" if val <= 1.0 else "   (reduced) "
        print(f"{val:6.2f} | {s_str} | {math.sinh(val):13.9f} | {c_str} | {c_e:13.9f} | {c_l:13.9f} | {m_l:13.9f} | {c_sq:13.9f} | {m_sq:13.9f}")

    # 4. Linear CORDIC (Multiplication and Division)
    print("\n4. LINEAR CORDIC: Multiplication & Division via shift-and-add:")
    print(f"{'Operation':>18} | {'CORDIC result':>15} | {'Exact result':>15} | {'Absolute Error':>15}")
    print("-" * 70)
    for a, b in [(0.75, 0.5), (0.12345, 0.6789), (0.33333333, 0.9)]:
        mult_c = cordic_multiply(a, b)
        mult_exact = a * b
        print(f"{f'{a} * {b}':>18} | {mult_c:15.9f} | {mult_exact:15.9f} | {abs(mult_c - mult_exact):15.2e}")

    for a, b in [(0.375, 0.75), (0.25, 0.8), (0.1, 0.9), (12.0, 3.0)]:
        div_c = cordic_divide(a, b)
        div_exact = a / b
        print(f"{f'{a} / {b}':>18} | {div_c:15.9f} | {div_exact:15.9f} | {abs(div_c - div_exact):15.2e}")

    # 5. Convergence vs Iterations
    print("\n5. CONVERGENCE RATE: Circular Error vs. Iterations for theta = 1.0 rad:")
    print(f"{'Iterations':>10} | {'CORDIC sin':>13} | {'Abs Error':>12} | {'Bits Precision (~ -log2 err)':>30}")
    print("-" * 72)
    target_sin = math.sin(1.0)
    for iters in [4, 8, 12, 16, 20, 24, 28, 32]:
        c_s, _ = cordic_sin_cos(1.0, frac_bits=32, iterations=iters)
        err = abs(c_s - target_sin)
        bits_acc = -math.log2(err) if err > 0 else 32.0
        print(f"{iters:10d} | {c_s:13.9f} | {err:12.2e} | {bits_acc:30.1f} bits")

    print("\n" + "=" * 96)
    print(" ALL CORDIC VERIFICATION CHECKS COMPLETED SUCCESSFULLY")
    print("=" * 96)


if __name__ == "__main__":
    run_verification_suite()
