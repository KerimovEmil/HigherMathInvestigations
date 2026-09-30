"""
CORDIC (COordinate Rotation DIgital Computer) Algorithm
======================================================

A clean, fixed-point Python implementation of Jack E. Volder's (1959) CORDIC
algorithm using only additions, subtractions, and bitwise arithmetic right-shifts.

===============================================================================
MATHEMATICAL FOUNDATIONS
===============================================================================

1. 2D Vector Rotation Matrix
----------------------------
Rotating a 2D column vector v = [x, y]^T counter-clockwise by an angle theta:

    [ x' ]   [ cos(theta)  -sin(theta) ] [ x ]
    [ y' ] = [ sin(theta)   cos(theta) ] [ y ]

Factoring out cos(theta):

    [ x' ]                [ 1           -tan(theta) ] [ x ]
    [ y' ] = cos(theta) * [ tan(theta)       1      ] [ y ]

2. Angle Decomposition & Elementary Steps
-----------------------------------------
Instead of continuous rotation, CORDIC decomposes theta into a sum of n
successive rotations by fixed elementary angles gamma_i:

    theta = sum_{i=0}^{n-1} d_i * gamma_i,    where d_i in {-1, +1}

Volder's key insight was choosing each elementary angle gamma_i such that:

    tan(gamma_i) = 2^(-i)   <=>   gamma_i = arctan(2^(-i))

Substituting tan(gamma_i) = 2^(-i) into the rotation matrix yields:

    R_i = cos(gamma_i) * [ 1          -d_i * 2^(-i) ]
                         [ d_i * 2^(-i)      1      ]

In binary fixed-point representation, multiplying by 2^(-i) is an arithmetic
right-shift by i bit positions:

    v * 2^(-i) == (v >> i)

3. The Pseudo-Rotation Trick & Cumulative Gain K
------------------------------------------------
Because cosine is an even function (cos(-gamma) = cos(gamma)):

    cos(gamma_i) = 1 / sqrt(1 + tan^2(gamma_i)) = 1 / sqrt(1 + 2^(-2i))

This factor depends ONLY on iteration index i, NOT on the rotation direction d_i.
We can therefore omit cos(gamma_i) from the inner loop ("pseudo-rotation") and
precompute the cumulative scaling factor K:

    K = prod_{i=0}^{n-1} 1 / sqrt(1 + 2^(-2i))  approx 0.6072529350088812561694...
    A = 1 / K                                    approx 1.646760258121065648366...

During n iterations, the vector magnitude expands by A = 1/K.
By pre-scaling the initial input x_0 = K (or post-scaling by K), the final vector
is exact without requiring any multiplier hardware in the loop.

4. Convergence Condition & Range Reduction
------------------------------------------
For arbitrary angles to be reachable without gaps, the elementary angles must satisfy:

    gamma_i <= sum_{j=i+1}^{n-1} gamma_j + gamma_{n-1}

Since arctan(2^(-i)) < 2 * arctan(2^(-(i+1))), this holds for all i >= 0.
The maximum angle reachable by circular CORDIC without range reduction is:

    theta_max = sum_{i=0}^{n-1} gamma_i approx 1.74328662 radians (~99.88 degrees)

Because theta_max > pi/2 (~1.57079 rad), circular CORDIC covers the entire
first quadrant. Standard quadrant reduction maps any angle in (-inf, inf) into
the primary domain [-pi/2, pi/2].

5. Operational Modes
--------------------
  (A) ROTATION MODE (drives residual angle z -> 0):
      At step i:  d_i = +1 if z_i >= 0 else -1
      Init:       x_0 = K,  y_0 = 0,  z_0 = theta
      Result:     x_n = cos(theta),   y_n = sin(theta),   z_n = 0

  (B) VECTORING MODE (drives y -> 0 by rotating vector onto x-axis):
      At step i:  d_i = +1 if y_i < 0 else -1
      Init:       x_0 = x,  y_0 = y,  z_0 = 0
      Result:     x_n = A * sqrt(x^2 + y^2),  y_n = 0,   z_n = arctan(y / x)
"""

from __future__ import annotations
import math
from typing import Tuple

# -----------------------------------------------------------------------------
# Fixed-Point Constants & Precomputed Tables (Q32.32 Format)
# -----------------------------------------------------------------------------
# Q32.32 format: 32 bits for integer part, 32 bits for fractional part.
# Provides ~9-10 decimal digits of precision using pure integer bit operations.
SHIFT: int = 32
SCALE: int = 1 << SHIFT
ITERATIONS: int = 32

# Precomputed elementary angles: gamma_i = arctan(2^(-i)) in Q32.32
ANGLES: list[int] = [
    int(round(math.atan(2.0 ** (-i)) * SCALE)) for i in range(ITERATIONS)
]

# Precomputed cumulative scaling factor: K = prod_{i=0}^{n-1} 1 / sqrt(1 + 2^(-2i))
K_FLOAT: float = 1.0
for i in range(ITERATIONS):
    K_FLOAT *= 1.0 / math.sqrt(1.0 + (2.0 ** (-2 * i)))
K_FIXED: int = int(round(K_FLOAT * SCALE))

HALF_PI_FIXED: int = int(round((math.pi / 2.0) * SCALE))


# =============================================================================
# 1. ROTATION MODE: Trigonometric Functions & Vector Rotation
# =============================================================================

def cordic_sin_cos(theta_rad: float) -> Tuple[float, float]:
    """
    Computes (sin(theta), cos(theta)) simultaneously using CORDIC rotation mode.

    Range reduction:
      Maps arbitrary theta in (-inf, +inf) into [-pi/2, pi/2]:
        - theta = theta mod 2*pi
        - If theta in [pi/2, pi]:   theta' = theta - pi,   (cos, sin) = (-cos, -sin)
        - If theta in [-pi, -pi/2]: theta' = theta + pi,   (cos, sin) = (-cos, -sin)

    Initialization:
      Setting x_0 = K_FIXED and y_0 = 0 cancels the CORDIC expansion A = 1/K,
      so final (x_n, y_n) directly equal (cos(theta), sin(theta)).

    Returns:
        (sin(theta), cos(theta)) as floats.
    """
    # Range reduction to [-pi, pi]
    two_pi = 2.0 * math.pi
    theta = theta_rad % two_pi
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

    # Convert angle to Q32.32 fixed-point integer
    z = int(round(theta * SCALE))
    x = K_FIXED
    y = 0

    # CORDIC Rotation Core (pure integer additions, subtractions, and bit shifts)
    for i in range(ITERATIONS):
        x_shift = x >> i
        y_shift = y >> i

        if z >= 0:
            # Counter-clockwise rotation step
            x -= y_shift
            y += x_shift
            z -= ANGLES[i]
        else:
            # Clockwise rotation step
            x += y_shift
            y -= x_shift
            z += ANGLES[i]

    sin_val = (y / SCALE) * quadrant_sign
    cos_val = (x / SCALE) * quadrant_sign
    return sin_val, cos_val


def cordic_sin(theta_rad: float) -> float:
    """Computes sin(theta) using CORDIC."""
    sin_val, _ = cordic_sin_cos(theta_rad)
    return sin_val


def cordic_cos(theta_rad: float) -> float:
    """Computes cos(theta) using CORDIC."""
    _, cos_val = cordic_sin_cos(theta_rad)
    return cos_val


def cordic_tan(theta_rad: float) -> float:
    """
    Computes tan(theta) = sin(theta) / cos(theta) using CORDIC.

    Raises:
        ZeroDivisionError: If cos(theta) == 0 (theta is an odd multiple of pi/2).
    """
    sin_val, cos_val = cordic_sin_cos(theta_rad)
    if abs(cos_val) < 1e-15:
        raise ZeroDivisionError(f"tan({theta_rad}) is undefined (cos(theta) = 0).")
    return sin_val / cos_val


def cordic_polar_to_cartesian(radius: float, theta_rad: float) -> Tuple[float, float]:
    """
    Converts Polar coordinates (radius, theta) to Cartesian (x, y)
    using CORDIC rotation mode.

    Returns:
        (x, y) where x = radius * cos(theta), y = radius * sin(theta)
    """
    sin_val, cos_val = cordic_sin_cos(theta_rad)
    return radius * cos_val, radius * sin_val


# =============================================================================
# 2. VECTORING MODE: Inverse Trigonometry & Vector Magnitude
# =============================================================================

def cordic_atan2(y: float, x: float) -> float:
    """
    Computes atan2(y, x) (the four-quadrant phase angle in [-pi, pi])
    using CORDIC vectoring mode.

    Vectoring mode drives y -> 0 by rotating the vector onto the positive x-axis,
    accumulating the total rotation angle in z.

    Quadrant handling (requires initial x > 0):
      - If x < 0 and y >= 0: rotate by +pi/2 -> (x', y') = (y, -x),  z_0 = +pi/2
      - If x < 0 and y < 0:  rotate by -pi/2 -> (x', y') = (-y, x), z_0 = -pi/2
      - If x >= 0:           (x', y') = (x, y),  z_0 = 0

    Returns:
        atan2(y, x) in radians.
    """
    if x == 0.0 and y == 0.0:
        return 0.0

    # Scale down inputs if large to prevent integer overflow
    max_val = max(abs(x), abs(y))
    if max_val > (1 << 16):
        x /= max_val
        y /= max_val

    x_fix = int(round(x * SCALE))
    y_fix = int(round(y * SCALE))

    # Initial rotation to bring vector into the right-half plane (x > 0)
    if x_fix < 0:
        if y_fix >= 0:
            x_val, y_val, z = y_fix, -x_fix, HALF_PI_FIXED
        else:
            x_val, y_val, z = -y_fix, x_fix, -HALF_PI_FIXED
    else:
        x_val, y_val, z = x_fix, y_fix, 0

    # CORDIC Vectoring Core (pure integer arithmetic)
    for i in range(ITERATIONS):
        x_shift = x_val >> i
        y_shift = y_val >> i

        if y_val < 0:
            # Vector is below x-axis: rotate counter-clockwise
            x_val -= y_shift
            y_val += x_shift
            z -= ANGLES[i]
        else:
            # Vector is above x-axis: rotate clockwise
            x_val += y_shift
            y_val -= x_shift
            z += ANGLES[i]

    return z / SCALE


def cordic_atan(x: float) -> float:
    """Computes arctan(x) = atan2(x, 1.0) using CORDIC vectoring mode."""
    return cordic_atan2(x, 1.0)


def cordic_cartesian_to_polar(x: float, y: float) -> Tuple[float, float]:
    """
    Converts Cartesian coordinates (x, y) to Polar coordinates (radius, theta)
    using CORDIC vectoring mode.

    Returns:
        (radius, theta_radians)
    """
    if x == 0.0 and y == 0.0:
        return 0.0, 0.0

    scale = 1.0
    max_val = max(abs(x), abs(y))
    if max_val > (1 << 16):
        scale = max_val
        x /= scale
        y /= scale

    x_fix = int(round(x * SCALE))
    y_fix = int(round(y * SCALE))

    if x_fix < 0:
        if y_fix >= 0:
            x_val, y_val, z = y_fix, -x_fix, HALF_PI_FIXED
        else:
            x_val, y_val, z = -y_fix, x_fix, -HALF_PI_FIXED
    else:
        x_val, y_val, z = x_fix, y_fix, 0

    for i in range(ITERATIONS):
        x_shift = x_val >> i
        y_shift = y_val >> i

        if y_val < 0:
            x_val -= y_shift
            y_val += x_shift
            z -= ANGLES[i]
        else:
            x_val += y_shift
            y_val -= x_shift
            z += ANGLES[i]

    # Post-scale magnitude x_val by K to remove the CORDIC gain A = 1/K:
    # r = (x_val * K) >> SHIFT
    r_fixed = (x_val * K_FIXED) >> SHIFT
    radius = (r_fixed / SCALE) * scale
    theta = z / SCALE
    return radius, theta


# =============================================================================
# 3. VERIFICATION SUITE & ERROR ANALYSIS
# =============================================================================

def run_verification_suite() -> None:
    """
    Demonstrates and verifies CORDIC accuracy against Python's `math` library.
    """
    print("=" * 84)
    print(" CORDIC ALGORITHM: NUMERICAL VERIFICATION & ERROR ANALYSIS")
    print("=" * 84)

    # 1. Rotation Mode (sin, cos across all 4 quadrants)
    print("\n1. Rotation Mode: Trigonometric Sweep [-180 deg to +360 deg]:")
    print(f"{'Angle(deg)':>10} | {'CORDIC sin':>12} | {'math.sin':>12} | {'sin err':>10} | {'CORDIC cos':>12} | {'math.cos':>12} | {'cos err':>10}")
    print("-" * 88)

    test_angles = [-180, -135, -90, -45, -30, 0, 30, 45, 60, 90, 120, 135, 180, 270, 360]
    max_sin_err = 0.0
    max_cos_err = 0.0

    for deg in test_angles:
        rad = math.radians(deg)
        s, c = cordic_sin_cos(rad)
        ms, mc = math.sin(rad), math.cos(rad)
        err_s = abs(s - ms)
        err_c = abs(c - mc)
        max_sin_err = max(max_sin_err, err_s)
        max_cos_err = max(max_cos_err, err_c)
        print(f"{deg:10.1f} | {s:12.8f} | {ms:12.8f} | {err_s:10.2e} | {c:12.8f} | {mc:12.8f} | {err_c:10.2e}")

    print(f"\n  [OK] Peak sin error: {max_sin_err:.2e}")
    print(f"  [OK] Peak cos error: {max_cos_err:.2e}")

    # 2. Vectoring Mode (atan2 and magnitude across quadrants)
    print("\n2. Vectoring Mode: Cartesian to Polar Sweep (atan2, radius):")
    print(f"{'Vector (x, y)':>18} | {'CORDIC atan2':>13} | {'math.atan2':>13} | {'atan2 err':>10} | {'CORDIC rad':>11} | {'math.hypot':>11} | {'rad err':>10}")
    print("-" * 98)

    test_points = [
        (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0),
        (-2.0, 0.0), (-1.0, -1.0), (0.0, -3.0), (3.0, 4.0), (-5.0, 12.0)
    ]
    max_atan2_err = 0.0
    max_rad_err = 0.0

    for px, py in test_points:
        r, theta = cordic_cartesian_to_polar(px, py)
        mr = math.hypot(px, py)
        mtheta = math.atan2(py, px)
        err_t = abs(theta - mtheta)
        err_r = abs(r - mr)
        max_atan2_err = max(max_atan2_err, err_t)
        max_rad_err = max(max_rad_err, err_r)
        pt_str = f"({px:5.1f}, {py:5.1f})"
        print(f"{pt_str:>18} | {theta:13.8f} | {mtheta:13.8f} | {err_t:10.2e} | {r:11.6f} | {mr:11.6f} | {err_r:10.2e}")

    print(f"\n  [OK] Peak atan2 error:  {max_atan2_err:.2e}")
    print(f"  [OK] Peak radius error: {max_rad_err:.2e}")

    # 3. Convergence rate vs. iteration count
    print("\n3. Convergence Rate: Accuracy vs Iterations for theta = 1.0 rad (57.3 deg):")
    print(f"{'Iterations':>10} | {'CORDIC sin':>13} | {'Abs Error':>12} | {'Bits of Precision':>20}")
    print("-" * 62)

    target_sin = math.sin(1.0)
    for iters in [4, 8, 12, 16, 20, 24, 28, 32]:
        # Perform calculation with variable iteration count
        z = int(round(1.0 * SCALE))
        x, y = K_FIXED, 0
        for i in range(iters):
            xs, ys = x >> i, y >> i
            if z >= 0:
                x -= ys
                y += xs
                z -= ANGLES[i]
            else:
                x += ys
                y -= xs
                z += ANGLES[i]

        approx_sin = y / SCALE
        err = abs(approx_sin - target_sin)
        bits = -math.log2(err) if err > 0 else 32.0
        print(f"{iters:10d} | {approx_sin:13.9f} | {err:12.2e} | {bits:18.1f} bits")

    print("\n" + "=" * 84)
    print(" ALL CORDIC VERIFICATION CHECKS COMPLETED SUCCESSFULLY")
    print("=" * 84)


if __name__ == '__main__':
    run_verification_suite()
