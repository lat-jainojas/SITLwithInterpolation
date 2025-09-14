#pragma once
#include <cstddef>

/*
 * interpolate_aerocoeffs
 *
 * Inputs (in this exact order):
 *   throttle, aileron_left, elevator, rudder, alpha, beta
 *
 * Output:
 *   aerocoeffs must be an array of length AERO_M (3). On success it will contain:
 *     aerocoeffs[0] = CL
 *     aerocoeffs[1] = CD
 *     aerocoeffs[2] = CMm_xcg_1.04
 *
 * Return value:
 *   true  -> interpolation succeeded (point found inside convex hull & valid simplex)
 *   false -> interpolation failed (aileron_left outside expected fixed value OR point outside hull)
 *
 * Notes:
 * - This implementation assumes the triangulation header `aero_table_takeoff.h` is available
 *   and that it declares the arrays AERO_X, AERO_Y, AERO_SIMPLICES, AERO_TINV, and the
 *   constants AERO_N, AERO_D, AERO_M, AERO_NSIMPLICES.
 *
 * - The table used for triangulation was built with Aileron_left fixed to 40 (so the
 *   reduced input dimension is D=5). This function checks that the passed aileron_left
 *   is close to 40 (tolerance default 1e-6). If not, it returns false.
 *
 * - The function uses the precomputed per-simplex T_inv matrices (AERO_TINV) and
 *   therefore only performs small matrix-vector multiplies at runtime (cheap).
 */
bool interpolate_aerocoeffs(double throttle,
                            double aileron_left,
                            double elevator,
                            double rudder,
                            double alpha,
                            double beta,
                            double aerocoeffs[] /*size >= AERO_M*/);
