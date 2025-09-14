// AeroInterp.cpp  -- corrected signedness & small safety checks
#include "AeroInterp.h"
#include "aero_table_takeoff_full.h" // generated header you already made

#include <cstddef>
#include <cmath>
#include <limits>
#include <cstdio>

// Tolerances
static constexpr double AILERON_TOL = 1e-6;      // allowed deviation for Aileron_left (was fixed to 40)
static constexpr double LAMBDA_EPS   = 1e-9;     // tolerance when checking lambda bounds
static constexpr double SUM_EPS      = 1e-6;     // tolerance for sum(lambda) ~ 1

// Safety: maximum (D+1) we expect (tuned for small D). Increase if you have larger D.
static constexpr std::size_t MAX_DPLUS1 = 16;

// Mapping original 6 inputs -> reduced inputs used by triangulation
// The original column order is: Throttle(0), Aileron_left(1), Elevator(2), Rudder(3), Alpha(4), Beta(5)
// The triangulation used the reduced columns in this order (as recorded in header comment):
//   Throttle, Elevator, Rudder, Alpha, Beta
// So the indices mapping are: 0, 2, 3, 4, 5
static const std::size_t ORIG_TO_REDUCED_IDX[] = { 0, 2, 3, 4, 5 };

bool interpolate_aerocoeffs(double throttle,
                            double aileron_left,
                            double elevator,
                            double rudder,
                            double alpha,
                            double beta,
                            double aerocoeffs[])
{
    // Basic consistency checks
    const std::size_t reducedD = 5;
    const std::size_t expected_mapping_len = sizeof(ORIG_TO_REDUCED_IDX)/sizeof(ORIG_TO_REDUCED_IDX[0]);
    if (reducedD != expected_mapping_len) {
        // Dimension mismatch between generated header and mapping in this code.
        std::fprintf(stderr, "AeroInterp: dimension mismatch AERO_D=%zu vs mapping=%zu\n", reducedD, expected_mapping_len);
        return false;
    }

    // Check bounds for small arrays
    const std::size_t Dplus1 = reducedD + 1;
    if (Dplus1 > MAX_DPLUS1) {
        std::fprintf(stderr, "AeroInterp: D+1 (%zu) exceeds MAX_DPLUS1 (%zu)\n", Dplus1, MAX_DPLUS1);
        return false;
    }

    // Check aileron is the fixed value used when building the table (40)
    if (std::fabs(aileron_left - 40.0) > AILERON_TOL) {
        // Fail gracefully if aileron is not 40 in 5-D table mode
        // For 6-D table (full), just accept any aileron value.
        return false;
    }

    // Build reduced point
    double orig[6] = { throttle, aileron_left, elevator, rudder, alpha, beta };
    double p_reduced_static[MAX_DPLUS1]; // reuse stack array
    for (std::size_t i = 0; i < reducedD; ++i) {
        p_reduced_static[i] = orig[ ORIG_TO_REDUCED_IDX[i] ];
    }

    // rhs = [p_reduced; 1.0]
    double rhs_static[MAX_DPLUS1];
    for (std::size_t i = 0; i < reducedD; ++i) rhs_static[i] = p_reduced_static[i];
    rhs_static[reducedD] = 1.0;

    // Iterate over simplices (AERO_FULL_NSIMPLICES is size_t, as you have full 6D table)
    for (std::size_t s = 0; s < AERO_FULL_NSIMPLICES; ++s) {
        const double* Tinv_flat = AERO_FULL_TINV[s]; // flattened row-major (D+1)x(D+1)

        // compute lambdas = Tinv * rhs
        double lambdas_static[MAX_DPLUS1];
        bool any_nan = false;
        for (std::size_t i = 0; i < Dplus1; ++i) {
            double acc = 0.0;
            const std::size_t row_off = i * Dplus1;
            for (std::size_t j = 0; j < Dplus1; ++j) {
                acc += Tinv_flat[row_off + j] * rhs_static[j];
            }
            if (!std::isfinite(acc)) { any_nan = true; break; }
            lambdas_static[i] = acc;
        }
        if (any_nan) continue;

        // Check lambda bounds and sum
        bool all_in = true;
        double sum_lambda = 0.0;
        for (std::size_t i = 0; i < Dplus1; ++i) {
            double L = lambdas_static[i];
            sum_lambda += L;
            if (L < -LAMBDA_EPS || L > 1.0 + LAMBDA_EPS) { all_in = false; break; }
        }
        if (!all_in) continue;
        if (std::fabs(sum_lambda - 1.0) > SUM_EPS) continue;

        // compute predicted outputs: Y_pred = sum_i lambdas[i] * Y[vertex_i]
        // safeguard simplex vertex indices and array bounds
        bool corrupt_simplex = false;
        for (std::size_t j = 0; j < (std::size_t)AERO_FULL_M; ++j) aerocoeffs[j] = 0.0;

        for (std::size_t k = 0; k < Dplus1; ++k) {
            int vid = AERO_FULL_SIMPLICES[s][k];
            if (vid < 0 || static_cast<std::size_t>(vid) >= AERO_FULL_N) { corrupt_simplex = true; break; }
            for (std::size_t j = 0; j < (std::size_t)AERO_FULL_M; ++j) {
                aerocoeffs[j] += lambdas_static[k] * AERO_FULL_Y[static_cast<std::size_t>(vid)][j];
            }
        }
        if (corrupt_simplex) continue;

        // success
        return true;
    }

    // no simplex found containing the point
    return false;
}
