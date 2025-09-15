/*
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
  very simple plane simulator class. Not aerodynamically accurate,
  just enough to be able to debug control logic for new frame types
*/

// Using the generated header file from barycentric interpolation
#include <sstream>
#include <iomanip>
#include <sys/select.h>
#include <AP_HAL/AP_HAL.h>

#include "AeroInterp.h"
#include "aero_table_takeoff.h"
#include "SIM_Plane.h"
#include<vector>
#include <stdio.h>
#include <AP_Filesystem/AP_Filesystem_config.h>
#include <AP_Filesystem/AP_Filesystem.h>
#include <AP_AHRS/AP_AHRS.h>
#include<stdio.h>
// #include "aero_data.h"
// #include "aero_knn.cpp"
// #include "aero_knn.h"

// The following code is for EQUINOX's SITL Model
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// clamp helper
static inline double clamp_double(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

/*
 * compute_n_from_throttle
 * - Uses expression: n = (19771*throttle - 1896.4)/60
 * - throttle expected in [0,1]. We clamp n to a safe minimum to avoid negative or zero rpm.
 */
static inline double compute_n_from_throttle(double throttle) {
    double n = (19771.0 * throttle - 1896.4) / 60.0;
    // guard against negative or zero spin; choose 1.0 rps as safe minimum
    if (n < 1.0) n = 1.0;
    return n;
}

/*
 * compute_Cmu
 * J = V_inf / (n * D)
 * D is in meters (user specified D = 0.120 m)
 * Cmu = 0.376284/(J^2.05254) + 0.035986
 *
 * We guard small airspeed by forcing a minimum V_inf for numeric stability.
 */
static inline double compute_Cmu(double V_inf, double n, double D = 0.120) {
    // avoid division by zero / extremely small J
    double Vsafe = (V_inf < 0.1) ? 0.1 : V_inf;       // 0.1 m/s min
    double J = Vsafe / (n * D);
    if (J < 1e-6) J = 1e-6;
    double Cmu = 0.376284 / pow(J, 2.05254) + 0.035986;
    return Cmu;
}

/*
 * compute_CL_from_model
 * Implements exact formula:
 *   CL = CL0 + CL_alpha*Alpha + CL_del_e*del_e + CL_del_f*del_f + CL_Cmu*Cmu + CL_del_a*del_a
 * with:
 *   CL0 = -1.36422*Cmu + 0.2979
 *   CL_alpha = 2*PI*(1+0.151*(Cmu)^0.5 + 0.219*Cmu)*PI/180
 *   CL_del_e = 0.0141
 *   CL_Cmu = 3.168
 *   CL_del_f = 0.0459*Cmu + 0.0135
 *   CL_del_a = 0
 *
 * Notes:
 * - alpha_deg, del_*_deg must be in degrees (the CL_alpha factor includes a deg->rad scaling in formula).
 * - throttle only needed to compute Cmu via n.
 */
static inline double compute_CL_from_model(double alpha_deg,
                                           double del_e_deg,
                                           double del_f_deg,
                                           double del_a_deg,
                                           double throttle,
                                           double V_inf) 
{
    // compute rotor speed n (rps)
    double n = compute_n_from_throttle(throttle);

    // compute Cmu
    double Cmu = compute_Cmu(V_inf, n, 0.120);

    // coefficients
    double CL0 = -1.36422 * Cmu + 0.2979;
    double CL_alpha = 2.0 * M_PI * (1.0 + 0.151 * sqrt(Cmu) + 0.219 * Cmu) * (M_PI / 180.0);
    double CL_del_e = 0.0141;
    double CL_Cmu = 3.168;
    double CL_del_f = 0.0459 * Cmu + 0.0135;
    double CL_del_a = 0.0;

    double CL = CL0
                + CL_alpha * alpha_deg
                + CL_del_e * del_e_deg
                + CL_del_f * del_f_deg
                + CL_Cmu * Cmu
                + CL_del_a * del_a_deg;

    return CL;
}


// ----------------- compute_CD_from_model -----------------
// CL_linear is the dimensionless lift coefficient (from analytic model).
// alpha_deg, del_*_deg are in DEGREES.
// throttle in [0,1], V_inf in m/s.
static inline double compute_CD_from_model(double CL_linear,
                                           double alpha_deg,
                                           double del_e_deg,
                                           double del_f_deg,
                                           double del_a_deg,
                                           double throttle,
                                           double V_inf)
{
    // 1) compute Cmu (re-uses compute_n_from_throttle / compute_Cmu)
    double n = compute_n_from_throttle(throttle);
    double Cmu = compute_Cmu(V_inf, n, 0.120);

    // 2) CD0 base pieces
    // CD0c = -0.51403*Cmu + 0.107575
    double CD0c = -0.51403 * Cmu + 0.107575;

    // CD0_del_f term ( "0.01285*del_f + 0.0021")
    // Interpret: component = 0.01285 * del_f_deg + 0.0021
    double CD0_del_f_term = 0.01285 * del_f_deg + 0.0021;

    // CD0 = CD0c + CD0_alpha*alpha + CD0_del_e*del_e + CD0_del_f_term + CD0_del_a*del_a + CD0_Cmu*Cmu
    // CD0_alpha = 0, CD0_del_e = 0.0038, CD0_del_a = 0.0015, CD0_Cmu = 1.322175
    double CD0_alpha = 0.0;
    double CD0_del_e = 0.0038;
    double CD0_del_a = 0.0015;
    double CD0_Cmu   = 1.322175;

    double CD0 = CD0c
                 + CD0_alpha * alpha_deg
                 + CD0_del_e * del_e_deg
                 + CD0_del_f_term
                 + CD0_del_a * del_a_deg
                 + CD0_Cmu * Cmu;

    // 3) k interpolation vs flap angle del_f_deg
    // k = 0.046 at del_f = 0 deg
    // k = 0.056 at del_f = 40 deg
    double del_f_clamped = del_f_deg;
    if (del_f_clamped < 0.0) del_f_clamped = 0.0;
    if (del_f_clamped > 40.0) del_f_clamped = 40.0;
    const double k0 = 0.046;
    const double k40 = 0.056;
    double k = k0 + (del_f_clamped / 40.0) * (k40 - k0);

    // 4) final CD
    double CD = CD0 + k * (CL_linear * CL_linear);

    return CD;
}



using namespace SITL;

Plane::Plane(const char *frame_str) :
    Aircraft(frame_str)
{

    const char *colon = strchr(frame_str, ':');
    size_t slen = strlen(frame_str);
    // The last 5 letters are ".json"
    if (colon != nullptr && slen > 5 && strcmp(&frame_str[slen-5], ".json") == 0) {
        load_coeffs(colon+1);
    } else {
        coefficient = default_coefficients;
    }

    //mass = 2.0f;
    mass = 65.0f;
    

    //inertias = 3,4,5;

    /*
       scaling from motor power to Newtons. Allows the plane to hold
       vertically against gravity when the motor is at hover_throttle
    */
    float hover_const = 1.3;
    thrust_scale = (mass * GRAVITY_MSS) / hover_const;
    frame_height = 0.1f;

    ground_behavior = GROUND_BEHAVIOR_FWD_ONLY;
    lock_step_scheduled = true;

    if (strstr(frame_str, "-heavy")) {
        mass = 8;
    }
    if (strstr(frame_str, "-jet")) {
        // a 22kg "jet", level top speed is 102m/s
        mass = 22;
        thrust_scale = (mass * GRAVITY_MSS) / hover_throttle;
    }
    if (strstr(frame_str, "-revthrust")) {
        reverse_thrust = true;
    }
    if (strstr(frame_str, "-elevon")) {
        elevons = true;
    } else if (strstr(frame_str, "-vtail")) {
        vtail = true;
    } else if (strstr(frame_str, "-dspoilers")) {
        dspoilers = true;
    } else if (strstr(frame_str, "-redundant")) {
        redundant = true;
    }
    if (strstr(frame_str, "-elevrev")) {
        reverse_elevator_rudder = true;
    }
    if (strstr(frame_str, "-catapult")) {
        have_launcher = true;
        launch_accel = 15;
        launch_time = 2;
    }
    if (strstr(frame_str, "-bungee")) {
        have_launcher = true;
        launch_accel = 7;
        launch_time = 4;
    }
    if (strstr(frame_str, "-throw")) {
        have_launcher = true;
        launch_accel = 25;
        launch_time = 0.4;
    }
    if (strstr(frame_str, "-tailsitter")) {
        tailsitter = true;
        ground_behavior = GROUND_BEHAVIOR_TAILSITTER;
        thrust_scale *= 1.5;
    }
    if (strstr(frame_str, "-steering")) {
        have_steering = true;
    }

#if AP_FILESYSTEM_FILE_READING_ENABLED
    if (strstr(frame_str, "-3d")) {
        aerobatic = true;
        thrust_scale *= 1.5;
        // setup parameters for plane-3d
        AP_Param::load_defaults_file("@ROMFS/models/plane.parm", false);
        AP_Param::load_defaults_file("@ROMFS/models/plane-3d.parm", false);
    }
#endif

    if (strstr(frame_str, "-ice")) {
        ice_engine = true;
    }

    if (strstr(frame_str, "-soaring")) {
        mass = 2.0;
        coefficient.c_drag_p = 0.05;
    }
}

void Plane::load_coeffs(const char *model_json)
{
    char *fname = nullptr;
    struct stat st;
    if (AP::FS().stat(model_json, &st) == 0) {
        fname = strdup(model_json);
    } else {
        IGNORE_RETURN(asprintf(&fname, "@ROMFS/models/%s", model_json));
        if (AP::FS().stat(model_json, &st) != 0) {
            AP_HAL::panic("%s failed to load", model_json);
        }
    }
    if (fname == nullptr) {
        AP_HAL::panic("%s failed to load", model_json);
    }
    AP_JSON::value *obj = AP_JSON::load_json(model_json);
    if (obj == nullptr) {
        AP_HAL::panic("%s failed to load", model_json);
    }

    enum class VarType {
        FLOAT,
        VECTOR3F,
    };

    struct json_search {
        const char *label;
        void *ptr;
        VarType t;
    };
    
    json_search vars[] = {
#define COFF_FLOAT(s) { #s, &coefficient.s, VarType::FLOAT }
        COFF_FLOAT(s),
        COFF_FLOAT(b),
        COFF_FLOAT(c),
        COFF_FLOAT(c_lift_0),
        COFF_FLOAT(c_lift_deltae),
        COFF_FLOAT(c_lift_a),
        COFF_FLOAT(c_lift_q),
        COFF_FLOAT(mcoeff),
        COFF_FLOAT(oswald),
        COFF_FLOAT(alpha_stall),
        COFF_FLOAT(c_drag_q),
        COFF_FLOAT(c_drag_deltae),
        COFF_FLOAT(c_drag_p),
        COFF_FLOAT(c_y_0),
        COFF_FLOAT(c_y_b),
        COFF_FLOAT(c_y_p),
        COFF_FLOAT(c_y_r),
        COFF_FLOAT(c_y_deltaa),
        COFF_FLOAT(c_y_deltar),
        COFF_FLOAT(c_l_0),
        COFF_FLOAT(c_l_p),
        COFF_FLOAT(c_l_b),
        COFF_FLOAT(c_l_r),
        COFF_FLOAT(c_l_deltaa),
        COFF_FLOAT(c_l_deltar),
        COFF_FLOAT(c_m_0),
        COFF_FLOAT(c_m_a),
        COFF_FLOAT(c_m_q),
        COFF_FLOAT(c_m_deltae),
        COFF_FLOAT(c_n_0),
        COFF_FLOAT(c_n_b),
        COFF_FLOAT(c_n_p),
        COFF_FLOAT(c_n_r),
        COFF_FLOAT(c_n_deltaa),
        COFF_FLOAT(c_n_deltar),
        COFF_FLOAT(deltaa_max),
        COFF_FLOAT(deltae_max),
        COFF_FLOAT(deltar_max),
        { "CGOffset", &coefficient.CGOffset, VarType::VECTOR3F },
    };

    for (uint8_t i=0; i<ARRAY_SIZE(vars); i++) {
        auto v = obj->get(vars[i].label);
        if (v.is<AP_JSON::null>()) {
            // use default value
            continue;
        }
        if (vars[i].t == VarType::FLOAT) {
            parse_float(v, vars[i].label, *((float *)vars[i].ptr));

        } else if (vars[i].t == VarType::VECTOR3F) {
            parse_vector3(v, vars[i].label, *(Vector3f *)vars[i].ptr);

        }
    }

    delete obj;

    ::printf("Loaded plane aero coefficients from %s\n", model_json);
}

void Plane::parse_float(AP_JSON::value val, const char* label, float &param) {
    if (!val.is<double>()) {
        AP_HAL::panic("Bad json type for %s: %s", label, val.to_str().c_str());
    }
    param = val.get<double>();
}

void Plane::parse_vector3(AP_JSON::value val, const char* label, Vector3f &param) {
    if (!val.is<AP_JSON::value::array>() || !val.contains(2) || val.contains(3)) {
        AP_HAL::panic("Bad json type for %s: %s", label, val.to_str().c_str());
    }
    for (uint8_t j=0; j<3; j++) {
        parse_float(val.get(j), label, param[j]);
    }
}

/*
  the following functions are from last_letter
  https://github.com/Georacer/last_letter/blob/master/last_letter/src/aerodynamicsLib.cpp
  many thanks to Georacer!
 */
// float Plane::liftCoeff(float alpha) const
// {
//     // const float alpha0 = coefficient.alpha_stall;
//     // const float M = coefficient.mcoeff;
//     // const float c_lift_0 = coefficient.c_lift_0;
//     // const float c_lift_a0 = coefficient.c_lift_a;

//     // // clamp the value of alpha to avoid exp(90) in calculation of sigmoid
//     // const float max_alpha_delta = 0.8f;
//     // if (alpha-alpha0 > max_alpha_delta) {
//     //     alpha = alpha0 + max_alpha_delta;
//     // } else if (alpha0-alpha > max_alpha_delta) {
//     //     alpha = alpha0 - max_alpha_delta;
//     // }
// 	// double sigmoid = ( 1+exp(-M*(alpha-alpha0))+exp(M*(alpha+alpha0)) ) / (1+exp(-M*(alpha-alpha0))) / (1+exp(M*(alpha+alpha0)));
// 	// double linear = (1.0-sigmoid) * (c_lift_0 + c_lift_a0*alpha); //Lift at small AoA
// 	// double flatPlate = sigmoid*(2*copysign(1,alpha)*pow(sin(alpha),2)*cos(alpha)); //Lift beyond stall

// 	// float result  = linear+flatPlate;
// 	// return result;
//     const float c_lift_0 = coefficient.c_lift_0;
//     const float c_lift_a = coefficient.c_lift_a;

//     // old linear model (per-radian since alpha is radians here)
//     float CL_linear_old = c_lift_0 + c_lift_a * alpha;

//     // call the new blender
//     return liftCoeff(alpha, CL_linear_old);
    
// }

/* New: blending lift function
   - alpha in radians
   - CL_linear is the dimensionless linear CL to be used pre-stall
   Returns blended CL (linear pre-stall smoothly blended to flat-plate post-stall).
*/
// float Plane::liftCoeff(float alpha, float CL_linear) const
// {
//     const float alpha0 = coefficient.alpha_stall;
//     const float M = coefficient.mcoeff;

//     // clamp alpha to avoid huge exp()
//     const float max_alpha_delta = 0.8f;
//     float a = alpha;
//     if (a - alpha0 > max_alpha_delta) a = alpha0 + max_alpha_delta;
//     else if (alpha0 - a > max_alpha_delta) a = alpha0 - max_alpha_delta;

//     double sigmoid = ( 1.0 + exp(-M*(a - alpha0)) + exp(M*(a + alpha0)) )
//                      / ((1.0 + exp(-M*(a - alpha0))) * (1.0 + exp(M*(a + alpha0))));

//     double linear = (1.0 - sigmoid) * (double)CL_linear;
//     double flatPlate = sigmoid * (2.0 * copysign(1.0, alpha) * pow(sin(alpha), 2) * cos(alpha));

//     return (float)(linear + flatPlate);
// }




// float Plane::dragCoeff(float alpha) const
// {
//     const float b = coefficient.b;
//     const float s = coefficient.s;
//     const float c_drag_p = coefficient.c_drag_p;
//     const float c_lift_0 = coefficient.c_lift_0;
//     const float c_lift_a0 = coefficient.c_lift_a;
//     const float oswald = coefficient.oswald;
    
// 	double AR = pow(b,2)/s;
// 	double c_drag_a = c_drag_p + pow(c_lift_0+c_lift_a0*alpha,2)/(M_PI*oswald*AR);

// 	return c_drag_a;
// }



// SITL::Wrench Plane::getForcesAndMoments(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const struct sitl_input &input, bool fm)
// {   float alpha = angle_of_attack;
//     float radtodeg = 180/(3.14);
//     printf("alpha = %.3f",alpha*radtodeg);
// 	//calculate aerodynamic torque
//     float effective_airspeed = airspeed;

//     if (tailsitter || aerobatic) {
//         /*
//           tailsitters get airspeed from prop-wash
//          */
//         effective_airspeed += inputThrust * 20;

//         // reduce effective angle of attack as thrust increases
//         alpha *= constrain_float(1 - inputThrust, 0, 1);
//     }

//     const float c_drag_q = coefficient.c_drag_q;
//     const float c_lift_q = coefficient.c_lift_q;
//     const float s = coefficient.s;
//     const float c = coefficient.c;
//     const float b = coefficient.b;
//     const float c_drag_deltae = coefficient.c_drag_deltae;
//     const float c_lift_deltae = coefficient.c_lift_deltae;
//     const float c_y_0 = coefficient.c_y_0;
//     const float c_y_b = coefficient.c_y_b;
//     const float c_y_p = coefficient.c_y_p;
//     const float c_y_r = coefficient.c_y_r;
//     const float c_y_deltaa = coefficient.c_y_deltaa;
//     const float c_y_deltar = coefficient.c_y_deltar;
//     const float c_drag_0 = coefficient.c_drag_0;
//     const float c_lift_0 = coefficient.c_lift_0;
//     const float c_l_0 = coefficient.c_l_0;
//     const float c_l_b = coefficient.c_l_b;
//     const float c_l_p = coefficient.c_l_p;
//     const float c_l_r = coefficient.c_l_r;
//     const float c_l_deltaa = coefficient.c_l_deltaa;
//     const float c_l_deltar = coefficient.c_l_deltar;
//     const float c_m_0 = coefficient.c_m_0;
//     // const float c_m_a = coefficient.c_m_a;
//     float c_m_a;
//     if (alpha<0){
//         c_m_a = 0.06*radtodeg;
//     }
//     else if (alpha>0 && alpha<(5/radtodeg))
//     {
//         c_m_a = 0.15*radtodeg;
//     }
//     else{
//         c_m_a = -0.1*radtodeg;
//     }
//     printf("Cm_alpha=%0.3f",c_m_a);
//     const float c_m_q = coefficient.c_m_q;
//     const float c_m_deltae = coefficient.c_m_deltae;
//     const float c_n_0 = coefficient.c_n_0;
//     const float c_n_b = coefficient.c_n_b;
//     const float c_n_p = coefficient.c_n_p;
//     const float c_n_r = coefficient.c_n_r;
//     const float c_n_deltaa = coefficient.c_n_deltaa;
//     const float c_n_deltar = coefficient.c_n_deltar;
//     const float Lf = 0.858f;    // CG to nose gear (+x) 
//     const float Lb = 0.122f;    // CG to main gear (aft)

//     float rho = air_density;

//     float throttle;
//     if (reverse_thrust) {
//         throttle = filtered_servo_angle(input, 2);
//     } else {
//         throttle = filtered_servo_range(input, 2);
//     }
    
//     float thrust     = throttle;
//     thrust *= thrust_scale;
// 	// //request lift and drag alpha-coefficients from the corresponding functions
// 	// double c_lift_a = liftCoeff(alpha);
// 	// double c_drag_a = dragCoeff(alpha);

//     // Adding code here to compute lift coefficient for EQUINOX

//     double alpha_deg  = alpha * (180.0 / M_PI);
//     double del_e_deg  = inputElevator * (180.0 / M_PI);
//     double del_a_deg  = inputAileron * (180.0 / M_PI);
//     double del_f_deg  = 0.0; // replace with actual flap angle in degrees if available

//         double throttle_for_model = clamp_double((double)throttle, 0.0, 1.0);
//     // Debug print the control inputs
//     printf("Control Inputs: alpha=%.1f° del_e=%.1f° del_a=%.1f° del_f=%.1f° throttle=%.2f airspeed=%.1f\n",
//         alpha_deg,
//         del_e_deg, 
//         del_a_deg,
//         del_f_deg,
//         throttle_for_model,
//         airspeed);

//     // compute the analytic linear CL (degrees inputs)
//     double CL_linear = compute_CL_from_model(alpha_deg, del_e_deg, del_f_deg, del_a_deg, throttle_for_model, (double)airspeed);

//     // get blended (smooth) CL from the new liftCoeff blender (alpha in radians)
//     double c_lift_a = liftCoeff(alpha, (float)CL_linear);

//     // compute drag consistently using CL_linear for induced portion

//     double c_drag_a = compute_CD_from_model(CL_linear, alpha_deg, del_e_deg, del_f_deg, del_a_deg, throttle_for_model, (double)airspeed);




// 	// //convert coefficients to the body frame
//     double c_x_0 = -c_drag_0*cos(alpha)+c_lift_0*sin(alpha);
// 	double c_x_a = -c_drag_a*cos(alpha)+c_lift_a*sin(alpha);
// 	double c_x_q = -c_drag_q*cos(alpha)+c_lift_q*sin(alpha);
//     double c_z_0 = -c_drag_0*sin(alpha)-c_lift_0*cos(alpha);
// 	double c_z_a = -c_drag_a*sin(alpha)-c_lift_a*cos(alpha);
// 	double c_z_q = -c_drag_q*sin(alpha)-c_lift_q*cos(alpha);






//     // double aerocfs[6];
    
//     // float deg_inputa = inputAileron*radtodeg;
//     // float deg_inpute = inputElevator*radtodeg;
//     // float deg_inputr = inputRudder*radtodeg;
//     // float deg_beta = beta*radtodeg;
//     // float deg_alpha = alpha*radtodeg;

//     //printf("thrust = %.3f, throttle = %.3f,inputAileron=%.3f,inputElevator=%.3f",thrust, throttle*100,deg_inputa,deg_inpute);
//     // fm = 0;
//     // aero_interpolate(throttle*100, deg_inputa, deg_inpute, deg_inputr, deg_beta, deg_alpha, aerocfs,fm);

//     // float CL_direct = aerocfs[0];
// 	// float CD_direct = aerocfs[1];
// 	// float CY_direct = aerocfs[2];
// 	// float Cl_direct = aerocfs[3];
// 	// float Cm_direct = aerocfs[4];
// 	// float Cn_direct = aerocfs[5];
//     const float phi = AP::ahrs().get_pitch();
//     // float CX_direct = CL_direct*sin(alpha)-CD_direct*cos(alpha)*cos(beta)-CY_direct*cos(alpha)*sin(beta);
//     // float CZ_direct = -CL_direct*cos(alpha)-CD_direct*sin(alpha)*cos(beta)-CY_direct*sin(alpha)*sin(beta);    
//     // printf("CL=%.3f,CD=%.3f,CY=%.3f,Cl=%.3f,Cm=%.3f,Cn=%.3f",CL_direct,CD_direct,CY_direct,Cl_direct,Cm_direct,Cn_direct);
//     // printf("throttle = %.3f, aileron = %.3f, elevator = %.3f, rudder = %.3f, beta = %.3f, alpha = %.3f",throttle*100, deg_inputa, deg_inpute, deg_inputr, deg_beta, deg_alpha);
// 	//read angular rates
// 	double p = gyro.x;
// 	double q = gyro.y;
// 	double r = gyro.z;
    
//     float thrust_offset = 0.091;
// 	//calculate aerodynamic force
// 	double qbar = 1.0/2.0*rho*pow(airspeed,2)*s; //Calculate dynamic pressure
//     float ax = 0.0f, ay = 0.0f, az = 0.0f;   // body forces
//     float la = 0.0f, ma = 0.0f, na = 0.0f;   // body moments

//     //double Nf,Nb;
//     if (is_zero(airspeed))
// 	{
// 		ax = 0;
// 		ay = 0;
// 		az = 0;
//         la = 0;
// 		ma = 0;
// 		na = 0;
// 	}
//     else{
//         float fx_aero_b = qbar*(c_x_0 + c_x_a + c_x_q*c*q/(2*airspeed) - c_drag_deltae*cos(alpha)*fabs(inputElevator) + c_lift_deltae*sin(alpha)*inputElevator);
// 		// split c_x_deltae to include "abs" term
// 		float fy_aero_b = qbar*(c_y_0 + c_y_b*beta + c_y_p*b*p/(2*airspeed) + c_y_r*b*r/(2*airspeed) + c_y_deltaa*inputAileron + c_y_deltar*inputRudder);
// 		float fz_aero_b = qbar*(c_z_0 + c_z_a + c_z_q*c*q/(2*airspeed) - c_drag_deltae*sin(alpha)*fabs(inputElevator) - c_lift_deltae*cos(alpha)*inputElevator);
//         // float fx_aero_b = qbar*CX_direct;
//         // float fy_aero_b = qbar*CY_direct;
//         // float fz_aero_b = qbar*CZ_direct;
//         float fz_aero_e = sin(phi)*fx_aero_b + cos(phi)*fz_aero_b;
//         float fz_thrust_e = -thrust*sin(phi);
//         float l_aero = qbar*b*(c_l_0 + c_l_b*beta + c_l_p*b*p/(2*effective_airspeed) + c_l_r*b*r/(2*effective_airspeed) + c_l_deltaa*inputAileron + c_l_deltar*inputRudder);
// 		float m_aero = qbar*c*(c_m_0 + c_m_a*alpha + c_m_q*c*q/(2*effective_airspeed) + c_m_deltae*inputElevator);
// 		float n_aero = qbar*b*(c_n_0 + c_n_b*beta + c_n_p*b*p/(2*effective_airspeed) + c_n_r*b*r/(2*effective_airspeed) + c_n_deltaa*inputAileron + c_n_deltar*inputRudder);
//         // float l_aero = qbar*b*Cl_direct;
//         // float m_aero = qbar*c*Cm_direct;
//         // float n_aero = qbar*b*Cn_direct;
//         float m_thrust = thrust*thrust_offset;
    

//         // -------------------- Liftoff/landing state machine --------------------
//        // Tunables
//         const float FORCE_MARGIN         = 0.02f;   // 2% margin above/below mg
//         const int   FORCE_HOLD_FRAMES    = 5;       // consecutive frames to confirm force trigger
//         const float ALT_HYSTERESIS       = 0.010f;  // 1 cm AGL to latch airborne (wheel clearance)
//         const float PEN_ENGAGE           = 0.005f;  // 5 mm engage window (up-positive)
//         const float MIN_LOAD_N           = 0.1f;    // treat tiny normals as zero (reduced for debug)

//         // altitude-stability detector (frames/tolerance)
//         const int   ALT_STABLE_FRAMES    = 5000;       // consecutive frames altitude must be steady
//         const float ALT_STABLE_TOL       = 0.002f;   // meters tolerance (2 mm)
//         const float VZ_TOUCH_THRESH      = 0.6f;    // m/s; vertical speed threshold for landing detection

//         // 1) Geometry
//         const Vector3f rNose_b(+Lf, 0.0f, frame_height);
//         const Vector3f rMain_b(-Lb, 0.0f, frame_height);
//         const Vector3f rNose_e = dcm.transposed() * rNose_b;   // body -> earth (NED)
//         const Vector3f rMain_e = dcm.transposed() * rMain_b;

//         // 1a) relative position (NED: z down-positive)
//         Vector3f locned;
//         bool alt_gotten = AP::ahrs().get_relative_position_NED_origin_float(locned);
//         // origin_hagl_up: up-positive aircraft origin height (if alt_gotten false -> 0)
//         const float origin_hagl_up = alt_gotten ? -locned.z : 0.0f;

//         // wheel world Z (down-positive)
//         const float wheelNose_world_z = locned.z + rNose_e.z;   // down-positive
//         const float wheelMain_world_z = locned.z + rMain_e.z;   // down-positive

//         // latched touchdown wheel world Z (down-positive)
//         static float ground_wheel_z_touch = 0.0f;
//         static bool  ground_wheel_z_latched = false;

//         // compute wheel clearance (up-positive) relative to latched wheel world z if latched,
//         // otherwise relative to immediate world z=0
//         float wheelNose_AGL_up, wheelMain_AGL_up;
//         if (ground_wheel_z_latched) {
//             wheelNose_AGL_up = ground_wheel_z_touch - wheelNose_world_z; // up-positive
//             wheelMain_AGL_up = ground_wheel_z_touch - wheelMain_world_z;
//         } else {
//             wheelNose_AGL_up = -wheelNose_world_z;
//             wheelMain_AGL_up = -wheelMain_world_z;
//         }
//         wheelNose_AGL_up = std::max(0.0f, wheelNose_AGL_up);
//         wheelMain_AGL_up = std::max(0.0f, wheelMain_AGL_up);
//         const float minWheel_AGL_up = std::min(wheelNose_AGL_up, wheelMain_AGL_up);
//         const float AGL_up = minWheel_AGL_up;

//         // 2) S1/S2 provisional split (unchanged)
//         const float S1 = mass*GRAVITY_MSS + fz_aero_e + fz_thrust_e; // down-positive
//         float cosph_raw = std::cos(phi);
//         float cosph = (std::fabs(cosph_raw) < 0.02f) ? (copysign(0.02f, cosph_raw)) : cosph_raw;
//         const float S2 = -(m_aero + m_thrust) / cosph;
//         float Nf_prov = ( S2 + Lb*S1 ) / (Lf + Lb);
//         float Nb_prov = ( Lf*S1 - S2 ) / (Lf + Lb);
//         Nf_prov = std::max(0.0f, Nf_prov);
//         Nb_prov = std::max(0.0f, Nb_prov);

//         // 3) Force trigger (unchanged)
//         const float lift_up   = -fz_aero_e;    // up-positive
//         const float thrust_up = -fz_thrust_e;  // up-positive
//         const float Fz_up     = lift_up + thrust_up; // up-positive net upward force

//         static float Fz_up_f = 0.0f;
//         const float EMA = 0.3f;
//         Fz_up_f = (1.0f - EMA) * Fz_up_f + EMA * Fz_up;

//         const bool force_exceeds = (Fz_up_f > mass*GRAVITY_MSS * (1.0f + FORCE_MARGIN));
//         const bool force_below   = (Fz_up_f < mass*GRAVITY_MSS * (1.0f - FORCE_MARGIN));

//         // 4) altitude-stability bookkeeping (shared across states)
//         static float alt_last = 0.0f;
//         static int   alt_stable_cnt = 0;
//         static bool has_taken_off = false;   
//         // update altitude-stability counter if we have valid HAGL
//         if (alt_gotten) {
//             const float dalt = fabsf(origin_hagl_up - alt_last);
//             if (dalt <= ALT_STABLE_TOL) {
//                 alt_stable_cnt++;
//             } else {
//                 alt_stable_cnt = 0;
//             }
//             alt_last = origin_hagl_up;
//         } else {
//             // no alt available, reset stability counter
//             alt_stable_cnt = 0;
//         }

//         // 5) State machine
//         enum class AirState : uint8_t { GROUND, FORCE_OK, AIRBORNE };
//         static AirState air_state = AirState::GROUND;
//         static int force_cnt = 0;

//         float Nf = 0.0f, Nb = 0.0f;

//         switch (air_state) {
//         case AirState::GROUND:
//         {
//             // If not latched yet (start-on-ground), latch now
//             if (!ground_wheel_z_latched && alt_gotten) {
//                 // choose the wheel with the smallest world z (closest to ground in down-positive)
//                 ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
//                 ground_wheel_z_latched = true;
//                 std::printf("DBG: initial-ground latch wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//             }

//             Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//             Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;

//             // Advance to FORCE_OK only if net upward force persistently > mg
//             force_cnt = force_exceeds ? (force_cnt + 1) : 0;
//             if (force_cnt >= FORCE_HOLD_FRAMES) {
//                 air_state = AirState::FORCE_OK;
//                 force_cnt = 0;
//             }
//             break;
//         }

//         case AirState::FORCE_OK:
//         {
//             // still treat as ground until clearance observed
//             Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//             Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//             fm = 1;
//             // commit to AIRBORNE only if we see sustained wheel clearance
//             if (AGL_up >= ALT_HYSTERESIS) {
//                 air_state = AirState::AIRBORNE;
//                 has_taken_off = true; 
//                 // clear latched touchdown reference so next landing re-latches
//                 ground_wheel_z_latched = false;
//                 Nf = Nb = 0.0f;
//                 break;
//             }

//             // if upward force collapses before clearance, return to GROUND
//             if (!has_taken_off && force_below) {
//                 air_state = AirState::GROUND;
//                 force_cnt = 0;
//                 // keep ground_wheel_z_latched as-is (do not clear)
//             }
//             break;
//         }

//         case AirState::AIRBORNE:
//         {
//             // normals zero while airborne
//             Nf = Nb = 0.0f;
//             fm = 1;

//             // Wheel-penetration touchdown detection (fast)
//             const bool nose_touch = (wheelNose_AGL_up <= PEN_ENGAGE);
//             const bool main_touch = (wheelMain_AGL_up <= PEN_ENGAGE);
//             if (nose_touch || main_touch) {
//                 air_state = AirState::GROUND;
//                 printf("CONDITION1\n");
//                 //exit(0);
//                 Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//                 Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//                 // latch world wheel z at touchdown (prefer main wheel)
//                 ground_wheel_z_touch = main_touch ? wheelMain_world_z : wheelNose_world_z;
//                 ground_wheel_z_latched = true;
//                 has_taken_off = false; // landed -> allow next ground latch on next takeoff
//                 std::printf("DBG: wheel-touch touchdown latched wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//                 force_cnt = 0;
//                 break;
//             }

//             // Robust fallback: force_below + low vertical speed + stable altitude for several frames
//             const float vz_down = velocity_ef.z; // NED down-positive
//             static int touchdown_force_cnt = 0;
//             if (Fz_up_f<10 && fabsf(vz_down) < VZ_TOUCH_THRESH && alt_stable_cnt >= ALT_STABLE_FRAMES) {
//                 touchdown_force_cnt++;
//             } else {
//                 touchdown_force_cnt = 0;
//             }
//             if (touchdown_force_cnt >= FORCE_HOLD_FRAMES) {
//                 printf("CONDITION2\n");
//                 //exit(0);
//                 air_state = AirState::GROUND;
//                 fm=0;
//                 Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//                 Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//                 ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
//                 ground_wheel_z_latched = true;
//                 has_taken_off = false;
//                 touchdown_force_cnt = 0;
//                 std::printf("DBG: force+stable-alt touchdown latched wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//                 break;
//             }
//             break;
//         } // end AIRBORNE
//         } // end switch

//         // debug print
//         // std::printf("DBG: alt_gotten=%d origin_hagl_up=%.3f ground_latched=%d wheel_z_touch=%.3f "
//         //             "wheelNose_z=%.3f wheelMain_z=%.3f wheelNose_AGL=%.3f wheelMain_AGL=%.3f AGL=%.3f "
//         //             "Fz_up_f=%.2f force_exceeds=%d force_below=%d alt_stable_cnt=%d state=%d\n",
//         //             int(alt_gotten),
//         //             origin_hagl_up,
//         //             int(ground_wheel_z_latched),
//         //             ground_wheel_z_touch,
//         //             wheelNose_world_z,
//         //             wheelMain_world_z,
//         //             wheelNose_AGL_up,
//         //             wheelMain_AGL_up,
//         //             AGL_up,
//         //             Fz_up_f,
//         //             int(force_exceeds),
//         //             int(force_below),
//         //             alt_stable_cnt,
//         //             int(air_state)
//         // );

//         // 5) Finish as before: project normals back to BODY and add to aero-only forces
//         std::printf("state=%d fm =%d AGL=%.3f Fz_up=%.1f Nf=%.1f Nb=%.1f\n",
//                     int(air_state),int(fm), AGL_up, Fz_up, Nf, Nb);
//         //std::printf("flightmode = %d,CL=%.5f,CD=%.5f,CM=%.5f,CY=%.5f\n",fm,CL_direct,CD_direct,Cm_direct,CY_direct);

//         const float fx_norm_b =  (Nf + Nb) * std::sin(phi);
//         const float fz_norm_b = -(Nf + Nb) * std::cos(phi);

//         ax = fx_aero_b + fx_norm_b;   // you return aero + ground only
//         ay = fy_aero_b;
//         ay = 0;
//         az = fz_aero_b + fz_norm_b;

//         la = l_aero;
//         ma = m_aero + Nb*Lb*std::cos(phi) - Nf*Lf*std::cos(phi) + m_thrust;
//         na = n_aero;
//         la = 0;
//         na = 0;


//     }
//     return Wrench{ Vector3f(ax, ay, az), Vector3f(la, ma, na) };

// }

// SITL::Wrench Plane::getForcesAndMoments(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const struct sitl_input &input, bool fm)
// {
//     // --- SELECTOR FLAG ---
//     // 0 -> use interpolator (preferred). non-zero -> use analytical model.
//     // If interpolator fails for a particular point, we automatically fall back to analytical.
//     const int use_interpolator_flag = 0;
//     const float phi = AP::ahrs().get_pitch();

//     // --- existing setup & local copies from your original function ---
//     float alpha = angle_of_attack;
//     const double radtodeg = 180.0 / M_PI;
//     printf("alpha = %.3f\n", alpha * radtodeg);

//     float effective_airspeed = airspeed;

//     if (tailsitter || aerobatic) {
//         effective_airspeed += inputThrust * 20;
//         alpha *= constrain_float(1 - inputThrust, 0, 1);
//     }

//     const float c_drag_q = coefficient.c_drag_q;
//     const float c_lift_q = coefficient.c_lift_q;
//     const float s = coefficient.s;
//     const float c = coefficient.c;
//     const float b = coefficient.b;
//     const float c_drag_deltae = coefficient.c_drag_deltae;
//     const float c_lift_deltae = coefficient.c_lift_deltae;
//     const float c_y_0 = coefficient.c_y_0;
//     const float c_y_b = coefficient.c_y_b;
//     const float c_y_p = coefficient.c_y_p;
//     const float c_y_r = coefficient.c_y_r;
//     const float c_y_deltaa = coefficient.c_y_deltaa;
//     const float c_y_deltar = coefficient.c_y_deltar;
//     const float c_drag_0 = coefficient.c_drag_0;
//     const float c_lift_0 = coefficient.c_lift_0;
//     const float c_l_0 = coefficient.c_l_0;
//     const float c_l_b = coefficient.c_l_b;
//     const float c_l_p = coefficient.c_l_p;
//     const float c_l_r = coefficient.c_l_r;
//     const float c_l_deltaa = coefficient.c_l_deltaa;
//     const float c_l_deltar = coefficient.c_l_deltar;
//     const float c_m_0 = coefficient.c_m_0;
//     const float c_m_q = coefficient.c_m_q;
//     const float c_m_deltae = coefficient.c_m_deltae;
//     const float c_n_0 = coefficient.c_n_0;
//     const float c_n_b = coefficient.c_n_b;
//     const float c_n_p = coefficient.c_n_p;
//     const float c_n_r = coefficient.c_n_r;
//     const float c_n_deltaa = coefficient.c_n_deltaa;
//     const float c_n_deltar = coefficient.c_n_deltar;
//     const float Lf = 0.858f;    // CG to nose gear (+x)
//     const float Lb = 0.122f;    // CG to main gear (aft)

//     float rho = air_density;

//     float throttle;
//     if (reverse_thrust) {
//         throttle = filtered_servo_angle(input, 2);
//     } else {
//         throttle = filtered_servo_range(input, 2);
//     }

//     float thrust = throttle;
//     thrust *= thrust_scale;

//     // Convert angles and controls to degrees for model / interpolator as needed
//     double alpha_deg  = alpha * (180.0 / M_PI);
//     double del_e_deg  = inputElevator * (180.0 / M_PI);
//     double del_a_deg  = inputAileron * (180.0 / M_PI);
//     double del_r_deg  = inputRudder * (180.0 / M_PI);
//     double del_f_deg  = 40.0; // replace with actual flap angle in degrees if available
//     double beta_deg   = beta * (180.0 / M_PI);   // assume `beta` is a member variable in radians

//     double throttle_for_model = clamp_double((double)throttle, 0.0, 1.0);

    
//     // Placeholders for aerodynamic coefficients (body-frame conversion later uses c_lift_a & c_drag_a)
//     double c_lift_a = 0.0;
//     double c_drag_a = 0.0;
//     double CM_pred  = 0.0;
//     bool   used_interpolator_success = false;

//     // --- Attempt interpolation (if selected) ---
//     if (use_interpolator_flag == 0) {
//         // Our interpolator expects:
//         //   throttle in percent (0..100), aileron,elevator,rudder,alpha,beta in degrees.
//             // ---------------------------------------------------------------------
//         // Interpolation attempt (with hardcoded flaps & aileron = 40 deg)
//         // ---------------------------------------------------------------------
//         double interp_aerocoeffs[ (size_t)AERO_M > 0 ? AERO_M : 3 ];
//         // throttle_for_model is 0..1; interpolator expects 0..100
//         double throttle_pct = throttle_for_model * 100.0;

//         // Use actual elevator, rudder, alpha, beta (in degrees) but FORCE aileron & flaps to 40°
//         double interp_del_a_deg = 40.0;      // hardcoded for interpolation only
//         double interp_del_f_deg = 40.0;   
//         del_f_deg = interp_del_f_deg;
//         del_a_deg = interp_del_a_deg;
//         // hardcoded flap and aileron angle (if used by model)
//         // note: we pass interp_del_f_deg only if your interpolator expects a flap column.
//         // your current header uses no flap column, so only aileron needs overriding there.

//         bool interp_ok = interpolate_aerocoeffs(
//                             throttle_pct,
//                             interp_del_a_deg,    // aileron_left (deg) — forced to 40
//                             del_e_deg,           // elevator (deg) — actual
//                             del_r_deg,           // rudder (deg)   — actual
//                             alpha_deg,           // alpha (deg)
//                             beta_deg,            // beta (deg)
//                             interp_aerocoeffs);

//         if (interp_ok) {
//             // Use interpolated CL, CD, CM
//             c_lift_a = interp_aerocoeffs[0];    // CL
//             c_drag_a = interp_aerocoeffs[1];    // CD
//             CM_pred  = interp_aerocoeffs[2];    // CMm_xcg_1.04
//             used_interpolator_success = true;
//             printf("Interpolator OK (with aileron/flaps=40): CL=%.6f CD=%.6f CM=%.6f\n", c_lift_a, c_drag_a, CM_pred);
//         } else {
//             // fallback handled later (use analytical path)
//             printf("Interpolator failed or out of hull (falling back to analytical model).\n");
//             used_interpolator_success = false;
//         }

//     }

//     // --- If interpolation not used or failed, use analytical method (original code) ---
//     double CL_linear = 0.0;
//     if (!used_interpolator_success) {
//         // compute the analytic linear CL (degrees inputs)
//         CL_linear = compute_CL_from_model(alpha_deg, del_e_deg, del_f_deg, del_a_deg, throttle_for_model, (double)airspeed);
//         // get blended (smooth) CL from the liftCoeff blender (alpha in radians)
//         c_lift_a = liftCoeff(alpha, (float)CL_linear);
//         // compute drag consistently using CL_linear for induced portion
//         c_drag_a = compute_CD_from_model(CL_linear, alpha_deg, del_e_deg, del_f_deg, del_a_deg, throttle_for_model, (double)airspeed);

//         // Compute CM via analytical form (used previously)
//         // // Determine c_m_a (Cm_alpha) like original code
//         // double c_m_a;
//         // if (alpha < 0) {
//         //     c_m_a = 0.06 * radtodeg;
//         // } else if (alpha > 0 && alpha < (5.0 / radtodeg)) {
//         //     c_m_a = 0.15 * radtodeg;
//         // } else {
//         //     c_m_a = -0.1 * radtodeg;
//         // }
//         // Keep existing CM calculation form
//         // m_aero computed later using c_m_0 + c_m_a*alpha + c_m_q*... + c_m_deltae*inputElevator
//         // Here CM_pred left unused for interpolator pathway
//     }

//     // Convert coefficients to body frame (these use c_drag_a and c_lift_a, computed above)
//     double c_x_0 = -c_drag_0 * cos(alpha) + c_lift_0 * sin(alpha);
//     double c_x_a = -c_drag_a * cos(alpha) + c_lift_a * sin(alpha);
//     double c_x_q = -c_drag_q * cos(alpha) + c_lift_q * sin(alpha);
//     double c_z_0 = -c_drag_0 * sin(alpha) - c_lift_0 * cos(alpha);
//     double c_z_a = -c_drag_a * sin(alpha) - c_lift_a * cos(alpha);
//     double c_z_q = -c_drag_q * sin(alpha) - c_lift_q * cos(alpha);

//     // angular rates
//     double p = gyro.x;
//     double q = gyro.y;
//     double r = gyro.z;

//     // Debug print
//     printf("Control Inputs: alpha=%.3f° del_e=%.3f° del_a=%.3f° del_r=%.3f° throttle=%.3f airspeed=%.3f\n",
//            alpha_deg, del_e_deg, del_a_deg, del_r_deg, throttle_for_model, airspeed);

           
//     float thrust_offset = 0.091f;

//     // dynamic pressure and initialize force/moments
//     double qbar = 0.5 * rho * pow(airspeed, 2.0) * s; // dynamic pressure * s
//     float ax = 0.0f, ay = 0.0f, az = 0.0f;   // body forces
//     float la = 0.0f, ma = 0.0f, na = 0.0f;   // body moments

//     if (is_zero(airspeed)) {
//         ax = ay = az = la = ma = na = 0.0f;
//     } else {
//         // aerodynamic forces (body frame)
//         float fx_aero_b = (float)( qbar * (c_x_0 + c_x_a + c_x_q * c * q / (2.0 * airspeed) - c_drag_deltae * cos(alpha) * fabs(inputElevator) + c_lift_deltae * sin(alpha) * inputElevator) );
//         float fy_aero_b = (float)( qbar * (c_y_0 + c_y_b * beta + c_y_p * b * p / (2.0 * airspeed) + c_y_r * b * r / (2.0 * airspeed) + c_y_deltaa * inputAileron + c_y_deltar * inputRudder) );
//         float fz_aero_b = (float)( qbar * (c_z_0 + c_z_a + c_z_q * c * q / (2.0 * airspeed) - c_drag_deltae * sin(alpha) * fabs(inputElevator) - c_lift_deltae * cos(alpha) * inputElevator) );

//         float fz_aero_e = sin(phi) * fx_aero_b + cos(phi) * fz_aero_b;
//         float fz_thrust_e = -thrust * sin(phi);

//         // moments: l, m, n
//         float l_aero;
//         float m_aero;
//         float n_aero;

//         if (used_interpolator_success) {
//             // Use interpolated CM directly to compute moment: m = qbar * c * CM
//             // CL & CD already used for force terms above
//             m_aero = (float)( qbar * c * CM_pred );
//             // roll and yaw moments: use conventional analytical forms for now (since table contained only CL,CD,CM)
//             l_aero = (float)( qbar * b * (c_l_0 + c_l_b * beta + c_l_p * b * p / (2.0 * effective_airspeed) + c_l_r * b * r / (2.0 * effective_airspeed) + c_l_deltaa * inputAileron + c_l_deltar * inputRudder) );
//             n_aero = (float)( qbar * b * (c_n_0 + c_n_b * beta + c_n_p * b * p / (2.0 * effective_airspeed) + c_n_r * b * r / (2.0 * effective_airspeed) + c_n_deltaa * inputAileron + c_n_deltar * inputRudder) );
//         } else {
//             // Analytical moments: keep your original formulas
//             l_aero = (float)( qbar * b * (c_l_0 + c_l_b * beta + c_l_p * b * p / (2.0 * effective_airspeed) + c_l_r * b * r / (2.0 * effective_airspeed) + c_l_deltaa * inputAileron + c_l_deltar * inputRudder) );

//             // compute c_m_a same as earlier
//             double c_m_a;
//             if (alpha < 0) {
//                 c_m_a = 0.06 * radtodeg;
//             } else if (alpha > 0 && alpha < (5.0 / radtodeg)) {
//                 c_m_a = 0.15 * radtodeg;
//             } else {
//                 c_m_a = -0.1 * radtodeg;
//             }
//             m_aero = (float)( qbar * c * (c_m_0 + c_m_a * alpha + c_m_q * c * q / (2.0 * effective_airspeed) + c_m_deltae * inputElevator) );
//             n_aero = (float)( qbar * b * (c_n_0 + c_n_b * beta + c_n_p * b * p / (2.0 * effective_airspeed) + c_n_r * b * r / (2.0 * effective_airspeed) + c_n_deltaa * inputAileron + c_n_deltar * inputRudder) );
//         }

//         float m_thrust = thrust * thrust_offset;

//         // --- Liftoff/landing state machine (unchanged) ---
//         // ... (keep entire state machine exactly as in your original code)
//         // For brevity we will reuse your existing state machine code verbatim below.
//         // (Paste the same state machine content you already have here.)
//         // NOTE: in this snippet we've placed the same state machine content earlier in your code.
//         // The remainder of your original state-machine & normal-force projection code follows.

//         // -------------------- Liftoff/landing state machine --------------------
//         // Tunables
//         const float FORCE_MARGIN         = 0.02f;   // 2% margin above/below mg
//         const int   FORCE_HOLD_FRAMES    = 5;       // consecutive frames to confirm force trigger
//         const float ALT_HYSTERESIS       = 0.010f;  // 1 cm AGL to latch airborne (wheel clearance)
//         const float PEN_ENGAGE           = 0.005f;  // 5 mm engage window (up-positive)
//         const float MIN_LOAD_N           = 0.1f;    // treat tiny normals as zero (reduced for debug)

//         const int   ALT_STABLE_FRAMES    = 5000;       // consecutive frames altitude must be steady
//         const float ALT_STABLE_TOL       = 0.002f;   // meters tolerance (2 mm)
//         const float VZ_TOUCH_THRESH      = 0.6f;    // m/s; vertical speed threshold for landing detection

//         const Vector3f rNose_b(+Lf, 0.0f, frame_height);
//         const Vector3f rMain_b(-Lb, 0.0f, frame_height);
//         const Vector3f rNose_e = dcm.transposed() * rNose_b;   // body -> earth (NED)
//         const Vector3f rMain_e = dcm.transposed() * rMain_b;

//         Vector3f locned;
//         bool alt_gotten = AP::ahrs().get_relative_position_NED_origin_float(locned);
//         const float origin_hagl_up = alt_gotten ? -locned.z : 0.0f;

//         const float wheelNose_world_z = locned.z + rNose_e.z;   // down-positive
//         const float wheelMain_world_z = locned.z + rMain_e.z;   // down-positive

//         static float ground_wheel_z_touch = 0.0f;
//         static bool  ground_wheel_z_latched = false;

//         float wheelNose_AGL_up, wheelMain_AGL_up;
//         if (ground_wheel_z_latched) {
//             wheelNose_AGL_up = ground_wheel_z_touch - wheelNose_world_z; // up-positive
//             wheelMain_AGL_up = ground_wheel_z_touch - wheelMain_world_z;
//         } else {
//             wheelNose_AGL_up = -wheelNose_world_z;
//             wheelMain_AGL_up = -wheelMain_world_z;
//         }
//         wheelNose_AGL_up = std::max(0.0f, wheelNose_AGL_up);
//         wheelMain_AGL_up = std::max(0.0f, wheelMain_AGL_up);
//         const float minWheel_AGL_up = std::min(wheelNose_AGL_up, wheelMain_AGL_up);
//         const float AGL_up = minWheel_AGL_up;

//         const float S1 = mass*GRAVITY_MSS + fz_aero_e + fz_thrust_e; // down-positive
//         float cosph_raw = std::cos(phi);
//         float cosph = (std::fabs(cosph_raw) < 0.02f) ? (copysign(0.02f, cosph_raw)) : cosph_raw;
//         const float S2 = -(m_aero + m_thrust) / cosph;
//         float Nf_prov = ( S2 + Lb*S1 ) / (Lf + Lb);
//         float Nb_prov = ( Lf*S1 - S2 ) / (Lf + Lb);
//         Nf_prov = std::max(0.0f, Nf_prov);
//         Nb_prov = std::max(0.0f, Nb_prov);

//         const float lift_up   = -fz_aero_e;    // up-positive
//         const float thrust_up = -fz_thrust_e;  // up-positive
//         const float Fz_up     = lift_up + thrust_up; // up-positive net upward force

//         static float Fz_up_f = 0.0f;
//         const float EMA = 0.3f;
//         Fz_up_f = (1.0f - EMA) * Fz_up_f + EMA * Fz_up;

//         const bool force_exceeds = (Fz_up_f > mass*GRAVITY_MSS * (1.0f + FORCE_MARGIN));
//         const bool force_below   = (Fz_up_f < mass*GRAVITY_MSS * (1.0f - FORCE_MARGIN));

//         static float alt_last = 0.0f;
//         static int   alt_stable_cnt = 0;
//         static bool has_taken_off = false;
//         if (alt_gotten) {
//             const float dalt = fabsf(origin_hagl_up - alt_last);
//             if (dalt <= ALT_STABLE_TOL) alt_stable_cnt++; else alt_stable_cnt = 0;
//             alt_last = origin_hagl_up;
//         } else {
//             alt_stable_cnt = 0;
//         }

//         enum class AirState : uint8_t { GROUND, FORCE_OK, AIRBORNE };
//         static AirState air_state = AirState::GROUND;
//         static int force_cnt = 0;

//         float Nf = 0.0f, Nb = 0.0f;

//         switch (air_state) {
//         case AirState::GROUND:
//         {
//             if (!ground_wheel_z_latched && alt_gotten) {
//                 ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
//                 ground_wheel_z_latched = true;
//                 std::printf("DBG: initial-ground latch wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//             }
//             Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//             Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//             force_cnt = force_exceeds ? (force_cnt + 1) : 0;
//             if (force_cnt >= FORCE_HOLD_FRAMES) { air_state = AirState::FORCE_OK; force_cnt = 0; }
//             break;
//         }
//         case AirState::FORCE_OK:
//         {
//             Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//             Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//             fm = 1;
//             if (AGL_up >= ALT_HYSTERESIS) {
//                 air_state = AirState::AIRBORNE;
//                 has_taken_off = true;
//                 ground_wheel_z_latched = false;
//                 Nf = Nb = 0.0f;
//                 break;
//             }
//             if (!has_taken_off && force_below) {
//                 air_state = AirState::GROUND;
//                 force_cnt = 0;
//             }
//             break;
//         }
//         case AirState::AIRBORNE:
//         {
//             Nf = Nb = 0.0f;
//             fm = 1;
//             const bool nose_touch = (wheelNose_AGL_up <= PEN_ENGAGE);
//             const bool main_touch = (wheelMain_AGL_up <= PEN_ENGAGE);
//             if (nose_touch || main_touch) {
//                 air_state = AirState::GROUND;
//                 printf("CONDITION1\n");
//                 Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//                 Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//                 ground_wheel_z_touch = main_touch ? wheelMain_world_z : wheelNose_world_z;
//                 ground_wheel_z_latched = true;
//                 has_taken_off = false;
//                 std::printf("DBG: wheel-touch touchdown latched wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//                 force_cnt = 0;
//                 break;
//             }
//             static int touchdown_force_cnt = 0;
//             if (Fz_up_f < 10 && fabsf(velocity_ef.z) < VZ_TOUCH_THRESH && alt_stable_cnt >= ALT_STABLE_FRAMES) touchdown_force_cnt++; else touchdown_force_cnt = 0;
//             if (touchdown_force_cnt >= FORCE_HOLD_FRAMES) {
//                 printf("CONDITION2\n");
//                 air_state = AirState::GROUND;
//                 fm = 0;
//                 Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
//                 Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
//                 ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
//                 ground_wheel_z_latched = true;
//                 has_taken_off = false;
//                 touchdown_force_cnt = 0;
//                 std::printf("DBG: force+stable-alt touchdown latched wheel_z_touch=%.3f (down-pos)\n", ground_wheel_z_touch);
//                 break;
//             }
//             break;
//         }
//         } // end switch

//         // std::printf("state=%d fm=%d AGL=%.3f Fz_up=%.1f Nf=%.1f Nb=%.1f\n",
//         //             int(air_state), int(fm), AGL_up, Fz_up, Nf, Nb);

//         const float fx_norm_b =  (Nf + Nb) * std::sin(phi);
//         const float fz_norm_b = -(Nf + Nb) * std::cos(phi);

//         ax = fx_aero_b + fx_norm_b;   // aero + ground only
//         ay = fy_aero_b;
//         ay = 0;
//         az = fz_aero_b + fz_norm_b;

//         la = l_aero;
//         ma = m_aero + Nb * Lb * std::cos(phi) - Nf * Lf * std::cos(phi) + m_thrust;
//         na = n_aero;
//         la = 0;
//         na = 0;
//     } // end else airspeed non-zero

//     return Wrench{ Vector3f(ax, ay, az), Vector3f(la, ma, na) };
// }



// Python interpolation interface with multiple fallback levels
#include <sstream>
#include <iomanip>
#include <cstdio>

// Static struct and function (file-local only)
struct AeroCoefficients {
    float CL, CD, CMm;
    bool valid;
    AeroCoefficients() : CL(0), CD(0), CMm(0), valid(false) {}
    AeroCoefficients(float cl, float cd, float cmm) : CL(cl), CD(cd), CMm(cmm), valid(true) {}
};

static bool interpolate_aerocoeffs_python
(double throttle_pct, double aileron_deg, double elevator_deg, 
                                         double rudder_deg, double alpha_deg, double beta_deg, 
                                         double interp_aerocoeffs[3]) {
                                            

                                    
        // Static variables for caching and throttling
    static AeroCoefficients last_valid_coeffs(0.5f, 0.08f, 0.0f);
    static int consecutive_failures = 0;
    static bool has_valid_cache = false;
    static uint32_t last_python_call_ms = 0;
    static uint32_t total_calls = 0;
    static uint32_t successful_python_calls = 0;
    
    total_calls++;
    uint32_t current_time_ms = AP_HAL::millis();
    
    // CRITICAL: Rate limit Python calls to prevent blocking
    // Only call Python every 200ms (5Hz max) - SITL runs at 400-1000Hz
    const uint32_t PYTHON_CALL_INTERVAL_MS = 200;
    
    bool should_call_python = (current_time_ms - last_python_call_ms) >= PYTHON_CALL_INTERVAL_MS;
    
    if (!should_call_python) {
        // Use cached values for high-frequency calls (395+ calls out of 400Hz)
        if (has_valid_cache) {
            interp_aerocoeffs[0] = last_valid_coeffs.CL;
            interp_aerocoeffs[1] = last_valid_coeffs.CD;
            interp_aerocoeffs[2] = last_valid_coeffs.CMm;
            return true;
        } else {
            // No cache yet, use simple defaults until first Python call succeeds
            interp_aerocoeffs[0] = 0.2f + 0.08f * alpha_deg;  // Basic CL = CL0 + CLa*alpha
            interp_aerocoeffs[0] = constrain_float(interp_aerocoeffs[0], -0.5f, 1.5f);
            interp_aerocoeffs[1] = 0.05f + 0.3f * interp_aerocoeffs[0] * interp_aerocoeffs[0];  // Basic CD
            interp_aerocoeffs[2] = -0.05f * alpha_deg;  // Basic CMm
            return true;
        }
    }
    
    // Time for a Python call
    last_python_call_ms = current_time_ms;
    
    // Build command with aggressive timeout to prevent blocking
    std::stringstream cmd;
    cmd << "cd /home/ojasj/Ardupilot_KNN/ardupilot && timeout 0.1s python3 aero_interpolator.py " 
        << std::fixed << std::setprecision(3)  // Reduced precision for speed
        << throttle_pct << " " << aileron_deg << " " << elevator_deg << " " 
        << rudder_deg << " " << alpha_deg << " " << beta_deg << " 2>/dev/null";
    
    // Non-blocking Python call with timeout protection
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (pipe) {
        // Set pipe to non-blocking mode
        int fd = fileno(pipe);
        int flags = fcntl(fd, F_GETFL, 0);
        fcntl(fd, F_SETFL, flags | O_NONBLOCK);
        
        char buffer[256];  // Smaller buffer for speed
        std::string result = "";
        
        // Use select() with very short timeout to prevent blocking
        fd_set readfds;
        struct timeval timeout;
        
        FD_ZERO(&readfds);
        FD_SET(fd, &readfds);
        timeout.tv_sec = 0;
        timeout.tv_usec = 50000;  // 50ms max wait - critical for real-time performance
        
        int select_result = select(fd + 1, &readfds, NULL, NULL, &timeout);
        
        if (select_result > 0 && FD_ISSET(fd, &readfds)) {
            // Data available, read quickly
            size_t bytes_read = fread(buffer, 1, sizeof(buffer) - 1, pipe);
            if (bytes_read > 0) {
                buffer[bytes_read] = '\0';
                result = buffer;
            }
        } else if (select_result == 0) {
            // Timeout - Python took too long
            consecutive_failures++;
        }
        
        pclose(pipe);
        
        // Fast parsing - only look for SUCCESS and extract numbers
        if (result.find("SUCCESS") != std::string::npos) {
            size_t json_start = result.find("{");
            if (json_start != std::string::npos && json_start + 50 < result.length()) {
                float CL = 0, CD = 0, CMm = 0;
                const char* json_str = result.c_str() + json_start;
                
                // Fast JSON parsing using sscanf
                if (sscanf(json_str, "{\"CL\": %f, \"CD\": %f, \"CMm\": %f}", &CL, &CD, &CMm) == 3) {
                    // Sanity check values
                    if (CL >= -2.0f && CL <= 3.0f && CD >= 0.0f && CD <= 1.0f && 
                        CMm >= -1.0f && CMm <= 1.0f) {
                        
                        // Success! Update cache
                        last_valid_coeffs = AeroCoefficients(CL, CD, CMm);
                        has_valid_cache = true;
                        consecutive_failures = 0;
                        successful_python_calls++;
                        
                        interp_aerocoeffs[0] = CL;
                        interp_aerocoeffs[1] = CD;
                        interp_aerocoeffs[2] = CMm;
                        
                        // Periodic status (every 25 successful calls = ~5 seconds at 5Hz)
                        if (successful_python_calls % 25 == 0) {
                            printf("Python Interp: CL=%.3f CD=%.3f CMm=%.3f (success rate: %d/%d)\n", 
                                   CL, CD, CMm, successful_python_calls, total_calls);
                        }
                        
                        return true;
                    }
                }
            }
        }
    }
    
    // Python failed - use fallback strategy
    consecutive_failures++;
    
    // Level 1: Use cached values with gradual decay (up to 100 failures = 20 seconds)
    if (has_valid_cache && consecutive_failures < 100) {
        // Apply very gentle decay toward safe values
        float decay_factor = std::min(0.1f, consecutive_failures * 0.001f);  // Max 10% decay
        
        float safe_CL = 0.2f + 0.08f * alpha_deg;
        float safe_CD = 0.05f + 0.3f * safe_CL * safe_CL;
        float safe_CMm = -0.05f * alpha_deg;
        
        float CL_blend = last_valid_coeffs.CL * (1.0f - decay_factor) + safe_CL * decay_factor;
        float CD_blend = last_valid_coeffs.CD * (1.0f - decay_factor) + safe_CD * decay_factor;
        float CMm_blend = last_valid_coeffs.CMm * (1.0f - decay_factor) + safe_CMm * decay_factor;
        
        interp_aerocoeffs[0] = CL_blend;
        interp_aerocoeffs[1] = CD_blend;
        interp_aerocoeffs[2] = CMm_blend;
        
        if (consecutive_failures % 20 == 0) {  // Every 4 seconds
            printf("Python Interp CACHED (decay=%.1f%%, failures=%d): CL=%.3f CD=%.3f CMm=%.3f\n", 
                   decay_factor * 100, consecutive_failures, CL_blend, CD_blend, CMm_blend);
        }
        
        return true;
    }
    
    // Level 2: Ultimate fallback - physics-based defaults
    float safe_CL = 0.2f + 0.08f * alpha_deg;  // CL = CL0 + CLa*alpha
    safe_CL = constrain_float(safe_CL, -0.5f, 1.5f);
    
    float safe_CD = 0.05f + 0.5f * safe_CL * safe_CL;  // CD = CD0 + CDi
    safe_CD = constrain_float(safe_CD, 0.03f, 0.3f);
    
    float safe_CMm = -0.05f * alpha_deg;  // Basic pitch stability
    safe_CMm = constrain_float(safe_CMm, -0.2f, 0.2f);
    
    interp_aerocoeffs[0] = safe_CL;
    interp_aerocoeffs[1] = safe_CD;
    interp_aerocoeffs[2] = safe_CMm;
    
    if (consecutive_failures % 20 == 0) {
        printf("Python Interp FALLBACK (failures=%d): CL=%.3f CD=%.3f CMm=%.3f\n", 
               consecutive_failures, safe_CL, safe_CD, safe_CMm);
    }
    
    return true;
}

SITL::Wrench Plane::getForcesAndMoments(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const struct sitl_input &input, bool fm)
{
    const float phi = AP::ahrs().get_pitch();
    float alpha = angle_of_attack;
    const double radtodeg = 180.0 / M_PI;
    
    printf("alpha = %.3f°\n", alpha * radtodeg);
    
    float effective_airspeed = airspeed;
    if (tailsitter || aerobatic) {
        effective_airspeed += inputThrust * 20;
        alpha *= constrain_float(1 - inputThrust, 0, 1);
    }
    
    // Get aircraft parameters
    const float s = coefficient.s;
    const float c = coefficient.c;
    // const float b = coefficient.b;
    const float Lf = 0.858f;    // CG to nose gear (+x)
    const float Lb = 0.122f;    // CG to main gear (aft)
    float rho = air_density;
    
    // Get throttle
    float throttle;
    if (reverse_thrust) {
        throttle = filtered_servo_angle(input, 2);
    } else {
        throttle = filtered_servo_range(input, 2);
    }
    float thrust = throttle * thrust_scale;
    
    // Convert to degrees for interpolation
    double alpha_deg = alpha * radtodeg;
    double elevator_deg =-inputElevator * radtodeg;
    double aileron_deg = -inputAileron * radtodeg;
    double rudder_deg = -inputRudder * radtodeg;
    double beta_deg = beta * radtodeg;
    double throttle_pct = constrain_float(throttle * 100.0f, 10.0f, 100.0f);  // 10-100% range
    
    // Constrain inputs to your CFD data ranges
    alpha_deg = constrain_float(alpha_deg, -3.0, 9.0);
    elevator_deg = constrain_float(elevator_deg, -15.0, 15.0);
    
    printf("Inputs: α=%.2f° δe=%.2f° δa=%.2f° δr=%.2f° thr=%.1f%% β=%.2f°\n",
           alpha_deg, elevator_deg, aileron_deg, rudder_deg, throttle_pct, beta_deg);
    
    // Get aerodynamic coefficients via Python interpolation
    double aero_coeffs[3];
    bool interp_success = interpolate_aerocoeffs_python(
        throttle_pct, aileron_deg, elevator_deg, rudder_deg, alpha_deg, beta_deg, aero_coeffs);
    
    if(interp_success) {
        printf("Using interpolated coefficients: CL=%.4f CD=%.4f CMm=%.4f\n", 
               aero_coeffs[0], aero_coeffs[1], aero_coeffs[2]);
    } else {
        printf("Interpolation failed completely, using safe defaults.\n");
    }
    // Extract coefficients
    float CL = aero_coeffs[0];
    float CD = aero_coeffs[1]; 
    double q = gyro.y;
    float CMm = aero_coeffs[2] - 95.2255*q;
    
    // Calculate dynamic pressure
    double qbar = 0.5 * rho * airspeed * airspeed * s;
    float thrust_offset = 0.0;
    
    // Initialize forces and moments
    float ax = 0.0f, ay = 0.0f, az = 0.0f;
    float la = 0.0f, ma = 0.0f, na = 0.0f;
    
    if (!is_zero(airspeed)) {
        // Calculate aerodynamic forces (body frame)
        float fx_aero = -CD * qbar;  // Drag (negative X)
        // float fy_aero = 0.0f;        // No side force (longitudinal model)
        float fz_aero = -CL * qbar;  // Lift (negative Z in NED body frame)
        
        // Calculate aerodynamic moments
        // float l_aero = 0.0f;                        // No roll moment (longitudinal model)
        float m_aero = CMm * qbar * c;              // Pitch moment from interpolation
        // float n_aero = 0.0f;                        // No yaw moment (longitudinal model)
        
        float fz_aero_e = sin(phi) * fx_aero + cos(phi) * fz_aero;
        float fz_thrust_e = -thrust * sin(phi);
        float m_thrust = thrust * thrust_offset;
        
        // Landing gear ground reaction forces (keep your existing state machine)
        const float FORCE_MARGIN = 0.02f;
        const int FORCE_HOLD_FRAMES = 5;
        const float ALT_HYSTERESIS = 0.010f;
        const float PEN_ENGAGE = 0.005f;
        const float MIN_LOAD_N = 0.1f;
        const int ALT_STABLE_FRAMES = 5000;
        const float ALT_STABLE_TOL = 0.002f;
        const float VZ_TOUCH_THRESH = 0.6f;
        
        const Vector3f rNose_b(+Lf, 0.0f, frame_height);
        const Vector3f rMain_b(-Lb, 0.0f, frame_height);
        const Vector3f rNose_e = dcm.transposed() * rNose_b;
        const Vector3f rMain_e = dcm.transposed() * rMain_b;
        
        Vector3f locned;
        bool alt_gotten = AP::ahrs().get_relative_position_NED_origin_float(locned);
        const float origin_hagl_up = alt_gotten ? -locned.z : 0.0f;
        
        const float wheelNose_world_z = locned.z + rNose_e.z;
        const float wheelMain_world_z = locned.z + rMain_e.z;
        
        static float ground_wheel_z_touch = 0.0f;
        static bool ground_wheel_z_latched = false;
        
        float wheelNose_AGL_up, wheelMain_AGL_up;
        if (ground_wheel_z_latched) {
            wheelNose_AGL_up = ground_wheel_z_touch - wheelNose_world_z;
            wheelMain_AGL_up = ground_wheel_z_touch - wheelMain_world_z;
        } else {
            wheelNose_AGL_up = -wheelNose_world_z;
            wheelMain_AGL_up = -wheelMain_world_z;
        }
        wheelNose_AGL_up = std::max(0.0f, wheelNose_AGL_up);
        wheelMain_AGL_up = std::max(0.0f, wheelMain_AGL_up);
        const float minWheel_AGL_up = std::min(wheelNose_AGL_up, wheelMain_AGL_up);
        const float AGL_up = minWheel_AGL_up;
        
        const float S1 = mass*GRAVITY_MSS + fz_aero_e + fz_thrust_e;
        float cosph_raw = std::cos(phi);
        float cosph = (std::fabs(cosph_raw) < 0.02f) ? (copysign(0.02f, cosph_raw)) : cosph_raw;
        const float S2 = -(m_aero + m_thrust) / cosph;
        float Nf_prov = (S2 + Lb*S1) / (Lf + Lb);
        float Nb_prov = (Lf*S1 - S2) / (Lf + Lb);
        Nf_prov = std::max(0.0f, Nf_prov);
        Nb_prov = std::max(0.0f, Nb_prov);
        
        const float lift_up = -fz_aero_e;
        const float thrust_up = -fz_thrust_e;
        const float Fz_up = lift_up + thrust_up;
        
        static float Fz_up_f = 0.0f;
        const float EMA = 0.3f;
        Fz_up_f = (1.0f - EMA) * Fz_up_f + EMA * Fz_up;
        
        const bool force_exceeds = (Fz_up_f > mass*GRAVITY_MSS * (1.0f + FORCE_MARGIN));
        const bool force_below = (Fz_up_f < mass*GRAVITY_MSS * (1.0f - FORCE_MARGIN));
        
        static float alt_last = 0.0f;
        static int alt_stable_cnt = 0;
        static bool has_taken_off = false;
        if (alt_gotten) {
            const float dalt = fabsf(origin_hagl_up - alt_last);
            if (dalt <= ALT_STABLE_TOL) alt_stable_cnt++; else alt_stable_cnt = 0;
            alt_last = origin_hagl_up;
        } else {
            alt_stable_cnt = 0;
        }
        
        enum class AirState : uint8_t { GROUND, FORCE_OK, AIRBORNE };
        static AirState air_state = AirState::GROUND;
        static int force_cnt = 0;
        float Nf = 0.0f, Nb = 0.0f;
        
        // State machine (unchanged from your original)
        switch (air_state) {
        case AirState::GROUND: {
            if (!ground_wheel_z_latched && alt_gotten) {
                ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
                ground_wheel_z_latched = true;
            }
            Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
            Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
            force_cnt = force_exceeds ? (force_cnt + 1) : 0;
            if (force_cnt >= FORCE_HOLD_FRAMES) { air_state = AirState::FORCE_OK; force_cnt = 0; }
            break;
        }
        case AirState::FORCE_OK: {
            Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
            Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
            fm = 1;
            if (AGL_up >= ALT_HYSTERESIS) {
                air_state = AirState::AIRBORNE;
                has_taken_off = true;
                ground_wheel_z_latched = false;
                Nf = Nb = 0.0f;
                break;
            }
            if (!has_taken_off && force_below) {
                air_state = AirState::GROUND;
                force_cnt = 0;
            }
            break;
        }
        case AirState::AIRBORNE: {
            Nf = Nb = 0.0f;
            fm = 1;
            const bool nose_touch = (wheelNose_AGL_up <= PEN_ENGAGE);
            const bool main_touch = (wheelMain_AGL_up <= PEN_ENGAGE);
            if (nose_touch || main_touch) {
                air_state = AirState::GROUND;
                Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
                Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
                ground_wheel_z_touch = main_touch ? wheelMain_world_z : wheelNose_world_z;
                ground_wheel_z_latched = true;
                has_taken_off = false;
                force_cnt = 0;
                break;
            }
            static int touchdown_force_cnt = 0;
            if (Fz_up_f < 10 && fabsf(velocity_ef.z) < VZ_TOUCH_THRESH && alt_stable_cnt >= ALT_STABLE_FRAMES) touchdown_force_cnt++; else touchdown_force_cnt = 0;
            if (touchdown_force_cnt >= FORCE_HOLD_FRAMES) {
                air_state = AirState::GROUND;
                fm = 0;
                Nf = (Nf_prov > MIN_LOAD_N) ? Nf_prov : 0.0f;
                Nb = (Nb_prov > MIN_LOAD_N) ? Nb_prov : 0.0f;
                ground_wheel_z_touch = std::min(wheelNose_world_z, wheelMain_world_z);
                ground_wheel_z_latched = true;
                has_taken_off = false;
                touchdown_force_cnt = 0;
                break;
            }
            break;
        }
        }
        
        // Calculate ground reaction forces
        const float fx_norm_b = (Nf + Nb) * std::sin(phi);
        const float fz_norm_b = -(Nf + Nb) * std::cos(phi);
        
        // Final forces and moments (simplified - longitudinal only)
        ax = fx_aero + fx_norm_b;
        ay = 0.0f;  // No lateral force
        az = fz_aero + fz_norm_b;
        
        la = 0.0f;  // No roll moment
        ma = m_aero + Nb * Lb * std::cos(phi) - Nf * Lf * std::cos(phi) + m_thrust;
        na = 0.0f;  // No yaw moment
    }
    
    return Wrench{Vector3f(ax, ay, az), Vector3f(la, ma, na)};
}




void Plane::calculate_forces(const struct sitl_input &input, Vector3f &rot_accel, bool fm)
{
    float aileron  = filtered_servo_angle(input, 0);
    float elevator = filtered_servo_angle(input, 1);
    float rudder   = filtered_servo_angle(input, 3);
    bool launch_triggered = input.servos[6] > 1700;
    float throttle;
    if (reverse_elevator_rudder) {
        elevator = -elevator;
        rudder = -rudder;
    }
    if (elevons) {
        // fake an elevon plane
        float ch1 = aileron;
        float ch2 = elevator;
        aileron  = (ch2-ch1)/2.0f;
        // the minus does away with the need for RC2_REVERSED=-1
        elevator = -(ch2+ch1)/2.0f;

        // assume no rudder
        rudder = 0;
    } else if (vtail) {
        // fake a vtail plane
        float ch1 = elevator;
        float ch2 = rudder;
        // this matches VTAIL_OUTPUT==2
        elevator = (ch2-ch1)/2.0f;
        rudder   = (ch2+ch1)/2.0f;
    } else if (dspoilers) {
        // fake a differential spoiler plane. Use outputs 1, 2, 4 and 5
        float dspoiler1_left = filtered_servo_angle(input, 0);
        float dspoiler1_right = filtered_servo_angle(input, 1);
        float dspoiler2_left = filtered_servo_angle(input, 3);
        float dspoiler2_right = filtered_servo_angle(input, 4);
        float elevon_left  = (dspoiler1_left + dspoiler2_left)/2;
        float elevon_right = (dspoiler1_right + dspoiler2_right)/2;
        aileron  = (elevon_right-elevon_left)/2;
        elevator = (elevon_left+elevon_right)/2;
        rudder = fabsf(dspoiler1_right - dspoiler2_right)/2 - fabsf(dspoiler1_left - dspoiler2_left)/2;
    } else if (redundant) {
        // channels 1/9 are left/right ailierons
        // channels 2/10 are left/right elevators
        // channels 4/12 are top/bottom rudders
        aileron  = (filtered_servo_angle(input, 0) + filtered_servo_angle(input, 8)) / 2.0;
        elevator = (filtered_servo_angle(input, 1) + filtered_servo_angle(input, 9)) / 2.0;
        rudder   = (filtered_servo_angle(input, 3) + filtered_servo_angle(input, 11)) / 2.0;
    }
    //printf("Aileron: %.1f elevator: %.1f rudder: %.1f\n", aileron, elevator, rudder);

    if (reverse_thrust) {
        throttle = filtered_servo_angle(input, 2);
    } else {
        throttle = filtered_servo_range(input, 2);
    }
    
    float thrust     = throttle;

    battery_voltage = sitl->batt_voltage - 0.7*throttle;
    battery_current = (battery_voltage/sitl->batt_voltage)*50.0f*sq(throttle);

    if (ice_engine) {
        thrust = icengine.update(input);
    }

    // calculate angle of attack
    angle_of_attack = atan2f(velocity_air_bf.z, velocity_air_bf.x);
    beta = atan2f(velocity_air_bf.y,velocity_air_bf.x);

    if (tailsitter || aerobatic) {
        /*
          tailsitters get 4x the control surfaces
         */
        aileron *= 4;
        elevator *= 4;
        rudder *= 4;
    }
    
     // simulate engine RPM
    motor_mask |= (1U<<2);
    rpm[2] = thrust * 7000;

    // scale thrust to newtons
    thrust *= thrust_scale;
    //float thrust_offset = 0.091;

    // Vector3f force = getForce(aileron, elevator, rudder);
    // rot_accel = getTorque(aileron, elevator, rudder, thrust, force);

    //rot_accel[1] = rot_accel[1]+thrust*thrust_offset;
    Wrench w = getForcesAndMoments(aileron, elevator, rudder, thrust, input,fm);

    Vector3f force = w.F;
    rot_accel = w.M;
    rot_accel[0] = rot_accel[0]/3.0f;
    rot_accel[1] = rot_accel[1]/15.0f;
    rot_accel[2] = rot_accel[2]/5.0f;

    if (have_launcher) {
        /*
          simple simulation of a launcher
         */
        if (launch_triggered) {
            uint64_t now = AP_HAL::millis64();
            if (launch_start_ms == 0) {
                launch_start_ms = now;
            }
            if (now - launch_start_ms < launch_time*1000) {
                force.x += mass * launch_accel;
                force.z += mass * launch_accel/3;
            }
        } else {
            // allow reset of catapult
            launch_start_ms = 0;
        }
    }

    

    accel_body = Vector3f(thrust, 0, 0) + force;
    accel_body /= mass;

    // add some noise
    if (thrust_scale > 0) {
        add_noise(fabsf(thrust) / thrust_scale);
    }

    if (on_ground() && !tailsitter) {
        // add some ground friction
        Vector3f vel_body = dcm.transposed() * velocity_ef;
        accel_body.x -= vel_body.x * 0.3f;
    }
}
    
/*
  update the plane simulation by one time step
 */
void Plane::update(const struct sitl_input &input)
{
    Vector3f rot_accel;

    bool flightmode = 0;

    update_wind(input);
    
    calculate_forces(input, rot_accel,flightmode);
    
    update_dynamics(rot_accel);

    /*
      add in ground steering, this should be replaced with a proper
      calculation of a nose wheel effect
    */
    if (have_steering && on_ground()) {
        const float steering = filtered_servo_angle(input, 4);
        const Vector3f velocity_bf = dcm.transposed() * velocity_ef;
        const float steer_scale = radians(5);
        gyro.z += steering * velocity_bf.x * steer_scale;
    }

    update_external_payload(input);

    // update lat/lon/altitude
    update_position();
    time_advance();

    // update magnetic field
    update_mag_field_bf();
}















// Torque calculation function
// Vector3f Plane::getTorque(float inputAileron, float inputElevator, float inputRudder, float inputThrust, const Vector3f &force) const
// {
//     float alpha = angle_of_attack;

// 	//calculate aerodynamic torque
//     float effective_airspeed = airspeed;

//     if (tailsitter || aerobatic) {
//         /*
//           tailsitters get airspeed from prop-wash
//          */
//         effective_airspeed += inputThrust * 20;

//         // reduce effective angle of attack as thrust increases
//         alpha *= constrain_float(1 - inputThrust, 0, 1);
//     }
    
//     const float s = coefficient.s;
//     const float c = coefficient.c;
//     const float b = coefficient.b;
//     const float c_l_0 = coefficient.c_l_0;
//     const float c_l_b = coefficient.c_l_b;
//     const float c_l_p = coefficient.c_l_p;
//     const float c_l_r = coefficient.c_l_r;
//     const float c_l_deltaa = coefficient.c_l_deltaa;
//     const float c_l_deltar = coefficient.c_l_deltar;
//     const float c_m_0 = coefficient.c_m_0;
//     const float c_m_a = coefficient.c_m_a;
//     const float c_m_q = coefficient.c_m_q;
//     const float c_m_deltae = coefficient.c_m_deltae;
//     const float c_n_0 = coefficient.c_n_0;
//     const float c_n_b = coefficient.c_n_b;
//     const float c_n_p = coefficient.c_n_p;
//     const float c_n_r = coefficient.c_n_r;
//     const float c_n_deltaa = coefficient.c_n_deltaa;
//     const float c_n_deltar = coefficient.c_n_deltar;
//     const Vector3f &CGOffset = coefficient.CGOffset;
    
//     float rho = air_density;

// 	//read angular rates
// 	double p = gyro.x;
// 	double q = gyro.y;
// 	double r = gyro.z;

// 	double qbar = 1.0/2.0*rho*pow(effective_airspeed,2)*s; //Calculate dynamic pressure
// 	double la, na, ma;
// 	if (is_zero(effective_airspeed))
// 	{
// 		la = 0;
// 		ma = 0;
// 		na = 0;
// 	}
// 	else
// 	{
// 		la = qbar*b*(c_l_0 + c_l_b*beta + c_l_p*b*p/(2*effective_airspeed) + c_l_r*b*r/(2*effective_airspeed) + c_l_deltaa*inputAileron + c_l_deltar*inputRudder);
// 		ma = qbar*c*(c_m_0 + c_m_a*alpha + c_m_q*c*q/(2*effective_airspeed) + c_m_deltae*inputElevator);
// 		na = qbar*b*(c_n_0 + c_n_b*beta + c_n_p*b*p/(2*effective_airspeed) + c_n_r*b*r/(2*effective_airspeed) + c_n_deltaa*inputAileron + c_n_deltar*inputRudder);
// 	}


// 	// Add torque to force misalignment with CG
// 	// r x F, where r is the distance from CoG to CoL
// 	la +=  CGOffset.y * force.z - CGOffset.z * force.y;
// 	ma += -CGOffset.x * force.z + CGOffset.z * force.x;
// 	na += -CGOffset.y * force.x + CGOffset.x * force.y;

// 	return Vector3f(la, ma, na);
// }

// // Force calculation function from last_letter
// Vector3f Plane::getForce(float inputAileron, float inputElevator, float inputRudder) const
// {
//     const float alpha = angle_of_attack;
//     const float c_drag_q = coefficient.c_drag_q;
//     const float c_lift_q = coefficient.c_lift_q;
//     const float s = coefficient.s;
//     const float c = coefficient.c;
//     const float b = coefficient.b;
//     const float c_drag_deltae = coefficient.c_drag_deltae;
//     const float c_lift_deltae = coefficient.c_lift_deltae;
//     const float c_y_0 = coefficient.c_y_0;
//     const float c_y_b = coefficient.c_y_b;
//     const float c_y_p = coefficient.c_y_p;
//     const float c_y_r = coefficient.c_y_r;
//     const float c_y_deltaa = coefficient.c_y_deltaa;
//     const float c_y_deltar = coefficient.c_y_deltar;
//     const float c_drag_0 = coefficient.c_drag_0;
//     const float c_lift_0 = coefficient.c_lift_0;
    


//     float rho = air_density;

// 	//request lift and drag alpha-coefficients from the corresponding functions
// 	double c_lift_a = liftCoeff(alpha);
// 	double c_drag_a = dragCoeff(alpha);

// 	//convert coefficients to the body frame
//     double c_x_0 = -c_drag_0*cos(alpha)+c_lift_0*sin(alpha);
// 	double c_x_a = -c_drag_a*cos(alpha)+c_lift_a*sin(alpha);
// 	double c_x_q = -c_drag_q*cos(alpha)+c_lift_q*sin(alpha);
//     double c_z_0 = -c_drag_0*sin(alpha)-c_lift_0*cos(alpha);
// 	double c_z_a = -c_drag_a*sin(alpha)-c_lift_a*cos(alpha);
// 	double c_z_q = -c_drag_q*sin(alpha)-c_lift_q*cos(alpha);

// 	//read angular rates
// 	double p = gyro.x;
// 	double q = gyro.y;
// 	double r = gyro.z;

// 	//calculate aerodynamic force
// 	double qbar = 1.0/2.0*rho*pow(airspeed,2)*s; //Calculate dynamic pressure
// 	double ax, ay, az;
// 	if (is_zero(airspeed))
// 	{
// 		ax = 0;
// 		ay = 0;
// 		az = 0;
// 	}
// 	else
// 	{
// 		ax = qbar*(c_x_0 + c_x_a + c_x_q*c*q/(2*airspeed) - c_drag_deltae*cos(alpha)*fabs(inputElevator) + c_lift_deltae*sin(alpha)*inputElevator);
// 		// split c_x_deltae to include "abs" term
// 		ay = qbar*(c_y_0 + c_y_b*beta + c_y_p*b*p/(2*airspeed) + c_y_r*b*r/(2*airspeed) + c_y_deltaa*inputAileron + c_y_deltar*inputRudder);
// 		az = qbar*(c_z_0 + c_z_a + c_z_q*c*q/(2*airspeed) - c_drag_deltae*sin(alpha)*fabs(inputElevator) - c_lift_deltae*cos(alpha)*inputElevator);
// 		// split c_z_deltae to include "abs" term
// 	}
//     return Vector3f(ax, ay, az);
// }