import sys
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay
import json

def barycentric_interpolate_with_lambdas(point, tri, values_train, vertices):
    """
    Your existing barycentric interpolation function - unchanged
    """
    simplex_coords = tri.points[vertices]
    T = np.vstack((simplex_coords.T, np.ones((1, simplex_coords.shape[0]))))
    y = np.append(point, 1.0)
    lambdas = np.linalg.solve(T, y)
    interpolated_values = np.dot(lambdas, values_train[vertices])
    return interpolated_values, lambdas

class AeroInterpolator:
    def __init__(self, csv_file="Aero_Data.csv"):
        """
        Load CSV data with columns: Throttle,Aileron_left,Elevator,Rudder,Alpha,Beta,CL,CD,CMm_xcg_1.04,Cmyu
        """
        print(f"Loading aerodynamic data from {csv_file}")
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} rows from CSV")
        
        # Input columns: Throttle, Aileron_left, Elevator, Rudder, Alpha, Beta (indices 0,1,2,3,4,5)
        self.X = df.iloc[:, [0,1,2,3,4,5]].values  # 6D input space
        
        # Output columns: CL, CD, CMm_xcg_1.04 (indices 6,7,8) - ignore Cmyu (index 9)
        self.Y = df.iloc[:, [6,7,8]].values  # 3D output space
        
        print(f"Input space shape: {self.X.shape}")
        print(f"Output space shape: {self.Y.shape}")
        
        # Build Delaunay triangulation
        try:
            self.tri = Delaunay(self.X)
            print(f"Successfully built Delaunay triangulation with {len(self.tri.simplices)} simplices")
        except Exception as e:
            print(f"Error building Delaunay triangulation: {e}")
            raise
    
    def interpolate(self, throttle, aileron_left, elevator, rudder, alpha, beta):
        """
        Interpolate aerodynamic coefficients
        
        Args:
            throttle: Throttle setting
            aileron_left: Left aileron deflection 
            elevator: Elevator deflection
            rudder: Rudder deflection
            alpha: Angle of attack (degrees)
            beta: Sideslip angle (degrees)
            
        Returns:
            tuple: (CL, CD, CMm) or None if outside convex hull
        """
        # Create 6D input vector
        pt = np.array([throttle, aileron_left, elevator, rudder, alpha, beta])
        
        # Find simplex containing this point
        simplex = self.tri.find_simplex(pt)
        if simplex == -1:
            print(f"Point {pt} is outside convex hull")
            return None
        
        # Get vertices of the simplex
        verts_local = self.tri.simplices[simplex]
        
        try:
            # Perform barycentric interpolation
            pred, lambdas = barycentric_interpolate_with_lambdas(pt, self.tri, self.Y, verts_local)
            
            # Check if lambdas are reasonable (for debugging)
            lam_out_of_bounds = np.any((lambdas < -1e-8) | (lambdas > 1.0 + 1e-8))
            if lam_out_of_bounds:
                print(f"Warning: Some barycentric coordinates out of bounds: {lambdas}")
            
            return pred[0], pred[1], pred[2]  # CL, CD, CMm
            
        except np.linalg.LinAlgError as e:
            print(f"Linear algebra error during interpolation: {e}")
            return None

def main():
    if len(sys.argv) != 7:
        print("Usage: python3 aero_interpolator.py <throttle> <aileron_left> <elevator> <rudder> <alpha> <beta>")
        print("Example: python3 aero_interpolator.py 50.0 10.0 5.0 2.0 8.0 1.0")
        sys.exit(1)
    
    try:
        # Parse command line arguments
        throttle = float(sys.argv[1])
        aileron_left = float(sys.argv[2])
        elevator = float(sys.argv[3])
        rudder = float(sys.argv[4])
        alpha = float(sys.argv[5])
        beta = float(sys.argv[6])
        
        print(f"Input: throttle={throttle}, aileron_left={aileron_left}, elevator={elevator}, rudder={rudder}, alpha={alpha}, beta={beta}")
        
        # Create interpolator
        interpolator = AeroInterpolator("Aero_Data.csv")
        
        # Perform interpolation
        result = interpolator.interpolate(throttle, aileron_left, elevator, rudder, alpha, beta)
        
        if result is None:
            print("ERROR: Point outside convex hull or interpolation failed")
            sys.exit(1)
        
        # Output results as JSON for easy parsing in C++
        output = {
            "CL": float(result[0]),
            "CD": float(result[1]), 
            "CMm": float(result[2])
        }
        
        print("SUCCESS")
        print(json.dumps(output))
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()



## Previous version with never-fail logic - kept for reference but not used

# import sys
# import numpy as np
# import pandas as pd
# from scipy.spatial import Delaunay
# from scipy.spatial.distance import cdist
# import json


# def barycentric_interpolate_with_lambdas(point, tri, values_train, vertices):
#     """
#     Your existing barycentric interpolation function - unchanged
#     """
#     simplex_coords = tri.points[vertices]
#     T = np.vstack((simplex_coords.T, np.ones((1, simplex_coords.shape[0]))))
#     y = np.append(point, 1.0)
#     lambdas = np.linalg.solve(T, y)
#     interpolated_values = np.dot(lambdas, values_train[vertices])
#     return interpolated_values, lambdas


# class AeroInterpolator:
#     def __init__(self, csv_file="Aero_Data.csv"):
#         """
#         Load CSV data with columns: Throttle,Aileron_left,Elevator,Rudder,Alpha,Beta,CL,CD,CMm_xcg_1.04,Cmyu
#         """
#         print(f"Loading aerodynamic data from {csv_file}")
#         df = pd.read_csv(csv_file)
#         print(f"Loaded {len(df)} rows from CSV")
        
#         # Input columns: Throttle, Aileron_left, Elevator, Rudder, Alpha, Beta (indices 0,1,2,3,4,5)
#         self.X = df.iloc[:, [0,1,2,3,4,5]].values  # 6D input space
        
#         # Output columns: CL, CD, CMm_xcg_1.04 (indices 6,7,8) - ignore Cmyu (index 9)
#         self.Y = df.iloc[:, [6,7,8]].values  # 3D output space
        
#         print(f"Input space shape: {self.X.shape}")
#         print(f"Output space shape: {self.Y.shape}")
        
#         # Compute data bounds and centroid for fallback methods
#         self.data_min = self.X.min(axis=0)
#         self.data_max = self.X.max(axis=0)
#         self.data_centroid = self.X.mean(axis=0)
#         self.data_std = self.X.std(axis=0)
        
#         print(f"Data ranges: throttle=[{self.data_min[0]:.1f}, {self.data_max[0]:.1f}], "
#               f"aileron=[{self.data_min[1]:.1f}, {self.data_max[1]:.1f}], "
#               f"elevator=[{self.data_min[2]:.1f}, {self.data_max[2]:.1f}], "
#               f"rudder=[{self.data_min[3]:.1f}, {self.data_max[3]:.1f}], "
#               f"alpha=[{self.data_min[4]:.1f}, {self.data_max[4]:.1f}], "
#               f"beta=[{self.data_min[5]:.1f}, {self.data_max[5]:.1f}]")
        
#         # Build Delaunay triangulation
#         try:
#             self.tri = Delaunay(self.X)
#             print(f"Successfully built Delaunay triangulation with {len(self.tri.simplices)} simplices")
#         except Exception as e:
#             print(f"Error building Delaunay triangulation: {e}")
#             raise
    
#     def interpolate_with_fallbacks(self, throttle, aileron_left, elevator, rudder, alpha, beta):
#         """
#         Never-fail interpolation with progressive fallback strategies
#         """
#         original_point = np.array([throttle, aileron_left, elevator, rudder, alpha, beta])
        
#         # Strategy 1: Try original point (barycentric interpolation)
#         result = self._try_barycentric_interpolation(original_point)
#         if result is not None:
#             return result, "barycentric"
        
#         # Strategy 2: Smart perturbation toward data centroid
#         result = self._try_centroid_perturbation(original_point)
#         if result is not None:
#             return result, "centroid_perturbation"
        
#         # Strategy 3: Clamp to data bounds and interpolate
#         result = self._try_boundary_clamping(original_point)
#         if result is not None:
#             return result, "boundary_clamping"
        
#         # Strategy 4: Nearest neighbor approach (always works)
#         result = self._nearest_neighbor_fallback(original_point)
#         return result, "nearest_neighbor"
    
#     def _try_barycentric_interpolation(self, point):
#         """Strategy 1: Standard barycentric interpolation"""
#         try:
#             simplex = self.tri.find_simplex(point)
#             if simplex == -1:
#                 return None
            
#             verts_local = self.tri.simplices[simplex]
#             pred, lambdas = barycentric_interpolate_with_lambdas(point, self.tri, self.Y, verts_local)
            
#             # Check if lambdas are reasonable
#             lam_out_of_bounds = np.any((lambdas < -1e-6) | (lambdas > 1.0 + 1e-6))
#             if lam_out_of_bounds:
#                 return None
            
#             return pred[0], pred[1], pred[2]  # CL, CD, CMm
            
#         except Exception:
#             return None
    
#     def _try_centroid_perturbation(self, original_point):
#         """Strategy 2: Gradually move point toward data centroid until inside hull"""
#         current_point = original_point.copy()
#         direction = self.data_centroid - original_point
        
#         # Try different step sizes
#         for step_size in [0.1, 0.2, 0.3, 0.4, 0.5]:
#             perturbed_point = original_point + step_size * direction
#             result = self._try_barycentric_interpolation(perturbed_point)
#             if result is not None:
#                 return result
        
#         return None
    
#     def _try_boundary_clamping(self, original_point):
#         """Strategy 3: Clamp point to data bounds and try interpolation"""
#         # Clamp each dimension to data min/max with small margin
#         margin = 0.01  # 1% margin inside bounds
#         ranges = self.data_max - self.data_min
        
#         clamped_point = np.zeros_like(original_point)
#         for i in range(len(original_point)):
#             margin_val = margin * ranges[i]
#             clamped_point[i] = np.clip(original_point[i], 
#                                      self.data_min[i] + margin_val, 
#                                      self.data_max[i] - margin_val)
        
#         result = self._try_barycentric_interpolation(clamped_point)
#         if result is not None:
#             return result
        
#         # If still fails, try moving toward centroid from clamped position
#         for step_size in [0.1, 0.3, 0.5, 0.7]:
#             direction = self.data_centroid - clamped_point
#             adjusted_point = clamped_point + step_size * direction
#             result = self._try_barycentric_interpolation(adjusted_point)
#             if result is not None:
#                 return result
        
#         return None
    
#     def _nearest_neighbor_fallback(self, original_point):
#         # Normalize the query point using same scaling as training data
#         normalized_point = (original_point - self.data_min) / (self.data_max - self.data_min)
        
#         # Calculate distances using normalized inputs
#         distances = cdist([normalized_point], self.X_normalized, metric='euclidean')[0]
#         nearest_indices = np.argsort(distances)[:5]
#         nearest_distances = distances[nearest_indices]
        
#         # Inverse distance weighting
#         weights = 1.0 / (nearest_distances + 1e-10)
#         weights /= weights.sum()
        
#         # Weighted average using original (non-normalized) outputs
#         weighted_output = np.average(self.Y[nearest_indices], axis=0, weights=weights)
        
#         return weighted_output[0], weighted_output[1], weighted_output[2]

    
#     def interpolate(self, throttle, aileron_left, elevator, rudder, alpha, beta):
#         """
#         Main interpolation function - never fails, always returns a result
#         """
#         result, method = self.interpolate_with_fallbacks(throttle, aileron_left, elevator, rudder, alpha, beta)
        
#         # Log which method was used for debugging
#         if method != "barycentric":
#             print(f"Used fallback method: {method}")
        
#         return result


# def main():
#     if len(sys.argv) != 7:
#         print("Usage: python3 aero_interpolator.py <throttle> <aileron_left> <elevator> <rudder> <alpha> <beta>")
#         print("Example: python3 aero_interpolator.py 50.0 10.0 5.0 2.0 8.0 1.0")
#         sys.exit(1)
    
#     try:
#         # Parse command line arguments
#         throttle = float(sys.argv[1])
#         aileron_left = float(sys.argv[2])
#         elevator = float(sys.argv[3])
#         rudder = float(sys.argv[4])
#         alpha = float(sys.argv[5])
#         beta = float(sys.argv[6])
        
#         print(f"Input: throttle={throttle}, aileron_left={aileron_left}, elevator={elevator}, rudder={rudder}, alpha={alpha}, beta={beta}")
        
#         # Create interpolator
#         interpolator = AeroInterpolator("Aero_Data.csv")
        
#         # Perform never-fail interpolation
#         result = interpolator.interpolate(throttle, aileron_left, elevator, rudder, alpha, beta)
        
#         # Result is guaranteed to exist
#         output = {
#             "CL": float(result[0]),
#             "CD": float(result[1]), 
#             "CMm": float(result[2])
#         }
        
#         print("SUCCESS")
#         print(json.dumps(output))
        
#     except Exception as e:
#         # Even in case of catastrophic failure, return safe defaults
#         print(f"CRITICAL ERROR: {str(e)}")
#         print("Using emergency fallback values")
        
#         # Emergency physics-based defaults
#         try:
#             alpha_val = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
#             emergency_CL = 0.2 + 0.08 * alpha_val
#             emergency_CD = 0.05 + 0.5 * emergency_CL * emergency_CL
#             emergency_CMm = -0.05 * alpha_val
#         except:
#             emergency_CL, emergency_CD, emergency_CMm = 0.5, 0.08, 0.0
        
#         output = {
#             "CL": emergency_CL,
#             "CD": emergency_CD,
#             "CMm": emergency_CMm
#         }
        
#         print("SUCCESS")  # Still report success so C++ doesn't fail
#         print(json.dumps(output))


# if __name__ == "__main__":
#     main()
