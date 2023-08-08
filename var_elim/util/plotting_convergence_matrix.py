import matplotlib.pyplot as plt
import json
import numpy as np
path = r'/auto/nest/nest/u/sakshi/VarElem/variable-elimination/var_elim/data/convergence_analysis/gas_pipelines/'

with open(path+ "success_matrix_con_elim.txt") as f:
    data = f.read()

elim_data = json.loads(data)

with open(path+ "success_matrix.txt") as f:
    data = f.read()

full_data = {'43.0,63.0, ipopt_linsolve': 32.369, '43.0,63.0, nlp_eval_time': 0.671, '43.0,63.0, success': 1, '43.0,63.44444444444444, ipopt_linsolve': 31.398, '43.0,63.44444444444444, nlp_eval_time': 0.54, '43.0,63.44444444444444, success': 1, '43.0,63.888888888888886, ipopt_linsolve': 81.116, '43.0,63.888888888888886, nlp_eval_time': 1.062, '43.0,63.888888888888886, success': 1, '43.0,64.33333333333333, ipopt_linsolve': None, '43.0,64.33333333333333, nlp_eval_time': None, '43.0,64.33333333333333, success': 0, '43.0,64.77777777777777, ipopt_linsolve': 49.647, '43.0,64.77777777777777, nlp_eval_time': 0.86, '43.0,64.77777777777777, success': 1, '43.0,65.22222222222223, ipopt_linsolve': 24.112, '43.0,65.22222222222223, nlp_eval_time': 0.539, '43.0,65.22222222222223, success': 1, '43.0,65.66666666666667, ipopt_linsolve': None, '43.0,65.66666666666667, nlp_eval_time': None, '43.0,65.66666666666667, success': 0, '43.0,66.11111111111111, ipopt_linsolve': 31.673, '43.0,66.11111111111111, nlp_eval_time': 0.695, '43.0,66.11111111111111, success': 1, '43.0,66.55555555555556, ipopt_linsolve': None, '43.0,66.55555555555556, nlp_eval_time': None, '43.0,66.55555555555556, success': 0, '43.0,67.0, ipopt_linsolve': 45.583, '43.0,67.0, nlp_eval_time': 1.039, '43.0,67.0, success': 1, '43.44444444444444,63.0, ipopt_linsolve': None, '43.44444444444444,63.0, nlp_eval_time': None, '43.44444444444444,63.0, success': 0, '43.44444444444444,63.44444444444444, ipopt_linsolve': 45.124, '43.44444444444444,63.44444444444444, nlp_eval_time': 0.68, '43.44444444444444,63.44444444444444, success': 1, '43.44444444444444,63.888888888888886, ipopt_linsolve': 55.009, '43.44444444444444,63.888888888888886, nlp_eval_time': 0.776, '43.44444444444444,63.888888888888886, success': 1, '43.44444444444444,64.33333333333333, ipopt_linsolve': None, '43.44444444444444,64.33333333333333, nlp_eval_time': None, '43.44444444444444,64.33333333333333, success': 0, '43.44444444444444,64.77777777777777, ipopt_linsolve': 43.839, '43.44444444444444,64.77777777777777, nlp_eval_time': 0.823, '43.44444444444444,64.77777777777777, success': 1, '43.44444444444444,65.22222222222223, ipopt_linsolve': 37.231, '43.44444444444444,65.22222222222223, nlp_eval_time': 0.824, '43.44444444444444,65.22222222222223, success': 1, '43.44444444444444,65.66666666666667, ipopt_linsolve': 36.808, '43.44444444444444,65.66666666666667, nlp_eval_time': 0.793, '43.44444444444444,65.66666666666667, success': 1, '43.44444444444444,66.11111111111111, ipopt_linsolve': 31.993, '43.44444444444444,66.11111111111111, nlp_eval_time': 0.712, '43.44444444444444,66.11111111111111, success': 1, '43.44444444444444,66.55555555555556, ipopt_linsolve': 33.71, '43.44444444444444,66.55555555555556, nlp_eval_time': 0.781, '43.44444444444444,66.55555555555556, success': 1, '43.44444444444444,67.0, ipopt_linsolve': None, '43.44444444444444,67.0, nlp_eval_time': None, '43.44444444444444,67.0, success': 0, '43.888888888888886,63.0, ipopt_linsolve': None, '43.888888888888886,63.0, nlp_eval_time': None, '43.888888888888886,63.0, success': 0, '43.888888888888886,63.44444444444444, ipopt_linsolve': None, '43.888888888888886,63.44444444444444, nlp_eval_time': None, '43.888888888888886,63.44444444444444, success': 0, '43.888888888888886,63.888888888888886, ipopt_linsolve': 38.312, '43.888888888888886,63.888888888888886, nlp_eval_time': 0.618, '43.888888888888886,63.888888888888886, success': 1, '43.888888888888886,64.33333333333333, ipopt_linsolve': None, '43.888888888888886,64.33333333333333, nlp_eval_time': None, '43.888888888888886,64.33333333333333, success': 0, '43.888888888888886,64.77777777777777, ipopt_linsolve': 37.773, '43.888888888888886,64.77777777777777, nlp_eval_time': 0.628, '43.888888888888886,64.77777777777777, success': 1, '43.888888888888886,65.22222222222223, ipopt_linsolve': 45.654, '43.888888888888886,65.22222222222223, nlp_eval_time': 0.663, '43.888888888888886,65.22222222222223, success': 1, '43.888888888888886,65.66666666666667, ipopt_linsolve': 41.219, '43.888888888888886,65.66666666666667, nlp_eval_time': 0.754, '43.888888888888886,65.66666666666667, success': 1, '43.888888888888886,66.11111111111111, ipopt_linsolve': 37.168, '43.888888888888886,66.11111111111111, nlp_eval_time': 0.697, '43.888888888888886,66.11111111111111, success': 1, '43.888888888888886,66.55555555555556, ipopt_linsolve': 34.755, '43.888888888888886,66.55555555555556, nlp_eval_time': 0.731, '43.888888888888886,66.55555555555556, success': 1, '43.888888888888886,67.0, ipopt_linsolve': 28.565, '43.888888888888886,67.0, nlp_eval_time': 0.609, '43.888888888888886,67.0, success': 1, '44.333333333333336,63.0, ipopt_linsolve': 56.6, '44.333333333333336,63.0, nlp_eval_time': 0.827, '44.333333333333336,63.0, success': 1, '44.333333333333336,63.44444444444444, ipopt_linsolve': 41.556, '44.333333333333336,63.44444444444444, nlp_eval_time': 0.841, '44.333333333333336,63.44444444444444, success': 1, '44.333333333333336,63.888888888888886, ipopt_linsolve': 35.769, '44.333333333333336,63.888888888888886, nlp_eval_time': 0.83, '44.333333333333336,63.888888888888886, success': 1, '44.333333333333336,64.33333333333333, ipopt_linsolve': 33.715, '44.333333333333336,64.33333333333333, nlp_eval_time': 0.779, '44.333333333333336,64.33333333333333, success': 1, '44.333333333333336,64.77777777777777, ipopt_linsolve': 33.201, '44.333333333333336,64.77777777777777, nlp_eval_time': 0.691, '44.333333333333336,64.77777777777777, success': 1, '44.333333333333336,65.22222222222223, ipopt_linsolve': 33.193, '44.333333333333336,65.22222222222223, nlp_eval_time': 0.628, '44.333333333333336,65.22222222222223, success': 1, '44.333333333333336,65.66666666666667, ipopt_linsolve': 37.244, '44.333333333333336,65.66666666666667, nlp_eval_time': 0.596, '44.333333333333336,65.66666666666667, success': 1, '44.333333333333336,66.11111111111111, ipopt_linsolve': 31.444, '44.333333333333336,66.11111111111111, nlp_eval_time': 0.736, '44.333333333333336,66.11111111111111, success': 1, '44.333333333333336,66.55555555555556, ipopt_linsolve': None, '44.333333333333336,66.55555555555556, nlp_eval_time': None, '44.333333333333336,66.55555555555556, success': 0, '44.333333333333336,67.0, ipopt_linsolve': 60.532, '44.333333333333336,67.0, nlp_eval_time': 0.897, '44.333333333333336,67.0, success': 1, '44.77777777777778,63.0, ipopt_linsolve': None, '44.77777777777778,63.0, nlp_eval_time': None, '44.77777777777778,63.0, success': 0, '44.77777777777778,63.44444444444444, ipopt_linsolve': 46.984, '44.77777777777778,63.44444444444444, nlp_eval_time': 0.708, '44.77777777777778,63.44444444444444, success': 1, '44.77777777777778,63.888888888888886, ipopt_linsolve': 35.277, '44.77777777777778,63.888888888888886, nlp_eval_time': 0.812, '44.77777777777778,63.888888888888886, success': 1, '44.77777777777778,64.33333333333333, ipopt_linsolve': 27.012, '44.77777777777778,64.33333333333333, nlp_eval_time': 0.669, '44.77777777777778,64.33333333333333, success': 1, '44.77777777777778,64.77777777777777, ipopt_linsolve': 40.062, '44.77777777777778,64.77777777777777, nlp_eval_time': 0.851, '44.77777777777778,64.77777777777777, success': 1, '44.77777777777778,65.22222222222223, ipopt_linsolve': 34.8, '44.77777777777778,65.22222222222223, nlp_eval_time': 0.81, '44.77777777777778,65.22222222222223, success': 1, '44.77777777777778,65.66666666666667, ipopt_linsolve': 34.993, '44.77777777777778,65.66666666666667, nlp_eval_time': 0.779, '44.77777777777778,65.66666666666667, success': 1, '44.77777777777778,66.11111111111111, ipopt_linsolve': 27.574, '44.77777777777778,66.11111111111111, nlp_eval_time': 0.696, '44.77777777777778,66.11111111111111, success': 1, '44.77777777777778,66.55555555555556, ipopt_linsolve': None, '44.77777777777778,66.55555555555556, nlp_eval_time': None, '44.77777777777778,66.55555555555556, success': 0, '44.77777777777778,67.0, ipopt_linsolve': 41.681, '44.77777777777778,67.0, nlp_eval_time': 0.766, '44.77777777777778,67.0, success': 1, '45.22222222222222,63.0, ipopt_linsolve': 34.759, '45.22222222222222,63.0, nlp_eval_time': 0.59, '45.22222222222222,63.0, success': 1, '45.22222222222222,63.44444444444444, ipopt_linsolve': 43.359, '45.22222222222222,63.44444444444444, nlp_eval_time': 0.654, '45.22222222222222,63.44444444444444, success': 1, '45.22222222222222,63.888888888888886, ipopt_linsolve': 27.86, '45.22222222222222,63.888888888888886, nlp_eval_time': 0.594, '45.22222222222222,63.888888888888886, success': 1, '45.22222222222222,64.33333333333333, ipopt_linsolve': 35.848, '45.22222222222222,64.33333333333333, nlp_eval_time': 0.63, '45.22222222222222,64.33333333333333, success': 1, '45.22222222222222,64.77777777777777, ipopt_linsolve': 48.841, '45.22222222222222,64.77777777777777, nlp_eval_time': 0.87, '45.22222222222222,64.77777777777777, success': 1, '45.22222222222222,65.22222222222223, ipopt_linsolve': None, '45.22222222222222,65.22222222222223, nlp_eval_time': None, '45.22222222222222,65.22222222222223, success': 0, '45.22222222222222,65.66666666666667, ipopt_linsolve': 28.134, '45.22222222222222,65.66666666666667, nlp_eval_time': 0.669, '45.22222222222222,65.66666666666667, success': 1, '45.22222222222222,66.11111111111111, ipopt_linsolve': 37.889, '45.22222222222222,66.11111111111111, nlp_eval_time': 0.757, '45.22222222222222,66.11111111111111, success': 1, '45.22222222222222,66.55555555555556, ipopt_linsolve': None, '45.22222222222222,66.55555555555556, nlp_eval_time': None, '45.22222222222222,66.55555555555556, success': 0, '45.22222222222222,67.0, ipopt_linsolve': 41.077, '45.22222222222222,67.0, nlp_eval_time': 0.824, '45.22222222222222,67.0, success': 1, '45.666666666666664,63.0, ipopt_linsolve': 45.945, '45.666666666666664,63.0, nlp_eval_time': 0.833, '45.666666666666664,63.0, success': 1, '45.666666666666664,63.44444444444444, ipopt_linsolve': 7839.711, '45.666666666666664,63.44444444444444, nlp_eval_time': 70.62, '45.666666666666664,63.44444444444444, success': 1, '45.666666666666664,63.888888888888886, ipopt_linsolve': 32.105, '45.666666666666664,63.888888888888886, nlp_eval_time': 0.725, '45.666666666666664,63.888888888888886, success': 1, '45.666666666666664,64.33333333333333, ipopt_linsolve': 45.55, '45.666666666666664,64.33333333333333, nlp_eval_time': 0.678, '45.666666666666664,64.33333333333333, success': 1, '45.666666666666664,64.77777777777777, ipopt_linsolve': 37.423, '45.666666666666664,64.77777777777777, nlp_eval_time': 0.611, '45.666666666666664,64.77777777777777, success': 1, '45.666666666666664,65.22222222222223, ipopt_linsolve': None, '45.666666666666664,65.22222222222223, nlp_eval_time': None, '45.666666666666664,65.22222222222223, success': 0, '45.666666666666664,65.66666666666667, ipopt_linsolve': 42.733, '45.666666666666664,65.66666666666667, nlp_eval_time': 0.985, '45.666666666666664,65.66666666666667, success': 1, '45.666666666666664,66.11111111111111, ipopt_linsolve': None, '45.666666666666664,66.11111111111111, nlp_eval_time': None, '45.666666666666664,66.11111111111111, success': 0, '45.666666666666664,66.55555555555556, ipopt_linsolve': 28.426, '45.666666666666664,66.55555555555556, nlp_eval_time': 0.538, '45.666666666666664,66.55555555555556, success': 1, '45.666666666666664,67.0, ipopt_linsolve': 42.385, '45.666666666666664,67.0, nlp_eval_time': 0.68, '45.666666666666664,67.0, success': 1, '46.111111111111114,63.0, ipopt_linsolve': 38.186, '46.111111111111114,63.0, nlp_eval_time': 0.607, '46.111111111111114,63.0, success': 1, '46.111111111111114,63.44444444444444, ipopt_linsolve': None, '46.111111111111114,63.44444444444444, nlp_eval_time': None, '46.111111111111114,63.44444444444444, success': 0, '46.111111111111114,63.888888888888886, ipopt_linsolve': None, '46.111111111111114,63.888888888888886, nlp_eval_time': None, '46.111111111111114,63.888888888888886, success': 0, '46.111111111111114,64.33333333333333, ipopt_linsolve': 34.641, '46.111111111111114,64.33333333333333, nlp_eval_time': 0.795, '46.111111111111114,64.33333333333333, success': 1, '46.111111111111114,64.77777777777777, ipopt_linsolve': 33.27, '46.111111111111114,64.77777777777777, nlp_eval_time': 0.71, '46.111111111111114,64.77777777777777, success': 1, '46.111111111111114,65.22222222222223, ipopt_linsolve': 49.191, '46.111111111111114,65.22222222222223, nlp_eval_time': 0.903, '46.111111111111114,65.22222222222223, success': 1, '46.111111111111114,65.66666666666667, ipopt_linsolve': None, '46.111111111111114,65.66666666666667, nlp_eval_time': None, '46.111111111111114,65.66666666666667, success': 0, '46.111111111111114,66.11111111111111, ipopt_linsolve': 29.221, '46.111111111111114,66.11111111111111, nlp_eval_time': 0.661, '46.111111111111114,66.11111111111111, success': 1, '46.111111111111114,66.55555555555556, ipopt_linsolve': 27.138, '46.111111111111114,66.55555555555556, nlp_eval_time': 0.704, '46.111111111111114,66.55555555555556, success': 1, '46.111111111111114,67.0, ipopt_linsolve': 74.536, '46.111111111111114,67.0, nlp_eval_time': 0.98, '46.111111111111114,67.0, success': 1, '46.55555555555556,63.0, ipopt_linsolve': 43.494, '46.55555555555556,63.0, nlp_eval_time': 0.668, '46.55555555555556,63.0, success': 1, '46.55555555555556,63.44444444444444, ipopt_linsolve': 34.599, '46.55555555555556,63.44444444444444, nlp_eval_time': 0.617, '46.55555555555556,63.44444444444444, success': 1, '46.55555555555556,63.888888888888886, ipopt_linsolve': None, '46.55555555555556,63.888888888888886, nlp_eval_time': None, '46.55555555555556,63.888888888888886, success': 0, '46.55555555555556,64.33333333333333, ipopt_linsolve': 34.69, '46.55555555555556,64.33333333333333, nlp_eval_time': 0.878, '46.55555555555556,64.33333333333333, success': 1, '46.55555555555556,64.77777777777777, ipopt_linsolve': 29.189, '46.55555555555556,64.77777777777777, nlp_eval_time': 0.566, '46.55555555555556,64.77777777777777, success': 1, '46.55555555555556,65.22222222222223, ipopt_linsolve': 41.601, '46.55555555555556,65.22222222222223, nlp_eval_time': 0.647, '46.55555555555556,65.22222222222223, success': 1, '46.55555555555556,65.66666666666667, ipopt_linsolve': 78.692, '46.55555555555556,65.66666666666667, nlp_eval_time': 0.898, '46.55555555555556,65.66666666666667, success': 1, '46.55555555555556,66.11111111111111, ipopt_linsolve': 29.252, '46.55555555555556,66.11111111111111, nlp_eval_time': 0.622, '46.55555555555556,66.11111111111111, success': 1, '46.55555555555556,66.55555555555556, ipopt_linsolve': 62.574, '46.55555555555556,66.55555555555556, nlp_eval_time': 0.844, '46.55555555555556,66.55555555555556, success': 1, '46.55555555555556,67.0, ipopt_linsolve': 51.159, '46.55555555555556,67.0, nlp_eval_time': 0.644, '46.55555555555556,67.0, success': 1, '47.0,63.0, ipopt_linsolve': 28.529, '47.0,63.0, nlp_eval_time': 0.648, '47.0,63.0, success': 1, '47.0,63.44444444444444, ipopt_linsolve': None, '47.0,63.44444444444444, nlp_eval_time': None, '47.0,63.44444444444444, success': 0, '47.0,63.888888888888886, ipopt_linsolve': None, '47.0,63.888888888888886, nlp_eval_time': None, '47.0,63.888888888888886, success': 0, '47.0,64.33333333333333, ipopt_linsolve': 42.872, '47.0,64.33333333333333, nlp_eval_time': 0.833, '47.0,64.33333333333333, success': 1, '47.0,64.77777777777777, ipopt_linsolve': None, '47.0,64.77777777777777, nlp_eval_time': None, '47.0,64.77777777777777, success': 0, '47.0,65.22222222222223, ipopt_linsolve': 48.35, '47.0,65.22222222222223, nlp_eval_time': 1.092, '47.0,65.22222222222223, success': 1, '47.0,65.66666666666667, ipopt_linsolve': 32.881, '47.0,65.66666666666667, nlp_eval_time': 0.714, '47.0,65.66666666666667, success': 1, '47.0,66.11111111111111, ipopt_linsolve': 28.89, '47.0,66.11111111111111, nlp_eval_time': 0.623, '47.0,66.11111111111111, success': 1, '47.0,66.55555555555556, ipopt_linsolve': 55.939, '47.0,66.55555555555556, nlp_eval_time': 0.786, '47.0,66.55555555555556, success': 1, '47.0,67.0, ipopt_linsolve': 30.55, '47.0,67.0, nlp_eval_time': 0.798, '47.0,67.0, success': 1}
elim_data = {'43.0,63.0, ipopt_linsolve': 4.353, '43.0,63.0, nlp_eval_time': 0.094, '43.0,63.0, success': 1, '43.0,63.44444444444444, ipopt_linsolve': 4.379, '43.0,63.44444444444444, nlp_eval_time': 0.081, '43.0,63.44444444444444, success': 1, '43.0,63.888888888888886, ipopt_linsolve': 5.269, '43.0,63.888888888888886, nlp_eval_time': 0.096, '43.0,63.888888888888886, success': 1, '43.0,64.33333333333333, ipopt_linsolve': 5.146, '43.0,64.33333333333333, nlp_eval_time': 0.095, '43.0,64.33333333333333, success': 1, '43.0,64.77777777777777, ipopt_linsolve': 4.953, '43.0,64.77777777777777, nlp_eval_time': 0.088, '43.0,64.77777777777777, success': 1, '43.0,65.22222222222223, ipopt_linsolve': 4.458, '43.0,65.22222222222223, nlp_eval_time': 0.07, '43.0,65.22222222222223, success': 1, '43.0,65.66666666666667, ipopt_linsolve': 7.987, '43.0,65.66666666666667, nlp_eval_time': 0.133, '43.0,65.66666666666667, success': 1, '43.0,66.11111111111111, ipopt_linsolve': 4.804, '43.0,66.11111111111111, nlp_eval_time': 0.089, '43.0,66.11111111111111, success': 1, '43.0,66.55555555555556, ipopt_linsolve': 4.527, '43.0,66.55555555555556, nlp_eval_time': 0.09, '43.0,66.55555555555556, success': 1, '43.0,67.0, ipopt_linsolve': 5.663, '43.0,67.0, nlp_eval_time': 0.09, '43.0,67.0, success': 1, '43.44444444444444,63.0, ipopt_linsolve': 4.002, '43.44444444444444,63.0, nlp_eval_time': 0.066, '43.44444444444444,63.0, success': 1, '43.44444444444444,63.44444444444444, ipopt_linsolve': 4.52, '43.44444444444444,63.44444444444444, nlp_eval_time': 0.088, '43.44444444444444,63.44444444444444, success': 1, '43.44444444444444,63.888888888888886, ipopt_linsolve': 5.331, '43.44444444444444,63.888888888888886, nlp_eval_time': 0.092, '43.44444444444444,63.888888888888886, success': 1, '43.44444444444444,64.33333333333333, ipopt_linsolve': 5.879, '43.44444444444444,64.33333333333333, nlp_eval_time': 0.104, '43.44444444444444,64.33333333333333, success': 1, '43.44444444444444,64.77777777777777, ipopt_linsolve': 4.661, '43.44444444444444,64.77777777777777, nlp_eval_time': 0.082, '43.44444444444444,64.77777777777777, success': 1, '43.44444444444444,65.22222222222223, ipopt_linsolve': 6.289, '43.44444444444444,65.22222222222223, nlp_eval_time': 0.075, '43.44444444444444,65.22222222222223, success': 1, '43.44444444444444,65.66666666666667, ipopt_linsolve': 4.553, '43.44444444444444,65.66666666666667, nlp_eval_time': 0.087, '43.44444444444444,65.66666666666667, success': 1, '43.44444444444444,66.11111111111111, ipopt_linsolve': 6.415, '43.44444444444444,66.11111111111111, nlp_eval_time': 0.096, '43.44444444444444,66.11111111111111, success': 1, '43.44444444444444,66.55555555555556, ipopt_linsolve': 7.031, '43.44444444444444,66.55555555555556, nlp_eval_time': 0.11, '43.44444444444444,66.55555555555556, success': 1, '43.44444444444444,67.0, ipopt_linsolve': 4.858, '43.44444444444444,67.0, nlp_eval_time': 0.087, '43.44444444444444,67.0, success': 1, '43.888888888888886,63.0, ipopt_linsolve': 5.374, '43.888888888888886,63.0, nlp_eval_time': 0.079, '43.888888888888886,63.0, success': 1, '43.888888888888886,63.44444444444444, ipopt_linsolve': 4.389, '43.888888888888886,63.44444444444444, nlp_eval_time': 0.084, '43.888888888888886,63.44444444444444, success': 1, '43.888888888888886,63.888888888888886, ipopt_linsolve': 4.923, '43.888888888888886,63.888888888888886, nlp_eval_time': 0.097, '43.888888888888886,63.888888888888886, success': 1, '43.888888888888886,64.33333333333333, ipopt_linsolve': 4.675, '43.888888888888886,64.33333333333333, nlp_eval_time': 0.076, '43.888888888888886,64.33333333333333, success': 1, '43.888888888888886,64.77777777777777, ipopt_linsolve': 4.258, '43.888888888888886,64.77777777777777, nlp_eval_time': 0.069, '43.888888888888886,64.77777777777777, success': 1, '43.888888888888886,65.22222222222223, ipopt_linsolve': 5.389, '43.888888888888886,65.22222222222223, nlp_eval_time': 0.124, '43.888888888888886,65.22222222222223, success': 1, '43.888888888888886,65.66666666666667, ipopt_linsolve': 4.603, '43.888888888888886,65.66666666666667, nlp_eval_time': 0.076, '43.888888888888886,65.66666666666667, success': 1, '43.888888888888886,66.11111111111111, ipopt_linsolve': 5.179, '43.888888888888886,66.11111111111111, nlp_eval_time': 0.09, '43.888888888888886,66.11111111111111, success': 1, '43.888888888888886,66.55555555555556, ipopt_linsolve': 6.549, '43.888888888888886,66.55555555555556, nlp_eval_time': 0.101, '43.888888888888886,66.55555555555556, success': 1, '43.888888888888886,67.0, ipopt_linsolve': 5.631, '43.888888888888886,67.0, nlp_eval_time': 0.105, '43.888888888888886,67.0, success': 1, '44.333333333333336,63.0, ipopt_linsolve': 4.84, '44.333333333333336,63.0, nlp_eval_time': 0.086, '44.333333333333336,63.0, success': 1, '44.333333333333336,63.44444444444444, ipopt_linsolve': 5.118, '44.333333333333336,63.44444444444444, nlp_eval_time': 0.148, '44.333333333333336,63.44444444444444, success': 1, '44.333333333333336,63.888888888888886, ipopt_linsolve': 4.617, '44.333333333333336,63.888888888888886, nlp_eval_time': 0.101, '44.333333333333336,63.888888888888886, success': 1, '44.333333333333336,64.33333333333333, ipopt_linsolve': 4.557, '44.333333333333336,64.33333333333333, nlp_eval_time': 0.087, '44.333333333333336,64.33333333333333, success': 1, '44.333333333333336,64.77777777777777, ipopt_linsolve': 4.709, '44.333333333333336,64.77777777777777, nlp_eval_time': 0.09, '44.333333333333336,64.77777777777777, success': 1, '44.333333333333336,65.22222222222223, ipopt_linsolve': 4.575, '44.333333333333336,65.22222222222223, nlp_eval_time': 0.088, '44.333333333333336,65.22222222222223, success': 1, '44.333333333333336,65.66666666666667, ipopt_linsolve': 6.1, '44.333333333333336,65.66666666666667, nlp_eval_time': 0.099, '44.333333333333336,65.66666666666667, success': 1, '44.333333333333336,66.11111111111111, ipopt_linsolve': 4.651, '44.333333333333336,66.11111111111111, nlp_eval_time': 0.089, '44.333333333333336,66.11111111111111, success': 1, '44.333333333333336,66.55555555555556, ipopt_linsolve': 5.018, '44.333333333333336,66.55555555555556, nlp_eval_time': 0.094, '44.333333333333336,66.55555555555556, success': 1, '44.333333333333336,67.0, ipopt_linsolve': 4.801, '44.333333333333336,67.0, nlp_eval_time': 0.09, '44.333333333333336,67.0, success': 1, '44.77777777777778,63.0, ipopt_linsolve': 4.503, '44.77777777777778,63.0, nlp_eval_time': 0.092, '44.77777777777778,63.0, success': 1, '44.77777777777778,63.44444444444444, ipopt_linsolve': 5.172, '44.77777777777778,63.44444444444444, nlp_eval_time': 0.086, '44.77777777777778,63.44444444444444, success': 1, '44.77777777777778,63.888888888888886, ipopt_linsolve': 5.568, '44.77777777777778,63.888888888888886, nlp_eval_time': 0.124, '44.77777777777778,63.888888888888886, success': 1, '44.77777777777778,64.33333333333333, ipopt_linsolve': 5.391, '44.77777777777778,64.33333333333333, nlp_eval_time': 0.093, '44.77777777777778,64.33333333333333, success': 1, '44.77777777777778,64.77777777777777, ipopt_linsolve': 4.855, '44.77777777777778,64.77777777777777, nlp_eval_time': 0.113, '44.77777777777778,64.77777777777777, success': 1, '44.77777777777778,65.22222222222223, ipopt_linsolve': 5.404, '44.77777777777778,65.22222222222223, nlp_eval_time': 0.101, '44.77777777777778,65.22222222222223, success': 1, '44.77777777777778,65.66666666666667, ipopt_linsolve': 4.626, '44.77777777777778,65.66666666666667, nlp_eval_time': 0.089, '44.77777777777778,65.66666666666667, success': 1, '44.77777777777778,66.11111111111111, ipopt_linsolve': 4.626, '44.77777777777778,66.11111111111111, nlp_eval_time': 0.086, '44.77777777777778,66.11111111111111, success': 1, '44.77777777777778,66.55555555555556, ipopt_linsolve': 5.032, '44.77777777777778,66.55555555555556, nlp_eval_time': 0.106, '44.77777777777778,66.55555555555556, success': 1, '44.77777777777778,67.0, ipopt_linsolve': 5.223, '44.77777777777778,67.0, nlp_eval_time': 0.103, '44.77777777777778,67.0, success': 1, '45.22222222222222,63.0, ipopt_linsolve': 4.856, '45.22222222222222,63.0, nlp_eval_time': 0.097, '45.22222222222222,63.0, success': 1, '45.22222222222222,63.44444444444444, ipopt_linsolve': 4.809, '45.22222222222222,63.44444444444444, nlp_eval_time': 0.086, '45.22222222222222,63.44444444444444, success': 1, '45.22222222222222,63.888888888888886, ipopt_linsolve': 5.472, '45.22222222222222,63.888888888888886, nlp_eval_time': 0.108, '45.22222222222222,63.888888888888886, success': 1, '45.22222222222222,64.33333333333333, ipopt_linsolve': 5.273, '45.22222222222222,64.33333333333333, nlp_eval_time': 0.106, '45.22222222222222,64.33333333333333, success': 1, '45.22222222222222,64.77777777777777, ipopt_linsolve': 4.856, '45.22222222222222,64.77777777777777, nlp_eval_time': 0.088, '45.22222222222222,64.77777777777777, success': 1, '45.22222222222222,65.22222222222223, ipopt_linsolve': 5.296, '45.22222222222222,65.22222222222223, nlp_eval_time': 0.089, '45.22222222222222,65.22222222222223, success': 1, '45.22222222222222,65.66666666666667, ipopt_linsolve': 4.702, '45.22222222222222,65.66666666666667, nlp_eval_time': 0.083, '45.22222222222222,65.66666666666667, success': 1, '45.22222222222222,66.11111111111111, ipopt_linsolve': 4.503, '45.22222222222222,66.11111111111111, nlp_eval_time': 0.085, '45.22222222222222,66.11111111111111, success': 1, '45.22222222222222,66.55555555555556, ipopt_linsolve': 5.216, '45.22222222222222,66.55555555555556, nlp_eval_time': 0.094, '45.22222222222222,66.55555555555556, success': 1, '45.22222222222222,67.0, ipopt_linsolve': 5.704, '45.22222222222222,67.0, nlp_eval_time': 0.101, '45.22222222222222,67.0, success': 1, '45.666666666666664,63.0, ipopt_linsolve': 4.364, '45.666666666666664,63.0, nlp_eval_time': 0.083, '45.666666666666664,63.0, success': 1, '45.666666666666664,63.44444444444444, ipopt_linsolve': 4.621, '45.666666666666664,63.44444444444444, nlp_eval_time': 0.079, '45.666666666666664,63.44444444444444, success': 1, '45.666666666666664,63.888888888888886, ipopt_linsolve': 4.611, '45.666666666666664,63.888888888888886, nlp_eval_time': 0.08, '45.666666666666664,63.888888888888886, success': 1, '45.666666666666664,64.33333333333333, ipopt_linsolve': 4.789, '45.666666666666664,64.33333333333333, nlp_eval_time': 0.089, '45.666666666666664,64.33333333333333, success': 1, '45.666666666666664,64.77777777777777, ipopt_linsolve': 4.931, '45.666666666666664,64.77777777777777, nlp_eval_time': 0.09, '45.666666666666664,64.77777777777777, success': 1, '45.666666666666664,65.22222222222223, ipopt_linsolve': 4.631, '45.666666666666664,65.22222222222223, nlp_eval_time': 0.099, '45.666666666666664,65.22222222222223, success': 1, '45.666666666666664,65.66666666666667, ipopt_linsolve': 5.925, '45.666666666666664,65.66666666666667, nlp_eval_time': 0.111, '45.666666666666664,65.66666666666667, success': 1, '45.666666666666664,66.11111111111111, ipopt_linsolve': 5.347, '45.666666666666664,66.11111111111111, nlp_eval_time': 0.089, '45.666666666666664,66.11111111111111, success': 1, '45.666666666666664,66.55555555555556, ipopt_linsolve': 5.23, '45.666666666666664,66.55555555555556, nlp_eval_time': 0.091, '45.666666666666664,66.55555555555556, success': 1, '45.666666666666664,67.0, ipopt_linsolve': 4.447, '45.666666666666664,67.0, nlp_eval_time': 0.085, '45.666666666666664,67.0, success': 1, '46.111111111111114,63.0, ipopt_linsolve': 5.631, '46.111111111111114,63.0, nlp_eval_time': 0.089, '46.111111111111114,63.0, success': 1, '46.111111111111114,63.44444444444444, ipopt_linsolve': 4.693, '46.111111111111114,63.44444444444444, nlp_eval_time': 0.073, '46.111111111111114,63.44444444444444, success': 1, '46.111111111111114,63.888888888888886, ipopt_linsolve': 4.583, '46.111111111111114,63.888888888888886, nlp_eval_time': 0.091, '46.111111111111114,63.888888888888886, success': 1, '46.111111111111114,64.33333333333333, ipopt_linsolve': 4.97, '46.111111111111114,64.33333333333333, nlp_eval_time': 0.084, '46.111111111111114,64.33333333333333, success': 1, '46.111111111111114,64.77777777777777, ipopt_linsolve': 6.029, '46.111111111111114,64.77777777777777, nlp_eval_time': 0.08, '46.111111111111114,64.77777777777777, success': 1, '46.111111111111114,65.22222222222223, ipopt_linsolve': 4.699, '46.111111111111114,65.22222222222223, nlp_eval_time': 0.084, '46.111111111111114,65.22222222222223, success': 1, '46.111111111111114,65.66666666666667, ipopt_linsolve': 4.517, '46.111111111111114,65.66666666666667, nlp_eval_time': 0.076, '46.111111111111114,65.66666666666667, success': 1, '46.111111111111114,66.11111111111111, ipopt_linsolve': 4.776, '46.111111111111114,66.11111111111111, nlp_eval_time': 0.115, '46.111111111111114,66.11111111111111, success': 1, '46.111111111111114,66.55555555555556, ipopt_linsolve': 5.28, '46.111111111111114,66.55555555555556, nlp_eval_time': 0.08, '46.111111111111114,66.55555555555556, success': 1, '46.111111111111114,67.0, ipopt_linsolve': 5.803, '46.111111111111114,67.0, nlp_eval_time': 0.085, '46.111111111111114,67.0, success': 1, '46.55555555555556,63.0, ipopt_linsolve': 5.406, '46.55555555555556,63.0, nlp_eval_time': 0.106, '46.55555555555556,63.0, success': 1, '46.55555555555556,63.44444444444444, ipopt_linsolve': 5.269, '46.55555555555556,63.44444444444444, nlp_eval_time': 0.113, '46.55555555555556,63.44444444444444, success': 1, '46.55555555555556,63.888888888888886, ipopt_linsolve': 5.059, '46.55555555555556,63.888888888888886, nlp_eval_time': 0.076, '46.55555555555556,63.888888888888886, success': 1, '46.55555555555556,64.33333333333333, ipopt_linsolve': 4.584, '46.55555555555556,64.33333333333333, nlp_eval_time': 0.093, '46.55555555555556,64.33333333333333, success': 1, '46.55555555555556,64.77777777777777, ipopt_linsolve': 5.134, '46.55555555555556,64.77777777777777, nlp_eval_time': 0.094, '46.55555555555556,64.77777777777777, success': 1, '46.55555555555556,65.22222222222223, ipopt_linsolve': 4.568, '46.55555555555556,65.22222222222223, nlp_eval_time': 0.072, '46.55555555555556,65.22222222222223, success': 1, '46.55555555555556,65.66666666666667, ipopt_linsolve': 4.665, '46.55555555555556,65.66666666666667, nlp_eval_time': 0.084, '46.55555555555556,65.66666666666667, success': 1, '46.55555555555556,66.11111111111111, ipopt_linsolve': 4.495, '46.55555555555556,66.11111111111111, nlp_eval_time': 0.072, '46.55555555555556,66.11111111111111, success': 1, '46.55555555555556,66.55555555555556, ipopt_linsolve': 4.915, '46.55555555555556,66.55555555555556, nlp_eval_time': 0.121, '46.55555555555556,66.55555555555556, success': 1, '46.55555555555556,67.0, ipopt_linsolve': 6.354, '46.55555555555556,67.0, nlp_eval_time': 0.112, '46.55555555555556,67.0, success': 1, '47.0,63.0, ipopt_linsolve': 4.627, '47.0,63.0, nlp_eval_time': 0.098, '47.0,63.0, success': 1, '47.0,63.44444444444444, ipopt_linsolve': 4.526, '47.0,63.44444444444444, nlp_eval_time': 0.09, '47.0,63.44444444444444, success': 1, '47.0,63.888888888888886, ipopt_linsolve': 4.807, '47.0,63.888888888888886, nlp_eval_time': 0.118, '47.0,63.888888888888886, success': 1, '47.0,64.33333333333333, ipopt_linsolve': 4.975, '47.0,64.33333333333333, nlp_eval_time': 0.091, '47.0,64.33333333333333, success': 1, '47.0,64.77777777777777, ipopt_linsolve': 5.468, '47.0,64.77777777777777, nlp_eval_time': 0.118, '47.0,64.77777777777777, success': 1, '47.0,65.22222222222223, ipopt_linsolve': 5.866, '47.0,65.22222222222223, nlp_eval_time': 0.11, '47.0,65.22222222222223, success': 1, '47.0,65.66666666666667, ipopt_linsolve': 5.065, '47.0,65.66666666666667, nlp_eval_time': 0.101, '47.0,65.66666666666667, success': 1, '47.0,66.11111111111111, ipopt_linsolve': 5.235, '47.0,66.11111111111111, nlp_eval_time': 0.097, '47.0,66.11111111111111, success': 1, '47.0,66.55555555555556, ipopt_linsolve': 6.451, '47.0,66.55555555555556, nlp_eval_time': 0.158, '47.0,66.55555555555556, success': 1, '47.0,67.0, ipopt_linsolve': 4.942, '47.0,67.0, nlp_eval_time': 0.097, '47.0,67.0, success': 1}

full_data_list = []
elim_data_list = []
for key in full_data.keys():
    if key.endswith('success'):
        full_data_list.append(full_data[key])
        elim_data_list.append(elim_data[key])

full_data_array = np.array(full_data_list)
elim_data_array = np.array(elim_data_list)

full_data_reshaped = np.reshape(full_data_array, (10,10))
elim_data_reshaped = np.reshape(elim_data_array, (10,10))

demand_mid = np.linspace(43, 47, 10)
demand_end = np.linspace(63, 67, 10)

x_label = []
y_label = []
for d in demand_mid:
    y_label.append(str(round(d, 1)))
    
for p in demand_end:
    x_label.append(str(round(p, 1)))


plt.figure()
plt.spy(full_data_reshaped, markersize=15,color = 'coral')
plt.spy(full_data_reshaped == 0, markersize = 15,color = 'k')
plt.xticks(range(0, 10), x_label, rotation = 90)
plt.yticks(range(0,10), y_label)
plt.ylabel('demand parameter 1')
plt.xlabel('demand parameter 2')
plt.title('Convergence for full model')
plt.savefig("gas_pipelines_full_100_100.png", dpi = 300, bbox_inches = "tight")

plt.figure()
plt.spy(elim_data_reshaped, markersize=15, color = 'coral')
plt.spy(elim_data_reshaped == 0, markersize = 15, color ='k')
plt.xticks(range(0, 10), x_label, rotation = 90,)
plt.yticks(range(0,10), y_label)
plt.ylabel('demand parameter 1')
plt.xlabel('demand parameter 2')
plt.title('Convergence for reduced model')
plt.savefig("gas_pipelines_reduced_100_100.png", dpi = 300, bbox_inches = "tight")

#plt.xlabel(x_label)
