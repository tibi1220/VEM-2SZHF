from .data import get_data
from .matrix import ext_matrix, sub_matrix, ext_vector, sub_vector

import sympy as sp
import numpy as np


def get_Jacobian(coords, indices, diff_N):
    J_11 = 0
    J_12 = 0
    J_21 = 0
    J_22 = 0
    for i in range(4):
        J_11 += diff_N[i]["xi"] * coords[indices[i]]['x']
        J_12 += diff_N[i]["xi"] * coords[indices[i]]['y']
        J_21 += diff_N[i]['eta'] * coords[indices[i]]['x']
        J_22 += diff_N[i]["eta"] * coords[indices[i]]['y']
    J = sp.Matrix([[J_11, J_12], [J_21, J_22]])
    return J


def calculate(config):
    Data = get_data(config["code"])
    Data["name"] = config["name"]
    Data["neptun"] = config["neptun"]

    a = Data['a']
    b = Data['b']
    c = Data['c']
    t = Data['t']
    state = Data['state']
    p = Data['p']
    E = Data['E']
    nu = Data['nu']

    code = Data['code']

    parametric = Data['parametric']
    numeric = Data['numeric'](code[1], a, b, c)

    # Define symbols
    xi, eta = sp.symbols("ξ,η")

    # Shape functions
    N = [
        1/4*(1-xi)*(1-eta),
        1/4*(1+xi)*(1-eta),
        1/4*(1+xi)*(1+eta),
        1/4*(1-xi)*(1+eta)
    ]

    # Calculate connection between local and global coord systems
    x, y = sp.symbols('x, y')
    transform = []

    for i in range(3):
        # Left hand side of eqs
        l_1, l_2 = 0, 0
        for j in range(4):
            # N_i * x_i
            l_1 += N[j] * numeric[parametric['rectangles'][i][j]]["x"]
            # N_i * y_i
            l_2 += N[j] * numeric[parametric['rectangles'][i][j]]["y"]

        l_1 = sp.simplify(l_1)
        l_2 = sp.simplify(l_2)

        eq_1 = sp.Eq(l_1, x)
        eq_2 = sp.Eq(l_2, y)

        sol = sp.solve([eq_1, eq_2], (xi, eta))

        # Right hand side of eqs
        r_1, r_2 = None, None

        if type(sol) is dict:
            r_1 = sol[xi]
            r_2 = sol[eta]
        else:
            [(r_1, r_2)] = sp.solve((eq_1, eq_2), (xi, eta))

        transform.append({
            "xi": sp.simplify(r_1),
            "eta": sp.simplify(r_2),
            "x": l_1,
            "y": l_2,
        })

    # Shape function derivatives in local coordinate system
    diff_local_N = []

    for i in range(4):
        diff_local_N.append({
            "xi": sp.diff(N[i], xi),
            "eta": sp.diff(N[i], eta)
        })

    # Jacobian and inverse Jacibian matrices
    J = []
    inv_J = []

    for i in range(3):
        J.append(get_Jacobian(
            numeric, parametric["rectangles"][i], diff_local_N)
        )
        inv_J.append(J[i].inv())

    # Shape function derivatives in global coordinate system
    diff_global_N = [[], [], []]

    for i in range(3):
        for j in range(4):
            diff_global_N[i].append({
                "x": inv_J[i][0] * diff_local_N[j]["xi"]
                + inv_J[i][1] * diff_local_N[j]["eta"],
                "y": inv_J[i][2] * diff_local_N[j]["xi"]
                + inv_J[i][3] * diff_local_N[j]["eta"],
            })

    # Calculate B matrix
    B = []
    for i in range(3):
        d = diff_global_N[i]
        B.append(sp.simplify(sp.Matrix([
            [d[0]['x'], 0, d[1]['x'], 0, d[2]['x'], 0, d[3]['x'], 0],
            [0, d[0]['y'], 0, d[1]['y'], 0, d[2]['y'], 0, d[3]['y']],
            [d[0]['y'], d[0]['x'], d[1]['y'], d[1]['x'],
                d[2]['y'], d[2]['x'], d[3]['y'], d[3]['x']]
        ])))

    # Calculate D matrix
    if state == "SF":
        D = E / (1 - nu**2) * sp.Matrix([
            [1,  nu, 0],
            [nu, 1,  0],
            [0,  0, (1-nu)/2]
        ])
    else:
        D = E / (1 + nu) / (1 - 2*nu) * sp.Matrix([
            [1-nu, nu,   0],
            [nu,   1-nu, 0],
            [0,    0,    (1 - 2*nu) / 2]
        ])

    # Calculate K_i matrices
    K_i = []

    for i in range(3):
        K_sym = sp.Transpose(B[i]) * D * B[i] * J[i].det() * t
        K_i.append(
            K_sym.subs({xi: -1/np.sqrt(3), eta: -1/np.sqrt(3)}) +
            K_sym.subs({xi: +1/np.sqrt(3), eta: -1/np.sqrt(3)}) +
            K_sym.subs({xi: +1/np.sqrt(3), eta: +1/np.sqrt(3)}) +
            K_sym.subs({xi: -1/np.sqrt(3), eta: +1/np.sqrt(3)})
        )

    # Calculate DOF matrix
    r = sp.Matrix(parametric["rectangles"])

    DOF = [
        [2*r[0, 0]-1, 2*r[0, 0], 2*r[0, 1]-1, 2*r[0, 1],
            2*r[0, 2]-1, 2*r[0, 2], 2*r[0, 3]-1, 2*r[0, 3]],
        [2*r[1, 0]-1, 2*r[1, 0], 2*r[1, 1]-1, 2*r[1, 1],
            2*r[1, 2]-1, 2*r[1, 2], 2*r[1, 3]-1, 2*r[1, 3]],
        [2*r[2, 0]-1, 2*r[2, 0], 2*r[2, 1]-1, 2*r[2, 1],
            2*r[2, 2]-1, 2*r[2, 2], 2*r[2, 3]-1, 2*r[2, 3]]
    ]

    # Calculate global K matrix
    K = sp.Matrix(
        ext_matrix(K_i[0], DOF[0], 16) +
        ext_matrix(K_i[1], DOF[1], 16) +
        ext_matrix(K_i[2], DOF[2], 16)
    )

    # Create symbolic U and F
    U = sp.MatrixSymbol("U", 16, 1)
    F = sp.MatrixSymbol("F", 16, 1)

    u_1, u_2, u_3, u_4, u_5, u_6, u_7, u_8 = sp.symbols(
        "U_1, U_2, U_3, U_4, U_5, U_6, U_7, U_8"
    )
    v_1, v_2, v_3, v_4, v_5, v_6, v_7, v_8 = sp.symbols(
        "V_1, V_2, V_3, V_4, V_5, V_6, V_7, V_8"
    )
    X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8 = sp.symbols(
        "F_1_x, F_2_x, F_3_x, F_4_x, F_5_x, F_6_x, F_7_x, F_8_x"
    )
    Y_1, Y_2, Y_3, Y_4, Y_5, Y_6, Y_7, Y_8 = sp.symbols(
        "F_1_y, F_2_y, F_3_y, F_4_y, F_5_y, F_6_y, F_7_y, F_8_y"
    )

    U_sym = sp.Matrix([[u_1], [v_1], [u_2], [v_2], [u_3], [v_3], [u_4], [v_4],
                       [u_5], [v_5], [u_6], [v_6], [u_7], [v_7], [u_8], [v_8]])
    F_sym = sp.Matrix([[X_1], [Y_1], [X_2], [Y_2], [X_3], [Y_3], [X_4], [Y_4],
                       [X_5], [Y_5], [X_6], [Y_6], [X_7], [Y_7], [X_8], [Y_8]])

    # Solva KU = F
    free = parametric['free']
    not_free = parametric["not_free"]
    distributed = parametric["distributed"]

    p_s, a_s, t_s = sp.symbols("p, a, t")

    for d in not_free:
        U_sym[d - 1] = 0

    for d in free:
        F_sym[d - 1] = 0

    for d in range(2):
        F_sym[distributed[d] - 1] = distributed[2] * p_s / 2 * a_s * t_s

    tmp = sp.Matrix([free])
    K_kond = sub_matrix(K, tmp)
    F_kond = sub_vector(F_sym.subs({p_s: p, a_s: a, t_s: t}), tmp)
    U_sym_kond = sub_vector(U_sym, tmp)

    K_kond_inv = K_kond.inv()
    U_kond = K_kond_inv * F_kond

    U_calc = ext_vector(U_kond, tmp, 16)
    F_calc = K * U_calc

    F_base = ext_vector(F_kond, tmp, 16)
    F_reac = F_calc - F_base

    # Calculate Delta in mm and in um
    Delta_mm = []

    for i in range(8):
        Delta_mm.append([(U_calc[2*i, 0]**2 + U_calc[2*i+1, 0]**2)**(1/2)])

    Delta_mm = sp.Matrix(Delta_mm)

    # Calculate energy
    U_i = []
    E_i = []

    for i in range(3):
        U_i.append(sp.Matrix(sub_vector(U_calc, DOF[i])))
        E_i.append((1/2 * U_i[i].T * K_i[i] * U_i[i])[0] / 1000)

    # Caclulate centroid
    C = []

    for i in range(3):
        # Calculate area
        A = 0
        x_s, y_s = 0, 0
        vertices = parametric['rectangles'][i]

        for j in range(4):
            x_i = numeric[vertices[j]]["x"]  # x_i
            y_i = numeric[vertices[j]]["y"]  # y_i
            x_ipp = numeric[vertices[j != 3 and j+1 or 0]]["x"]  # x_i+1
            y_ipp = numeric[vertices[j != 3 and j+1 or 0]]["y"]  # y_i+1

            A += 1/2 * (x_i*y_ipp - x_ipp*y_i)
            x_s += (x_i + x_ipp) * (x_i*y_ipp - x_ipp*y_i)
            y_s += (y_i + y_ipp) * (x_i*y_ipp - x_ipp*y_i)

        x_s /= 6*A
        y_s /= 6*A

        C.append({
            "x": x_s,
            "y": y_s,
            "xi": transform[i]["xi"].subs({x: x_s, y: y_s}),
            "eta": transform[i]["eta"].subs({x: x_s, y: y_s}),
        })

    # Calculate sigma
    sigma_i = []

    for i in range(3):
        sigma_i.append(
            D * B[i].subs({xi: C[i]["xi"], eta: C[i]["eta"]}) * U_i[i]
        )

    # Merge newly created variables with original variables
    M = {
        "numeric": numeric,
        "parametric": parametric,

        "N": N,
        "transform": transform,

        "diff_local_N": diff_local_N,
        "diff_global_N": diff_global_N,
        "J": J,
        "inv_J": inv_J,
        "B": B,
        "D": D,
        "K_i": K_i,
        "DOF": DOF,

        "K": K,
        "U": U,
        "F": F,

        "U_sym": U_sym,
        "U_sym_kond": U_sym_kond,
        "F_sym": F_sym,

        "K_kond": K_kond,
        "K_kond_inv": K_kond_inv,
        "U_kond": U_kond,
        "F_kond": F_kond,

        "U_calc": U_calc,
        "F_calc": F_calc,
        "F_base": F_base,
        "F_reac": F_reac,

        "Delta_mm": Delta_mm,
        "Delta_um": Delta_mm * 1000,

        "U_i": U_i,
        "E_i": E_i,

        "C": C,
        "sigma_i": sigma_i,
        "sigma_2": sigma_i[1],
    }

    for k, v in M.items():
        Data[k] = v

    # Return variables
    return Data
