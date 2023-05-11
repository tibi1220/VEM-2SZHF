from sympy import latex, Number, Float, Matrix, preorder_traversal
import re


def round(num, numDecimalPlaces=None):
    if numDecimalPlaces is not None:
        num = float(f"{num:.{numDecimalPlaces}f}")
    return num


def prin_TeX(num, unit, dec=None):
    if dec is not None and dec != "":
        num = round(num, dec)
    print(f"\\SI{{{num}}}{{{unit}}}", end="")


def print_matrix(matrix, th=1e-16, mult=1, dec=""):
    (cols, rows) = matrix.shape
    for i in range(cols):
        for j in range(rows):
            num = matrix[i, j]
            if abs(num) < th:
                print(0)
            else:
                prin_TeX(num * mult, "", dec)
            if j == rows - 1:
                print("\\\\")
            else:
                print("&")


def print_raw_matrix(matrix):
    (cols, rows) = matrix.shape
    for i in range(cols):
        for j in range(rows):
            print(matrix[i, j])
            if j == rows - 1:
                print("\\\\")
            else:
                print("&")


# def round_expr(expr, num_digits):
#     return expr.xreplace({n: round(n, num_digits) for n in expr.atoms(Number)})


def round_expr(ex, digits):
    for a in preorder_traversal(ex):
        if isinstance(a, Float):
            ex = ex.subs(a, round(a, digits))
    return ex


def my_latex(var, digits=-1, **kwargs):
    if digits != -1:
        if type(var) == Matrix:
            (c, r) = var.shape

            for i in range(c):
                for j in range(r):
                    var[i, j] = round_expr(var[i, j], digits)
        else:
            var = round_expr(var, digits)

    return re.sub(
        r"\{,\}0(?=\D|$)",
        "",
        latex(
            var, decimal_separator="comma", mul_symbol=r"\,", full_prec=False, **kwargs
        ).replace("frac", "dfrac")
    )


def printer(t):
    v = t["variables"]

    def printSIDirect(num):
        prin_TeX(num["value"], num["unit"], num["dec"])

    def printSIVar(num):
        prin_TeX(v[num["name"]], num["unit"], num["dec"])

    def printVec(vec):
        print(v[vec["name"]][vec["index"]])

    return {
        "printSIDirect": printSIDirect,
        "printSIVar": printSIVar,
        "printVec": printVec,
    }
