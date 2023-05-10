Data = {
    # Non specific info about construction
    "parametric": {
        1: {
            "code": 1,
            "p": {
                "b": 7,
                "e": 8,
                "l": {"x": 0, "y": -1},
                "v": {"x": 3.5, "y": -1, "l": "left"},
            },
            "F": {"x": [2, 8], "y": [5, 6]},
            "rectangles": [
                [5, 6, 3, 1],
                [3, 4, 2, 1],
                [7, 8, 2, 4],
            ],
            "free": [1, 2, 4, 5, 6, 7, 8, 9, 11, 13, 14, 16],
            "not_free": [3, 10, 12, 15],
            "distributed": [14, 16, -1],
            1: {"x": "0", "y": "b"},
            2: {"x": "c", "y": "b"},
            3: {"x": "c / 3", "y": "b - a"},
            4: {"x": "c - a", "y": "b - a"},
            5: {"x": "0", "y": "0"},
            6: {"x": "c / 3", "y": "0"},
            7: {"x": "c - a", "y": "0"},
            8: {"x": "c", "y": "0"},
        },
    },

    # Specific info about construction
    "numeric": lambda code, a, b, c: {
        1: {
            1: {"x": 0, "y": b},
            2: {"x": c, "y": b},
            3: {"x": c / 3, "y": b - a},
            4: {"x": c - a, "y": b - a},
            5: {"x": 0, "y": 0},
            6: {"x": c / 3, "y": 0},
            7: {"x": c - a, "y": 0},
            8: {"x": c, "y": 0},
        },
    }[code],

    # State
    "state": ["SF", "SA", "SF", "SA"],

    # Dimentions
    "a": [30, 25, 20, 35],  # mm
    "b": [80, 75, 70, 85],  # mm
    "c": [130, 120, 100, 140],  # mm
    "t": [8, 6, 4, 10],  # mm

    # Distributed load
    "p": [450, 350, 250, 550],  # MPa
}


def get_data(code):
    return {
        "parametric": Data["parametric"][code[1]],
        "numeric": Data["numeric"],

        "E": 180000,  # MPa
        "nu": 0.4,  # Poisson

        # B
        "state": Data["state"][code[2]-1],
        "a": Data["a"][code[2]-1],

        # C
        "p": Data["p"][code[3]-1],
        "b": Data["b"][code[3]-1],

        # D
        "c": Data["c"][code[4]-1],
        "t": Data["t"][code[4]-1],

        # Code information
        "code": code,
    }
