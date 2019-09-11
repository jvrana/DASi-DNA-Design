import numpy as np
import pandas as pd


class Globals(object):
    time_cost = 10.0
    material_modifier = 1.0


class PrimerParams(object):
    # @title Parameters
    time_cost = Globals.time_cost  # @param {type:"number"}
    material_modifier = Globals.material_modifier  # @param {type:"number"}
    min_anneal = 16  # @param {type:"number"}
    min_span = -300  # @param {type:"number"}

    # efficiency table
    eff_df = pd.DataFrame(
        [
            [0, 10, 0.0],
            [10, 15, 0.3],
            [15, 20, 0.6],
            [20, 30, 0.8],
            [30, 40, 0.9],
            [40, 50, 0.8],
            [50, 100, 0.75],
            [100, 120, 0.5],
            [120, 150, 0.3],
            [150, 250, 0.1],
            [250, 300, 0.0],
        ],
        columns=["min", "max", "efficiency"],
    )

    # primer cost table
    primer_df = pd.DataFrame(
        [
            [
                "Non-primer (special case)",
                min_anneal,
                min_anneal + 1,
                0.0,
                0.0,
                0.0,
            ],  # special case
            ["IDTPrimer", 16.0, 60.0, 0.0, 0.15, 1.5],
            ["IDTUltramer", 45.0, 200.0, 0.0, 0.40, 1.5],
        ],
        columns=["name", "min", "max", "base cost", "cost per bp", "time (days)"],
    )


class SynthesisParams(object):
    synthesis_df = pd.DataFrame(
        [
            [0, 1, 0, 0],
            [1, 100, np.inf, np.inf],
            [100, 500, 89.0, 3.0],
            [500, 750, 129.0, 3.0],
            [750, 1000, 149.0, 4.0],
            [1000, 1250, 209.0, 7.0],
            [1250, 1500, 249.0, 7.0],
            [1500, 1750, 289.0, 7.0],
            [1750, 2000, 329.0, 7.0],
            [2000, 2250, 399.0, 7.0],
        ],
        columns=["min", "max", "base", "time"],
    )

    time_cost = Globals.time_cost
    material_modifier = Globals.material_modifier
    step_size = 10
    left_span_range = (-300, 300)
