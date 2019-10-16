```json
{
    "globals": {
        "time_cost": 10.0,
        "material_modified": 1.0
    },
    "primers": {
        "min_anneal": 16,
        "min_span": -300,
        "efficiency_table": [
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
            [250, 300, 0.0]
        ],
        "cost": [
            {
                "name": "IDTPrimer",
                "min": 16,
                "max": 60,
                "base_cost": 0,
                "cost_per_bp": 0.15,
                "time (days)": 1.5
            },
            {
                "name": "IDTUltramer",
                "min": 45,
                "max": 200,
                "base_cost": 0,
                "cost_per_bp": 0.40,
                "time (days)": 1.5
            }
        ]
    },
    "synthesis": {
        "step_size": 10,
        "left_span_range": [-300, 300],
        "cost": []
    }
}```