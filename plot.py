import matplotlib.pyplot as plt
import pandas as pd

experiments = {
    r"$\delta$ = 5e-03": "runs/actual_good_tests/spot_test_1/metrics.csv",
    r"$\delta$ = 1e-03": "runs/actual_good_tests/spot_test_3/metrics.csv",
    r"$\delta$ = 5e-04": "runs/actual_good_tests/spot_test_2/metrics.csv",
    r"$\delta$ = 1e-04": "runs/actual_good_tests/spot_test_4/metrics.csv",
    "cycle": "runs/actual_good_tests/spot_test_cycle/metrics.csv"
}

CYCLE = {
    "cycle": 10,
}


def per_cycle(data, length, how):
    """One row per finished cycle.

    how = "last": the value at the end of the cycle (the cleaned-up state: sphericity, min angle)
    how = "sum" : the total over the cycle (time)
    how = "mean": the average over the cycle (steps per walk, error)
    """
    data = data[data["iter"] > 0].copy()                  
    data["cycle"] = (data["iter"] - 1) // length          

    counts = data.groupby("cycle")["iter"].transform("size")
    data = data[counts == length]                         

    groups = data.groupby("cycle")
    ends = groups["iter"].max()

    if how == "last":
        values = groups.apply(lambda g: g.sort_values("iter").iloc[-1])
    elif how == "sum":
        values = groups.sum(numeric_only=True)
    else:
        values = groups.mean(numeric_only=True)

    values["iter"] = ends.values
    return values


def plot_column(column, ylabel, location, how="last"):
    plt.figure(figsize=(8, 5))                 # new window
    for label, filename in experiments.items():
        data = pd.read_csv(filename)

        if label in CYCLE:
            data = per_cycle(data, CYCLE[label], how)

        plt.plot(data["iter"], data[column], label=label, marker="o" if label in CYCLE else None, markersize=3)
    plt.xlabel("Iterations")
    plt.ylabel(ylabel)
    plt.legend(loc=location)
    plt.grid(True)


plot_column("sphericity", "Sphericity", "lower right", how="last")
plot_column("avg_walk_steps", "Steps per walk", "upper right", how="mean")
plot_column("time_sec", "Time per iteration", "upper right", how="sum")
plot_column("err_rel_bbox", "RMS", "upper right", how="mean")
plot_column("min_angle_deg", "min_deg", "upper right", how="last")

plt.show()