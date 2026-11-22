import glob
import os

import matplotlib.pyplot as plt
import numpy as np

experiments = {
    "lfsMin = 0.99": "runs/actual_good_tests/spot_test_lfs_1",
    "lfsMin = 0.5": "runs/actual_good_tests/spot_test_lfs_2",
    "lfsMin = 0.09": "runs/actual_good_tests/spot_test_lfs_3",
    "lfsMin = 0.05": "runs/actual_good_tests/spot_test_lfs_4",
    "lfsMin = 0.009": "runs/actual_good_tests/spot_test_lfs_5"
}

EVERY = {
    "cycle": (10, 0),
}

AS_PERCENT = True    # False: plot the number of flipped triangles instead of the percentage


def read_csv(path, dtype=float):
    with open(path) as f:
        sep = ";" if ";" in f.readline() else ","
    return np.loadtxt(path, delimiter=sep, skiprows=1, dtype=dtype)


def flipped(P, F):
    a, b, c = P[F[:, 0]], P[F[:, 1]], P[F[:, 2]]
    normal = np.cross(b - a, c - a)
    outward = (a + b + c) / 3 - P.mean(0)
    s = (normal * outward).sum(1)
    return s * np.sign(np.median(s)) < 0


def count_flipped(run, step=1, offset=0):
    F = read_csv(os.path.join(run, "faces.csv"), dtype=int)
    frames = sorted(glob.glob(os.path.join(run, "frames", "pos_*.csv")))

    rows = []
    for path in frames:
        it = int(os.path.basename(path)[4:-4])

        if (it - offset) % step != 0:
            continue

        P = read_csv(path)
        rows.append((it, flipped(P, F).sum()))

    return np.array(rows), len(F)


if __name__ == "__main__":
    plt.figure(figsize=(8, 5))

    for label, run in experiments.items():
        step, offset = EVERY.get(label, (1, 0))
        rows, n_tri = count_flipped(run, step, offset)

        values = 100.0 * rows[:, 1] / n_tri if AS_PERCENT else rows[:, 1]
        plt.plot(rows[:, 0], values, label=label)

        print(f"{label}: {int(rows[-1, 1])} flipped of {n_tri} at iteration {int(rows[-1, 0])}")

    plt.xlabel("Iterations")
    plt.ylabel("Flipped triangles (%)" if AS_PERCENT else "Flipped triangles")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.savefig("flipped_triangles.pdf", bbox_inches="tight")
    plt.show()