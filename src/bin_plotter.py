import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def load_bins(csv_path):
    df = pd.read_csv(csv_path)
    # Compute von Mises optional
    Sxx = df["Sxx"].values
    Syy = df["Syy"].values
    Szz = df["Szz"].values

    df["von_mises"] = np.sqrt(
        0.5 * ((Sxx - Syy)**2 + (Syy - Szz)**2 + (Szz - Sxx)**2)
    )

    return df

def set_axes_equal(ax):
    """
    Make x, y, z axes have equal scale so voxels appear as true cubes.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_bins_3d(ax, df, title, color_by="von_mises", cmap="inferno"):

    x = df["x"].values
    y = df["y"].values
    z = df["z"].values
    c = df[color_by].values

    sc = ax.scatter(
        x, y, z,
        c=c,
        cmap=cmap,
        s=40,
        marker='o',
        alpha=0.9,
        edgecolor='k',
        linewidths=0.3
    )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    plt.colorbar(sc, ax=ax, shrink=0.6, pad=0.1, label=color_by)


def compare_voxel_bins(truth_csv, pred_csv, color_by="von_mises"):
    truth = load_bins(truth_csv)
    pred = load_bins(pred_csv)

    fig = plt.figure(figsize=(12, 6))

    # Left plot: Ground truth
    ax1 = fig.add_subplot(121, projection="3d")
    plot_bins_3d(ax1, truth, "Ground Truth Voxels", color_by=color_by)
    set_axes_equal(ax1)

    # Right plot: Predicted
    ax2 = fig.add_subplot(122, projection="3d")
    plot_bins_3d(ax2, pred, "Predicted Voxels", color_by=color_by)
    set_axes_equal(ax2)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    compare_voxel_bins("~/code/research/residual_stress/residual_stress_estimator/build/calculated_chunks.csv", 
                      "~/code/research/residual_stress/residual_stress_estimator/build/ground_truth_chunks.csv")
