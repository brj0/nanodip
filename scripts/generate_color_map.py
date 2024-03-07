import os

import numpy as np
from numba import njit

from nanodip.config import TMP
from nanodip.data import Reference, Sample
from nanodip.epidip import top_variable_cpgs
from nanodip.plots import UMAPData

ONE_THIRD = 1.0 / 3.0
ONE_SIXTH = 1.0 / 6.0
TWO_THIRD = 2.0 / 3.0


@njit
def _v(m1, m2, hue):
    hue = hue % 1.0
    if hue < ONE_SIXTH:
        return m1 + (m2 - m1) * hue * 6.0
    if hue < 0.5:
        return m2
    if hue < TWO_THIRD:
        return m1 + (m2 - m1) * (TWO_THIRD - hue) * 6.0
    return m1


@njit
def hls_to_rgb(h, l, s):
    if s == 0.0:
        return l, l, l
    if l <= 0.5:
        m2 = l * (1.0 + s)
    else:
        m2 = l + s - (l * s)
    m1 = 2.0 * l - m2
    return (
        _v(m1, m2, h + ONE_THIRD),
        _v(m1, m2, h),
        _v(m1, m2, h - ONE_THIRD),
    )


@njit
def _random_color():
    """Pseudorandom color that is good visable.
    Args:
        var: string to hash
    Returns:
        Tripple of rgb color.
    """
    hue = np.random.randint(0, 365)
    saturation = np.random.randint(10, 101)
    lightness = np.random.randint(30, 101)
    # plotly needs rgb values, not hsl
    rgb_frac = hls_to_rgb(hue / 364, lightness / 100, saturation / 100)
    rgb = [int(255 * x) for x in rgb_frac]
    return rgb


@njit
def color_diff(col0, col1):
    rmean = (col0[0] + col1[0]) / 2
    dr = col0[0] - col1[0]
    dg = col0[1] - col1[1]
    db = col0[2] - col1[2]
    return (
        (512 + rmean) * dr * dr / 256
        + 4 * dg * dg
        + (767 - rmean) * db * db / 256
    )


def generate_discrete_color_map(
    names, coords, output_path=os.path.join(TMP, "discrete_color_map.py")
):
    """
    Generate a discrete color map to ensure neighboring points have
    well-distinguishable colors.
    The algorithm continuously improves the coloring and regularly saves
    results to disk.

    Args:
        names (list of str): List of methylation classes.
        coords (DataFrame): DataFrame containing x- and y-values
        output_path (str): Path of the output file.

    Returns:
        None
    """
    unique_names = list(set(names))
    all_dist = np.sum((coords[:, np.newaxis] - coords) ** 2, axis=-1)
    n_names = len(unique_names)
    dist = np.zeros([n_names, n_names])
    for i in range(n_names):
        index_i = names == unique_names[i]
        for j in range(i, n_names):
            index_j = names == unique_names[j]
            # d = np.min(all_dist[np.ix_(index_i, index_j)])
            d = np.min(all_dist[index_i, :][:, index_j])
            dist[i, j] = d
            dist[j, i] = d

    color_map = np.array(
        [_random_color() for i in range(n_names)], dtype=np.int32
    )
    diameter = np.max(np.sqrt(dist))

    AVG_DIST = 80000
    dc_min = AVG_DIST // 100
    delta_dc_min = (AVG_DIST / 10 - dc_min) // 100
    dist_min = diameter / 100
    delta_dist_min = (diameter / 3 - dist_min) / 100
    n_fail = n_names

    @njit
    def improve_color_map(color_map, dc_min, dist_min, n_fail):
        prev_color = np.zeros(3, dtype=np.int32)
        prev_idx = -1
        while True:
            too_similar = np.zeros(n_names, dtype=np.int32)
            for i in range(n_names):
                for j in range(i + 1, n_names):
                    if dist[i, j] < dist_min:
                        dc = color_diff(color_map[i], color_map[j])
                        if dc < dc_min:
                            too_similar[i] = 1
                            too_similar[j] = 1
            curr_n_fail = np.sum(too_similar)
            if curr_n_fail < n_fail:
                return curr_n_fail
            # Roll back changes if new color map is worse
            if curr_n_fail > n_fail and prev_idx != -1:
                color_map[prev_idx] = prev_color
            indices = np.where(too_similar)[0]
            rand_idx = np.random.choice(indices)
            # Save for next round
            prev_idx = rand_idx
            prev_color = color_map[rand_idx]
            color_map[rand_idx] = _random_color()

    def save_to_disk(color_map):
        result = {
            key: f"'rgb{color_map[key]}'" for key in sorted(color_map.keys())
        }
        with open(output_path, "w") as f:
            f.write("# Dictionary containing color data\n")
            f.write("color_map = {\n")
            for key, value in result.items():
                f.write(f"    '{key}': {value},\n")
            f.write("}\n")

    while True:
        n_fail = improve_color_map(color_map, dc_min, dist_min, n_fail)
        if n_fail == 0:
            color_map_dict = {
                x: tuple(y) for x, y in zip(unique_names, color_map)
            }
            save_to_disk(color_map_dict)
            if np.random.randint(2) == 0:
                dc_min += delta_dc_min
            else:
                dist_min += delta_dist_min
            print("New round: dc_min =", dc_min, "dist_min =", dist_min)
            n_fail = n_names
        else:
            print("    ", n_fail)


REFERENCE_NAME = "AllIDATv2_20210804_HPAP_Sarc"
NR_TOP_CPGS = 50000

reference = Reference(REFERENCE_NAME)
top_cpgs = top_variable_cpgs(reference, NR_TOP_CPGS)
sample = Sample.from_cpgs(top_cpgs)
umap = UMAPData(reference, sample)
umap.make_umap_plot()
umap.plot.show()


names = umap.umap_df.methylation_class.values
coords = umap.umap_df[["x", "y"]].values


generate_discrete_color_map(names, coords)
