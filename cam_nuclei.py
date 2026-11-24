"""CAM H&E section -> tissue area, nucleus count, nucleus density, nucleus size / shape statistics.

  python cam_nuclei.py image.jpg [more images ...] --um 0.6595 --out results/

Pipeline
  1. exclude : a burned-in scale bar (long, thin, pure-white horizontal bar) and a box around it for its text labels;
               images without such a bar exclude nothing
  2. nuclei  : InstanSeg brightfield_nuclei (Apache-2.0, 16 MB) on the whole image (eval_medium_image: contrast is
               normalised once for the image, then tiled internally) -> one instance mask per nucleus;
               seed threshold 0.3 instead of the model's 0.7 (fewer missed nuclei, chosen on public val splits);
               optional 8-view test-time augmentation (--tta)
  3. tissue  : (stained OR textured) AND >= NEAR_MIN detected nuclei within NEAR_UM
                 stained  = optical density against the LOCAL glass brightness (uneven illumination)
                 textured = local OD standard deviation (pale loose CAM mesenchyme has fibres; glass is flat)
                 near a nucleus: removes glass that looks textured (camera banding, dust) and has no cells
               small holes (< HOLE_UM2) count as tissue, fragments < MIN_TISSUE_UM2 are dropped
  4. per nucleus: area, equivalent diameter, axes, elongation, circularity, solidity, mean haematoxylin OD
  5. summary : tissue area, count, density per mm2 of tissue, size distribution, local density map

Known limits (see README): chicken red blood cells are nucleated, so RBC nuclei (in vessels, haemorrhage) are
counted; results depend on the pixel size being right; CAM accuracy has not been checked against expert labels.
"""
import argparse, json, math, os, time
import numpy as np, cv2

UM_DEFAULT = 0.6595          # default pixel size: 10x objective of the system this was developed on
NEAR_UM = 35.0               # tissue needs >= NEAR_MIN nuclei within this distance
NEAR_MIN = 2
HOLE_UM2 = 2000.0            # holes in the tissue mask smaller than this count as tissue
MIN_TISSUE_UM2 = 5000.0      # tissue fragments smaller than this are dropped
GRID_UM = 100.0              # cell size of the local density map
SEED_THRESHOLD = 0.3         # InstanSeg seed threshold (model default 0.7 undercounts; chosen on public val splits)

_MODEL = {}


def instanseg():
    if "m" not in _MODEL:
        from instanseg import InstanSeg
        _MODEL["m"] = InstanSeg("brightfield_nuclei", verbosity=0)
    return _MODEL["m"]


# ----------------------------------------------------------------------------------------------- masks
def exclude_mask(rgb):
    """burned-in scale bar: a pure-white component much wider than tall (>= 80 px, aspect >= 8); its text labels sit
    above / below it, so a box around the bar (bar width x 3.5 bar-widths tall... capped) is excluded"""
    H, W = rgb.shape[:2]
    white = (rgb.min(-1) >= 248).astype(np.uint8)
    n, lab, st, _ = cv2.connectedComponentsWithStats(white, 8)
    m = np.zeros((H, W), bool)
    bars = [i for i in range(1, n) if st[i, 2] >= 80 and st[i, 3] <= 30 and st[i, 2] >= 8 * st[i, 3]
            and st[i, 4] >= 0.6 * st[i, 2] * st[i, 3]]                         # solid, long, thin
    if bars:                                                                  # one scale bar per image: the longest
        x, y, w, h, _ = st[max(bars, key=lambda i: st[i, 2])]
        pad_x, up, down = int(0.15 * w), int(0.45 * w), int(0.35 * w)
        m[max(0, y - up):min(H, y + h + down), max(0, x - pad_x):min(W, x + w + pad_x)] = True
    return m


def background(rgb, um, exclude=None):
    """local glass brightness per channel (uneven illumination / vignetting): bright envelope over ~250 um"""
    f = 8; small = cv2.resize(rgb, (rgb.shape[1] // f, rgb.shape[0] // f), interpolation=cv2.INTER_AREA).astype(np.float32)
    if exclude is not None and exclude.any():
        ex = cv2.resize(exclude.astype(np.uint8), (small.shape[1], small.shape[0]), interpolation=cv2.INTER_NEAREST) > 0
        if (~ex).any():
            small[ex] = np.median(small[~ex], axis=0)
    k = max(3, int(round(250 / um / f)) | 1)
    env = cv2.dilate(small, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k)))      # brightest nearby = glass
    env = cv2.GaussianBlur(env, (0, 0), k / 2)
    return cv2.resize(env, (rgb.shape[1], rgb.shape[0]), interpolation=cv2.INTER_LINEAR).clip(1, 255)


def _disk(r):
    return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))


def tissue_mask(rgb, um, centres, exclude=None):
    """(stained OR textured) AND near a nucleus; centres = (N, 2) array of (y, x) px"""
    H, W = rgb.shape[:2]
    ok = ~exclude if exclude is not None else np.ones((H, W), bool)
    bg = background(rgb, um, exclude)
    od = -np.log10(np.clip((rgb.astype(np.float32) + 1) / (bg + 1), 1e-3, 1.0)).sum(-1)
    s = max(3, int(round(10 / um)) | 1)
    mu = cv2.blur(od, (s, s)); tex = np.sqrt(np.maximum(cv2.blur(od * od, (s, s)) - mu * mu, 0))
    ods = cv2.GaussianBlur(od, (0, 0), max(1.0, 2.0 / um))
    glass = ok & (ods < np.percentile(ods[ok], 30))                            # lowest-OD 30 %: mostly glass
    t_od = float(np.clip(np.percentile(ods[glass], 99) + 0.03, 0.04, 0.25)) if glass.any() else 0.08
    t_tex = float(np.clip(np.percentile(tex[glass], 95) * 1.5, 0.012, 0.05)) if glass.any() else 0.02
    raw = ((ods > t_od) | (tex > t_tex)) & ok
    # cell support: >= NEAR_MIN nuclei within NEAR_UM (a lone false detection on glass does not make tissue);
    # counted on a 4x coarser grid for speed
    f = 4; h4, w4 = (H + f - 1) // f, (W + f - 1) // f
    pts = np.zeros((h4, w4), np.float32)
    if len(centres):
        c = (np.round(centres) // f).astype(int)
        np.add.at(pts, (c[:, 0].clip(0, h4 - 1), c[:, 1].clip(0, w4 - 1)), 1)
    r4 = max(1, int(round(NEAR_UM / um / f)))
    support = cv2.filter2D(pts, -1, _disk(r4).astype(np.float32), borderType=cv2.BORDER_CONSTANT) >= NEAR_MIN - 0.5
    near = cv2.resize(support.astype(np.uint8), (W, H), interpolation=cv2.INTER_NEAREST)
    m = raw & (near > 0)
    m = cv2.morphologyEx(m.astype(np.uint8), cv2.MORPH_CLOSE, _disk(max(1, int(round(8.0 / um))))) > 0
    m = cv2.morphologyEx(m.astype(np.uint8), cv2.MORPH_OPEN, _disk(max(1, int(round(4.0 / um))))) > 0
    # fill small enclosed holes (gaps in loose tissue), drop small fragments
    n, lab, st, _ = cv2.connectedComponentsWithStats((~m).astype(np.uint8), 8)
    touches = np.zeros(n, bool); touches[np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))] = True
    small = np.zeros(n, bool); small[1:] = st[1:, 4] * um * um < HOLE_UM2
    m |= small[lab] & ~touches[lab]
    n, lab, st, _ = cv2.connectedComponentsWithStats(m.astype(np.uint8), 8)
    keep = np.zeros(n, bool); keep[1:] = st[1:, 4] * um * um >= MIN_TISSUE_UM2
    return keep[lab] & ok, {"od_threshold": round(t_od, 3), "texture_threshold": round(t_tex, 4)}


# ----------------------------------------------------------------------------------------------- nuclei
def nuclei_labels(rgb, um, scale=1.0, seed=SEED_THRESHOLD):
    """instance label map (int32) of the whole image; InstanSeg normalises the image once, then tiles internally.
    seed = InstanSeg seed threshold (model default 0.7). 0.3 was chosen on PanNuke + PUMA validation: on the test
    splits point F1 PanNuke 0.769 -> 0.794, PUMA 0.915 -> 0.925, NuInsSeg (unseen) 0.881 -> 0.898 and the count
    ratio n_pred / n_true 0.75 / 0.92 / 0.84 -> 0.84 / 0.98 / 0.90 (default InstanSeg misses nuclei)."""
    import torch
    m = instanseg()
    kw = {"seed_threshold": float(seed)} if seed is not None else {}
    res = m.eval_medium_image(np.ascontiguousarray(rgb), pixel_size=um * scale, tile_size=1024,
                              batch_size=2 if torch.cuda.is_available() else 1, return_image_tensor=False, target="nuclei", **kw)
    lab = res[0] if isinstance(res, (tuple, list)) else res
    lab = np.asarray(lab.cpu() if hasattr(lab, "cpu") else lab).squeeze()
    lab = lab[0] if lab.ndim == 3 else lab
    return lab.astype(np.int32)


def nuclei_labels_tta(rgb, um, scale=1.0, seed=SEED_THRESHOLD, fuse_um=3.0, progress=None):
    """8-view test-time augmentation (flips x 90-degree rotations), union: a nucleus found in ANY view is kept.
    Each view's label map is mapped back exactly (inverse array transform); a nucleus from a later view is added only
    if its centre is > fuse_um from every nucleus already kept, and only on pixels not yet taken.
    On public test sets this raised point F1 (NuInsSeg 0.881 -> 0.905, PanNuke 0.769 -> 0.795) and moved the count
    ratio from 0.84 to 0.96 (InstanSeg alone misses nuclei). Cost: 8x inference."""
    from scipy.spatial import cKDTree
    out = None; cen_all = None
    views = [(f, k) for f in (False, True) for k in range(4)]
    for n, (flip, k) in enumerate(views):
        v = rgb[:, ::-1] if flip else rgb
        lab = nuclei_labels(np.ascontiguousarray(np.rot90(v, k)), um, scale, seed)
        lab = np.rot90(lab, -k); lab = np.ascontiguousarray(lab[:, ::-1] if flip else lab)
        ids, cen = centres_of(lab)
        if out is None:
            out, cen_all = lab, cen
        elif len(ids):
            d, _ = cKDTree(cen_all).query(cen, distance_upper_bound=fuse_um / um) if len(cen_all) else (np.full(len(cen), np.inf), None)
            new = ~np.isfinite(d)
            if new.any():
                base = int(out.max())
                lut = np.zeros(lab.max() + 1, np.int32); lut[ids[new]] = np.arange(base + 1, base + 1 + new.sum())
                add = lut[lab]; put = (add > 0) & (out == 0); out[put] = add[put]
                cen_all = np.concatenate([cen_all, cen[new]])
        if progress:
            progress((n + 1) / len(views))
    return out


def centres_of(lab):
    from scipy import ndimage as ndi
    ids = np.unique(lab[lab > 0])
    return ids, (np.array(ndi.center_of_mass(np.ones_like(lab), lab, ids)) if len(ids) else np.zeros((0, 2)))


def haematoxylin_od(rgb):
    from skimage.color import rgb2hed
    return rgb2hed(rgb)[..., 0].astype(np.float32)


def nucleus_table(lab, rgb, um):
    from skimage.measure import regionprops
    hema = haematoxylin_od(rgb); rows = []
    for p in regionprops(lab, intensity_image=hema):
        a = p.area * um * um; per = max(p.perimeter_crofton, 1e-6) * um
        rows.append({"id": int(p.label), "x_um": round(p.centroid[1] * um, 2), "y_um": round(p.centroid[0] * um, 2),
                     "area_um2": round(a, 2), "eq_diam_um": round(2 * math.sqrt(a / math.pi), 2),
                     "major_axis_um": round(p.axis_major_length * um, 2), "minor_axis_um": round(p.axis_minor_length * um, 2),
                     "elongation": round(p.axis_major_length / max(p.axis_minor_length, 1e-6), 3),
                     "circularity": round(min(1.0, 4 * math.pi * a / per ** 2), 3), "solidity": round(p.solidity, 3),
                     "haematoxylin_od": round(float(p.intensity_mean), 4)})
    return rows


# ----------------------------------------------------------------------------------------------- analysis
def analyse(rgb, um=UM_DEFAULT, scale=1.0, tta=False, progress=None, seed=SEED_THRESHOLD):
    t0 = time.time()
    exc = exclude_mask(rgb)
    lab = nuclei_labels_tta(rgb, um, scale, seed, progress=progress) if tta else nuclei_labels(rgb, um, scale, seed)
    t_nuc = time.time() - t0
    ids, cen = centres_of(lab)
    if len(ids):                                                              # drop nuclei on the scale bar / text
        ci = np.round(cen).astype(int)
        keep = ~exc[ci[:, 0], ci[:, 1]]; ids, cen = ids[keep], cen[keep]
    tissue, thr = tissue_mask(rgb, um, cen, exc)
    if len(ids):                                                              # keep nuclei whose centre is in tissue
        ci = np.round(cen).astype(int)
        keep = tissue[ci[:, 0], ci[:, 1]]
        lut = np.zeros(lab.max() + 1, np.int32); lut[ids[keep]] = np.arange(1, keep.sum() + 1)
        lab = lut[lab]
    else:
        lab = np.zeros_like(lab)
    rows = nucleus_table(lab, rgb, um)
    area_mm2 = float(tissue.sum()) * (um / 1000) ** 2
    a = np.array([r["area_um2"] for r in rows]) if rows else np.zeros(0)
    q = lambda v, p: round(float(np.percentile(v, p)), 2) if len(v) else None
    summary = {
        "image_px": [int(rgb.shape[1]), int(rgb.shape[0])], "um_per_px": um,
        "tissue_area_mm2": round(area_mm2, 4),
        "tissue_fraction_of_image_pct": round(100 * float(tissue.mean()), 1),
        "nuclei_count": len(rows),
        "nuclei_density_per_mm2": round(len(rows) / area_mm2, 1) if area_mm2 > 0 else None,
        "nucleus_area_um2_median": q(a, 50), "nucleus_area_um2_p10": q(a, 10), "nucleus_area_um2_p90": q(a, 90),
        "nucleus_area_um2_mean": round(float(a.mean()), 2) if len(a) else None,
        "nucleus_eq_diam_um_median": q(np.array([r["eq_diam_um"] for r in rows]), 50) if rows else None,
        "nucleus_circularity_median": q(np.array([r["circularity"] for r in rows]), 50) if rows else None,
        "elongated_nuclei_pct": round(100 * float(np.mean([r["elongation"] > 2.0 for r in rows])), 1) if rows else None,
        "high_accuracy_tta": bool(tta), "seed_threshold": seed, "scale_bar_excluded": bool(exc.any()), **thr,
        "seconds_nuclei": round(t_nuc, 1), "seconds_total": round(time.time() - t0, 1),
    }
    return summary, rows, lab, tissue, exc


def density_map(rows, tissue, um, grid_um=GRID_UM):
    """nuclei per mm2 in grid cells of grid_um; NaN where a cell is < 50 % tissue"""
    g = max(1, int(round(grid_um / um))); H, W = tissue.shape
    gh, gw = math.ceil(H / g), math.ceil(W / g)
    cnt = np.zeros((gh, gw))
    for r in rows:
        cnt[min(gh - 1, int(r["y_um"] / um) // g), min(gw - 1, int(r["x_um"] / um) // g)] += 1
    t = np.pad(tissue, ((0, gh * g - H), (0, gw * g - W)))
    frac = t.reshape(gh, g, gw, g).mean((1, 3))
    area = frac * (g * um / 1000) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(frac >= 0.5, cnt / area, np.nan)


def overlay(rgb, lab, tissue, exc=None):
    o = rgb.copy()
    lf = lab.astype(np.float32)
    edge = (lab > 0) & ((cv2.erode(lf, np.ones((3, 3), np.uint8)) != lf) | (cv2.dilate(lf, np.ones((3, 3), np.uint8)) != lf))
    o[edge] = (255, 215, 0)                                                   # nuclei: yellow outlines
    cs, _ = cv2.findContours(tissue.astype(np.uint8), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(o, cs, -1, (0, 160, 255), max(2, rgb.shape[1] // 1000))  # tissue boundary: blue
    if exc is not None and exc.any():
        o[exc] = (o[exc] * 0.5 + np.array([128, 128, 128]) * 0.5).astype(np.uint8)   # excluded: greyed
    return o


def density_image(d, size):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    cmap = plt.get_cmap("magma").copy(); cmap.set_bad("#d9d9d9")
    im = ax.imshow(d, cmap=cmap); ax.set_axis_off()
    fig.colorbar(im, ax=ax, fraction=0.035, label="nuclei / mm²")
    ax.set_title(f"local nucleus density ({GRID_UM:.0f} µm cells; grey = < 50 % tissue)", fontsize=9, pad=6)
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.canvas.draw()
    img = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy(); plt.close(fig)
    return img


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("images", nargs="+"); ap.add_argument("--um", type=float, default=UM_DEFAULT, help="micrometres per pixel")
    ap.add_argument("--scale", type=float, default=1.0, help="InstanSeg pixel-size factor (advanced)")
    ap.add_argument("--seed", type=float, default=SEED_THRESHOLD, help="InstanSeg seed threshold (model default 0.7)")
    ap.add_argument("--tta", action="store_true", help="high-accuracy mode: 8-view union (fewer missed nuclei, 8x slower)")
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); allsum = []
    import csv
    for p in a.images:
        rgb = cv2.cvtColor(cv2.imread(p), cv2.COLOR_BGR2RGB); n = os.path.splitext(os.path.basename(p))[0]
        s, rows, lab, tissue, exc = analyse(rgb, a.um, a.scale, a.tta, seed=a.seed)
        s["image"] = os.path.basename(p); allsum.append(s)
        cv2.imwrite(os.path.join(a.out, n + "_overlay.jpg"), cv2.cvtColor(overlay(rgb, lab, tissue, exc), cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 90])
        if rows:
            with open(os.path.join(a.out, n + "_nuclei.csv"), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(json.dumps(s, ensure_ascii=False), flush=True)
    json.dump(allsum, open(os.path.join(a.out, "summary.json"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
