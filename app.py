"""Gradio web app (English / Deutsch): upload an H&E image -> tissue area, nucleus count, density and size statistics.

  python app.py            then open http://127.0.0.1:7860
"""
import csv, json, os, tempfile
import numpy as np
import gradio as gr
import cam_nuclei as cn

HERE = os.path.dirname(os.path.abspath(__file__))
LANGS = ["English", "Deutsch"]

T = {
    "English": {
        "intro": "## Nucleus count, density and size in H&E sections (CAM assay)\n"
                 "Upload an H&E image and enter its pixel size. **Yellow** = detected nuclei, **blue** = tissue boundary, "
                 "**grey** = excluded burned-in scale bar.\n\n"
                 "⚠️ Research demo, not a validated measurement tool. Chicken red blood cells are nucleated, so nuclei "
                 "inside vessels are counted too. Accuracy on CAM sections has not yet been checked against expert labels.",
        "image": "H&E image", "um": "Pixel size (µm / px)",
        "tta": "High-accuracy mode: detect in 8 orientations and merge (fewer missed nuclei, ~8× slower)",
        "run": "Analyse", "examples": "Examples (public PUMA melanoma images, CC0, 0.2263 µm/px)",
        "overlay": "Detected nuclei", "table": "Summary", "metric": "Metric", "value": "Value",
        "density": "Local nucleus density", "hist": "Nucleus size distribution",
        "files": "Download: per-nucleus measurements (CSV) and summary (JSON)",
        "err_img": "Please upload an image first.", "err_um": "The pixel size must be greater than 0.",
        "busy": "detecting nuclei",
        "rows": {
            "tissue_area_mm2": "Tissue area (mm²)", "tissue_fraction_of_image_pct": "Tissue fraction of image (%)",
            "nuclei_count": "Nuclei count", "nuclei_density_per_mm2": "Nucleus density (per mm² tissue)",
            "nucleus_area_um2_median": "Nucleus area, median (µm²)", "nucleus_area_um2_p10": "Nucleus area, 10th percentile (µm²)",
            "nucleus_area_um2_p90": "Nucleus area, 90th percentile (µm²)", "nucleus_eq_diam_um_median": "Equivalent diameter, median (µm)",
            "nucleus_circularity_median": "Circularity, median (0–1)", "elongated_nuclei_pct": "Elongated nuclei, long/short axis > 2 (%)",
            "um_per_px": "Pixel size (µm / px)", "high_accuracy_tta": "High-accuracy mode", "seed_threshold": "InstanSeg seed threshold",
            "scale_bar_excluded": "Burned-in scale bar excluded", "seconds_total": "Processing time (s)"},
    },
    "Deutsch": {
        "intro": "## Anzahl, Dichte und Größe von Zellkernen in H&E-Schnitten (CAM-Assay)\n"
                 "Laden Sie ein H&E-Bild hoch und geben Sie die Pixelgröße an. **Gelb** = erkannte Zellkerne, "
                 "**Blau** = Gewebegrenze, **Grau** = ausgeschlossener eingebrannter Maßstabsbalken.\n\n"
                 "⚠️ Forschungs-Demo, kein validiertes Messwerkzeug. Erythrozyten des Huhns sind kernhaltig, daher werden "
                 "auch Kerne in Gefäßen mitgezählt. Die Genauigkeit auf CAM-Schnitten wurde noch nicht mit "
                 "Experten-Annotationen geprüft.",
        "image": "H&E-Bild", "um": "Pixelgröße (µm / px)",
        "tta": "Hochgenauigkeitsmodus: Erkennung in 8 Orientierungen, zusammengeführt (weniger übersehene Kerne, ca. 8× langsamer)",
        "run": "Analysieren", "examples": "Beispiele (öffentliche PUMA-Melanom-Bilder, CC0, 0,2263 µm/px)",
        "overlay": "Erkannte Zellkerne", "table": "Zusammenfassung", "metric": "Kennwert", "value": "Wert",
        "density": "Lokale Kerndichte", "hist": "Verteilung der Kerngröße",
        "files": "Download: Messwerte je Zellkern (CSV) und Zusammenfassung (JSON)",
        "err_img": "Bitte zuerst ein Bild hochladen.", "err_um": "Die Pixelgröße muss größer als 0 sein.",
        "busy": "Zellkerne werden erkannt",
        "rows": {
            "tissue_area_mm2": "Gewebefläche (mm²)", "tissue_fraction_of_image_pct": "Gewebeanteil am Bild (%)",
            "nuclei_count": "Anzahl Zellkerne", "nuclei_density_per_mm2": "Kerndichte (pro mm² Gewebe)",
            "nucleus_area_um2_median": "Kernfläche, Median (µm²)", "nucleus_area_um2_p10": "Kernfläche, 10. Perzentil (µm²)",
            "nucleus_area_um2_p90": "Kernfläche, 90. Perzentil (µm²)", "nucleus_eq_diam_um_median": "Äquivalentdurchmesser, Median (µm)",
            "nucleus_circularity_median": "Rundheit, Median (0–1)", "elongated_nuclei_pct": "Längliche Kerne, Längs-/Querachse > 2 (%)",
            "um_per_px": "Pixelgröße (µm / px)", "high_accuracy_tta": "Hochgenauigkeitsmodus", "seed_threshold": "InstanSeg-Seed-Schwelle",
            "scale_bar_excluded": "Eingebrannter Maßstabsbalken ausgeschlossen", "seconds_total": "Rechenzeit (s)"},
    },
}


def size_hist(rows, lang):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 3))
    if rows:
        ax.hist([r["area_um2"] for r in rows], bins=50, color="#4c78a8")
    de = lang == "Deutsch"
    ax.set_xlabel("Kernfläche (µm²)" if de else "nucleus area (µm²)"); ax.set_ylabel("Anzahl" if de else "count")
    ax.set_title(T[lang]["hist"]); fig.tight_layout()
    return fig


def fmt(v, lang):
    if isinstance(v, bool):
        return ("ja" if v else "nein") if lang == "Deutsch" else ("yes" if v else "no")
    if isinstance(v, float) and lang == "Deutsch":
        return f"{v:g}".replace(".", ",")
    return v


def make_table(s, lang):
    import pandas as pd
    t = T[lang]
    rows = [[label, fmt(s[k], lang)] for k, label in t["rows"].items() if k in s] if s else []
    return pd.DataFrame(rows, columns=[t["metric"], t["value"]])


def run(img, um, tta, lang, progress=gr.Progress()):
    t = T[lang]
    if img is None:
        raise gr.Error(t["err_img"])
    if not um or um <= 0:
        raise gr.Error(t["err_um"])
    rgb = np.ascontiguousarray(img[..., :3].astype(np.uint8))
    progress(0.05, desc=t["busy"])
    s, rows, lab, tissue, exc = cn.analyse(rgb, float(um), tta=bool(tta), progress=lambda f: progress(0.05 + 0.85 * f, desc=t["busy"]))
    ov = cn.overlay(rgb, lab, tissue, exc)
    dm = cn.density_image(cn.density_map(rows, tissue, float(um)), (720, 540))
    table = make_table(s, lang)
    d = tempfile.mkdtemp(); csv_p = os.path.join(d, "nuclei.csv"); js_p = os.path.join(d, "summary.json")
    with open(csv_p, "w", newline="") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump(s, open(js_p, "w"), indent=1)
    return ov, dm, table, size_hist(rows, lang), [csv_p, js_p], s


EX_DIR = os.path.join(HERE, "examples")
EX_UM = 0.2263                      # PUMA regions of interest
examples = [os.path.join(EX_DIR, f) for f in sorted(os.listdir(EX_DIR)) if f.lower().endswith((".jpg", ".png"))] \
    if os.path.isdir(EX_DIR) else []


def ex_label(lang, i):
    return f"{'Beispiel' if lang == 'Deutsch' else 'Example'} {i + 1}"


def load_example(i):
    import cv2
    return cv2.cvtColor(cv2.imread(examples[i]), cv2.COLOR_BGR2RGB), EX_UM


def relabel(lang, last):
    t = T[lang]
    return (gr.update(value=t["intro"]), gr.update(label=t["image"]), gr.update(label=t["um"]), gr.update(label=t["tta"]),
            gr.update(value=t["run"]), gr.update(label=t["overlay"]), gr.update(label=t["table"], value=make_table(last, lang)),
            gr.update(label=t["density"]), gr.update(label=t["hist"]), gr.update(label=t["files"]), gr.update(value=f"**{t['examples']}**"),
            *[gr.update(value=ex_label(lang, i)) for i in range(len(examples))])

with gr.Blocks(title="CAM nuclei analysis") as demo:
    t0 = T["English"]
    lang = gr.Radio(LANGS, value="English", label="Language / Sprache")
    intro = gr.Markdown(t0["intro"])
    with gr.Row():
        with gr.Column(scale=1):
            inp = gr.Image(label=t0["image"], type="numpy")
            um = gr.Number(value=cn.UM_DEFAULT, label=t0["um"], precision=4)
            tta = gr.Checkbox(value=False, label=t0["tta"])
            btn = gr.Button(t0["run"], variant="primary")
            ex_title = gr.Markdown(f"**{t0['examples']}**")
            with gr.Row():
                ex_btns = [gr.Button(ex_label("English", i), size="sm") for i in range(len(examples))]
            for i, b in enumerate(ex_btns):
                b.click(lambda i=i: load_example(i), None, [inp, um])
        with gr.Column(scale=2):
            out_ov = gr.Image(label=t0["overlay"])
            out_tab = gr.Dataframe(value=make_table(None, "English"), label=t0["table"], wrap=True)
    with gr.Row():
        out_dm = gr.Image(label=t0["density"])
        out_h = gr.Plot(label=t0["hist"])
    out_f = gr.File(label=t0["files"], file_count="multiple")
    last = gr.State(None)                       # last summary, so a language switch can re-render the table
    btn.click(run, [inp, um, tta, lang], [out_ov, out_dm, out_tab, out_h, out_f, last])
    lang.change(relabel, [lang, last], [intro, inp, um, tta, btn, out_ov, out_tab, out_dm, out_h, out_f, ex_title, *ex_btns])

if __name__ == "__main__":
    demo.launch()
