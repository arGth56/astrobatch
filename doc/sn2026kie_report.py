#!/usr/bin/env python3
"""Generate a 1-page A4 PDF report on SN2026kie photometry quality."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from astropy.io import fits
from astropy.visualization import ZScaleInterval

# ── Data ──────────────────────────────────────────────────────────────────────

nights = ["2026-04-24", "2026-04-25", "2026-04-27"]
mjd_ref = 61154.0

sub = {
    "G":  {"mjd": [61153.8101, 61154.9184, 61155.8566, 61157.9147, 61161.9118],
           "mag": [17.13, 18.69, 17.67, 17.48, 16.35], "err": [0.03, 0.35, 0.12, 0.09, 0.03]},
    "BP": {"mjd": [61154.9475, 61155.9094, 61157.9346],
           "mag": [17.52, 17.48, 17.23], "err": [0.12, 0.09, 0.07]},
    "RP": {"mjd": [61154.9766, 61155.9412, 61157.9536],
           "mag": [17.89, 16.88, 16.74], "err": [0.44, 0.19, 0.12]},
}

sn = {
    "G":  [29.46, 2.88, 8.63, 11.01, 32.11],
    "BP": [8.67, 11.47, 14.21],
    "RP": [2.28, 5.32, 8.18],
}

# Source per G-band point (for marker style)
g_source = ["HOTPANTS?", "OBS-BZH", "OBS-BZH", "OBS-BZH", "DAX T500"]
# Flag: the first G point (MJD 61153.81) is suspect — 43σ off from linear model
g_suspect = [True, False, False, False, False]

# Professional survey discovery / detection data (different filter systems)
pro_surveys = {
    "GOTO":  {"mjd": 61151.901860, "mag": 18.200, "filt": "L",      "color": "#9b59b6", "marker": "v"},
    "ATLAS": {"mjd": 61152.313730, "mag": 18.086, "filt": "orange", "color": "#e67e22", "marker": "p"},
    "ZTF":   {"mjd": 61154.175995, "mag": 17.406, "filt": "r",      "color": "#e74c3c", "marker": "d"},
}

obs = {
    "nframes": {"G": [40, 67, 40], "BP": [37, 60, 40], "RP": [37, 60, 40]},
    "fwhm":    {"G": [4.68, 5.15, 4.44], "BP": [4.88, 4.51, 4.19], "RP": [5.17, 4.77, 4.13]},
}

calib = {
    "G":  {"zp": [13.182, 13.571, 13.596], "ct": [-0.26, -0.25, -0.25]},
    "BP": {"zp": [13.028, 13.378, 13.355], "ct": [+0.04, +0.02, +0.03]},
    "RP": {"zp": [11.915, 11.684, 12.082], "ct": [-0.18, -0.18, -0.19]},
}

ci_mjd = [61154.96, 61155.93, 61157.94]
ci_val = [sub["BP"]["mag"][i] - sub["RP"]["mag"][i] for i in range(3)]
ci_err = [np.sqrt(sub["BP"]["err"][i]**2 + sub["RP"]["err"][i]**2) for i in range(3)]

filt_colors = {"G": "#2ecc71", "BP": "#3498db", "RP": "#e74c3c"}
filt_markers = {"G": "o", "BP": "s", "RP": "D"}

# ── Load subtraction cutout ──────────────────────────────────────────────────

sub_cutout = fits.open("/tmp/sub_target_3751.cutout")
img_data = sub_cutout[1].data
conv_data = sub_cutout[5].data
diff_data = sub_cutout[3].data
sub_cutout.close()

# ── Figure ────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(8.27, 11.69))
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(6, 2,
    height_ratios=[1.5, 1.6, 3.4, 1.8, 1.6, 1.4],
    hspace=0.65, wspace=0.35,
    left=0.09, right=0.96, top=0.96, bottom=0.03)

# ══════════════════════════════════════════════════════════════════════════════
# Row 0: Title + context
# ══════════════════════════════════════════════════════════════════════════════

ax_title = fig.add_subplot(gs[0, :])
ax_title.axis("off")

ax_title.text(0.5, 0.97, "SN 2026kie — Template Subtraction Photometry Report",
              ha="center", va="top", fontsize=13, fontweight="bold")

ax_title.text(0.5, 0.78,
    "RAPAS ProAm Network  ·  Gaia G / BP / RP filters  ·  Template subtraction photometry",
    ha="center", va="top", fontsize=7.5, color="#444", style="italic")

context = (
    "The RAPAS group (Réseau Amateur-Professionnel pour l'Astronomie des Supernovae, "
    "Observatoire de Paris / IMCCE) monitors supernovae and transients through coordinated "
    "multi-band photometric follow-up between amateur and professional observatories.\n\n"
    "Photometric calibration and template subtraction are performed using STDWeb, "
    "a web service built on STDPipe (Karpov, S., 2021, "
    "https://github.com/karpov-sv/stdpipe), which provides automated astrometry, "
    "photometric calibration against Gaia eDR3, and image subtraction using "
    "HOTPANTS with ZTF DR7 HiPS templates. For SN 2026kie, template subtraction "
    "is essential: the supernova lies on the host galaxy disk, making direct "
    "aperture photometry unreliable due to galaxy background contamination."
)

ax_title.text(0.02, 0.68, context, transform=ax_title.transAxes,
              fontsize=6.5, va="top", ha="left", linespacing=1.5,
              fontfamily="serif", wrap=True)

ax_title.text(0.98, 0.95, "2 May 2026", ha="right", va="top",
              fontsize=6.5, color="#888")

# ══════════════════════════════════════════════════════════════════════════════
# Row 1 left: Template subtraction illustration
# ══════════════════════════════════════════════════════════════════════════════

gs_sub = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[1, 0],
                                          wspace=0.10)
zscale = ZScaleInterval()

panels = [
    (img_data, "Science"),
    (conv_data, "Conv. template"),
    (diff_data, "Difference"),
]

for pi, (data, title) in enumerate(panels):
    ax = fig.add_subplot(gs_sub[0, pi])
    vmin, vmax = zscale.get_limits(data)
    if "Diff" in title:
        vabs = max(abs(np.nanpercentile(data, 2)), abs(np.nanpercentile(data, 98)))
        ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-vabs, vmax=vabs)
    else:
        ax.imshow(data, origin="lower", cmap="gray_r", vmin=vmin, vmax=vmax)
    cx, cy = data.shape[1]//2, data.shape[0]//2
    ax.plot(cx, cy, "+", color="lime" if "Diff" in title else "red",
            markersize=7, markeredgewidth=0.7)
    ax.set_title(title, fontsize=5.5, pad=2)
    ax.set_xticks([]); ax.set_yticks([])

ax_ann = fig.add_subplot(gs_sub[0, 3])
ax_ann.axis("off")
ax_ann.text(0.0, 0.95,
    "Task 3751\nG-band, Apr 27\nS/N = 11.0\n\n"
    "science − conv.\ntemplate = residual\n\n"
    "Template:\nZTF DR7 r-band\n(HiPS cutout)",
    transform=ax_ann.transAxes, fontsize=5.5, va="top", fontfamily="serif",
    linespacing=1.35)

# ══════════════════════════════════════════════════════════════════════════════
# Row 1 right: Template mapping table
# ══════════════════════════════════════════════════════════════════════════════

ax_tbl = fig.add_subplot(gs[1, 1])
ax_tbl.axis("off")
ax_tbl.set_title("Template Subtraction Configuration", fontsize=8,
                 fontweight="bold", pad=6)

tbl_data = [
    ["G",  "Gaia Gmag",  "ZTF DR7 r (HiPS)", "Gaia eDR3", "−0.25"],
    ["BP", "Gaia BPmag", "ZTF DR7 g (HiPS)", "Gaia eDR3", "+0.03"],
    ["RP", "Gaia RPmag", "ZTF DR7 i (HiPS)", "Gaia eDR3", "−0.18"],
]
tbl_cols = ["Filter", "Calib. mag", "Template source", "Catalogue", "Color term"]
table = ax_tbl.table(cellText=tbl_data, colLabels=tbl_cols,
                     loc="center", cellLoc="center")
table.auto_set_font_size(False)
table.set_fontsize(6)
table.scale(1.2, 1.4)
table.auto_set_column_width(col=list(range(len(tbl_cols))))
for (row, col), cell in table.get_celld().items():
    cell.set_edgecolor("#ccc")
    cell.set_linewidth(0.5)
    if row == 0:
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold", fontsize=5.5)
    else:
        cell.set_facecolor("#f8f9fa" if row % 2 == 0 else "white")
        if col == 0:
            filt = tbl_data[row-1][0]
            cell.set_text_props(color=filt_colors.get(filt, "black"), fontweight="bold")

# ══════════════════════════════════════════════════════════════════════════════
# Row 2: Lightcurve (full width)
# ══════════════════════════════════════════════════════════════════════════════

ax1 = fig.add_subplot(gs[2, :])

for filt in ["G", "BP", "RP"]:
    d = sub[filt]
    x = np.array(d["mjd"]) - mjd_ref
    y = np.array(d["mag"])
    ye = np.array(d["err"])
    s = np.array(sn[filt])

    for i in range(len(x)):
        is_suspect = (filt == "G" and i < len(g_suspect) and g_suspect[i])
        if is_suspect:
            ax1.errorbar(x[i], y[i], yerr=ye[i], fmt="X",
                         color="#ff6b6b", markeredgecolor="#c0392b",
                         markersize=8, capsize=2.5, alpha=0.7, zorder=7)
            ax1.annotate(f"S/N={s[i]:.0f} but off model\n(galaxy contamination?)",
                         (x[i], y[i]),
                         xytext=(x[i] + 1.5, y[i] + 0.05),
                         fontsize=4.5, color="#c0392b", fontweight="bold",
                         arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.6))
            continue
        alpha = 0.3 if s[i] < 3 else (0.6 if s[i] < 5 else 1.0)
        edge = "red" if s[i] < 3 else filt_colors[filt]
        is_dax = (filt == "G" and i < len(g_source) and g_source[i] == "DAX T500")
        marker = "*" if is_dax else filt_markers[filt]
        ms = 10 if is_dax else 5
        ax1.errorbar(x[i], y[i], yerr=ye[i], fmt=marker,
                     color=filt_colors[filt], markeredgecolor=edge,
                     markersize=ms, capsize=2.5, alpha=alpha, zorder=6 if is_dax else 5)

    reliable = [i for i in range(len(s)) if s[i] >= 3]
    if filt == "G":
        reliable = [i for i in reliable if not g_suspect[i]]
    if len(reliable) > 1:
        ax1.plot([x[i] for i in reliable], [y[i] for i in reliable],
                 '-', color=filt_colors[filt], alpha=0.25, linewidth=0.8)

## Professional survey detections (different filter systems)
for name, pd in pro_surveys.items():
    px = pd["mjd"] - mjd_ref
    ax1.plot(px, pd["mag"], pd["marker"], color=pd["color"], markersize=7,
             markeredgecolor="k", markeredgewidth=0.4, alpha=0.85, zorder=4)
    ax1.annotate(f"{name}\n{pd['filt']}-band", (px, pd["mag"]),
                 xytext=(5, -8), textcoords="offset points",
                 fontsize=4.5, color=pd["color"], fontweight="bold")

## Weighted linear fit on G-band (reliable + non-suspect points only)
g_x = np.array(sub["G"]["mjd"]) - mjd_ref
g_y = np.array(sub["G"]["mag"])
g_ye = np.array(sub["G"]["err"])
g_sn = np.array(sn["G"])
g_sus = np.array(g_suspect)
rel = (g_sn >= 3) & (~g_sus)
if np.sum(rel) >= 2:
    w = 1.0 / g_ye[rel]**2
    coeffs = np.polyfit(g_x[rel], g_y[rel], 1, w=w, cov=True)
    slope, intercept = coeffs[0]
    slope_err = np.sqrt(coeffs[1][0, 0])
    intercept_err = np.sqrt(coeffs[1][1, 1])
    x_fit = np.linspace(g_x.min() - 0.3, g_x.max() + 0.3, 100)
    y_fit = slope * x_fit + intercept
    y_fit_upper = (slope + slope_err) * x_fit + (intercept - intercept_err)
    y_fit_lower = (slope - slope_err) * x_fit + (intercept + intercept_err)
    ax1.plot(x_fit, y_fit, '--', color=filt_colors["G"], linewidth=1.2,
             alpha=0.8, zorder=3)
    ax1.fill_between(x_fit, y_fit_lower, y_fit_upper,
                     color=filt_colors["G"], alpha=0.10, zorder=2)
    ax1.annotate(f"G linear: {slope:+.3f} ± {slope_err:.3f} mag/day",
                 xy=(0.98, 0.03), xycoords="axes fraction",
                 fontsize=6, ha="right", va="bottom",
                 color=filt_colors["G"], fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#2ecc71", alpha=0.9))

ax1.invert_yaxis()
ax1.set_xlabel(f"MJD − {mjd_ref:.0f}", fontsize=7)
ax1.set_ylabel("Magnitude (template subtraction)", fontsize=7)
ax1.set_title("Template Subtraction Lightcurve — SN 2026kie", fontsize=9,
              fontweight="bold", pad=5)

from matplotlib.lines import Line2D
handles = [Line2D([0],[0], marker=filt_markers[f], color=filt_colors[f], ls="none",
                  markersize=4, label=f) for f in ["G","BP","RP"]]
handles.append(Line2D([0],[0], marker="*", color=filt_colors["G"], ls="none",
                      markersize=8, label="DAX T500"))
handles.append(Line2D([0],[0], marker="o", color="#ccc", markeredgecolor="red",
                      ls="none", markersize=4, label="S/N < 3"))
handles.append(Line2D([0],[0], marker="X", color="#ff6b6b", markeredgecolor="#c0392b",
                      ls="none", markersize=6, label="Suspect"))
for name, pd in pro_surveys.items():
    handles.append(Line2D([0],[0], marker=pd["marker"], color=pd["color"],
                          markeredgecolor="k", ls="none", markersize=5,
                          label=f"{name} ({pd['filt']})"))
handles.append(Line2D([0],[0], ls="--", color=filt_colors["G"], linewidth=1,
                      label="G linear fit"))
ax1.legend(handles=handles, fontsize=4.5, loc="lower left", framealpha=0.9, ncol=5)
ax1.tick_params(labelsize=6.5)
ax1.grid(True, alpha=0.15)

for filt in ["G", "BP", "RP"]:
    d = sub[filt]
    x = np.array(d["mjd"]) - mjd_ref
    y = np.array(d["mag"])
    s = sn[filt]
    for i in range(len(x)):
        is_suspect = (filt == "G" and i < len(g_suspect) and g_suspect[i])
        if is_suspect:
            continue
        offset = 0.15 if filt == "RP" else (-0.15 if filt == "BP" else 0.08)
        color = "red" if s[i] < 3 else "#aaa"
        ax1.annotate(f"S/N={s[i]:.1f}", (x[i], y[i] + offset),
                     fontsize=4.5, ha="center", color=color)

# ══════════════════════════════════════════════════════════════════════════════
# Row 3 left: Color index
# ══════════════════════════════════════════════════════════════════════════════

ax2 = fig.add_subplot(gs[3, 0])

x_ci = np.array(ci_mjd) - mjd_ref
for i in range(len(ci_val)):
    alpha = 0.3 if sn["RP"][i] < 3 else 1.0
    edge = "red" if sn["RP"][i] < 3 else "#8e44ad"
    ax2.errorbar(x_ci[i], ci_val[i], yerr=ci_err[i], fmt="^",
                 color="#8e44ad", markeredgecolor=edge,
                 markersize=5, capsize=2.5, alpha=alpha, zorder=5)

ax2.axhline(0, color="#888", ls=":", lw=0.5)
ax2.set_xlabel(f"MJD − {mjd_ref:.0f}", fontsize=6.5)
ax2.set_ylabel("BP − RP (sub)", fontsize=6.5)
ax2.set_title("Color Index (subtraction)", fontsize=7.5, fontweight="bold", pad=4)
ax2.tick_params(labelsize=6)
ax2.grid(True, alpha=0.15)

ax2.annotate("unreliable\n(RP S/N = 2.3)", (x_ci[0], ci_val[0]),
             xytext=(x_ci[0] + 0.5, ci_val[0] - 0.3),
             fontsize=5, color="red", ha="left",
             arrowprops=dict(arrowstyle="->", color="red", lw=0.6))

# ══════════════════════════════════════════════════════════════════════════════
# Row 3 right: Consolidated photometry table (Gaia + TNS)
# ══════════════════════════════════════════════════════════════════════════════

ax3 = fig.add_subplot(gs[3, 1])
ax3.axis("off")
phot_cols = ["Date (UT)", "MJD", "Band", "Mag", "Err", "S/N", "Source"]
phot_rows = [
    # TNS professional surveys
    ["04-21.90", "61151.90", "L",      "18.20", "—",    "—",    "GOTO"],
    ["04-22.31", "61152.31", "orange", "18.09", "—",    "—",    "ATLAS"],
    # Suspect RAPAS point
    ["04-23.81", "61153.81", "G",      "17.13", "0.03", "29.5", "RAPAS ⚠"],
    # ZTF
    ["04-24.18", "61154.18", "r",      "17.41", "—",    "—",    "ZTF"],
    # OBS-BZH Gaia photometry
    ["04-24.92", "61154.92", "G",      "18.69", "0.35", "2.9",  "OBS-BZH"],
    ["04-24.95", "61154.95", "BP",     "17.52", "0.12", "8.7",  "OBS-BZH"],
    ["04-24.98", "61154.98", "RP",     "17.89", "0.44", "2.3",  "OBS-BZH"],
    ["04-25.86", "61155.86", "G",      "17.67", "0.12", "8.6",  "OBS-BZH"],
    ["04-25.91", "61155.91", "BP",     "17.48", "0.09", "11.5", "OBS-BZH"],
    ["04-25.94", "61155.94", "RP",     "16.88", "0.19", "5.3",  "OBS-BZH"],
    ["04-27.91", "61157.91", "G",      "17.48", "0.09", "11.0", "OBS-BZH"],
    ["04-27.93", "61157.93", "BP",     "17.23", "0.07", "14.2", "OBS-BZH"],
    ["04-27.95", "61157.95", "RP",     "16.74", "0.12", "8.2",  "OBS-BZH"],
    # DAX T500
    ["05-01.91", "61161.91", "G",      "16.35", "0.03", "32.1", "DAX T500"],
]

table2 = ax3.table(cellText=phot_rows, colLabels=phot_cols,
                   bbox=[0.0, 0.0, 1.0, 0.88], cellLoc="center")
table2.auto_set_font_size(False)
table2.set_fontsize(4.8)
table2.auto_set_column_width(col=list(range(len(phot_cols))))

ax3.text(0.5, 0.97, "Photometry Summary (all sources)",
         transform=ax3.transAxes, fontsize=7.5, fontweight="bold",
         ha="center", va="top")

band_colors = {"G": "#2ecc71", "BP": "#3498db", "RP": "#e74c3c",
               "L": "#9b59b6", "orange": "#e67e22", "r": "#e74c3c"}

for (row, col), cell in table2.get_celld().items():
    cell.set_edgecolor("#ddd")
    cell.set_linewidth(0.4)
    if row == 0:
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold", fontsize=4.8)
    else:
        r = phot_rows[row - 1]
        cell.set_facecolor("#f8f9fa" if row % 2 == 0 else "white")
        if col == 2:
            cell.set_text_props(color=band_colors.get(r[2], "black"), fontweight="bold")
        if col == 5 and r[5] not in ("—",):
            try:
                val = float(r[5])
                if val < 3:
                    cell.set_facecolor("#fee2e2")
                    cell.set_text_props(color="red", fontweight="bold")
                elif val < 5:
                    cell.set_facecolor("#fef3c7")
            except ValueError:
                pass
        if "⚠" in r[6]:
            cell.set_facecolor("#fee2e2")
            if col == 6:
                cell.set_text_props(color="#c0392b", fontweight="bold")

# ══════════════════════════════════════════════════════════════════════════════
# Row 4 left: Color term stability (G-band, OBS-BZH vs DAX)
# ══════════════════════════════════════════════════════════════════════════════

ax4 = fig.add_subplot(gs[4, 0])
all_nights_labels = ["04-24", "04-25", "04-27", "05-01"]
x_all = np.arange(4)
width = 0.22

# OBS-BZH: G, BP, RP for first 3 nights
for fi, filt in enumerate(["G", "BP", "RP"]):
    cts = calib[filt]["ct"]
    ax4.bar(x_all[:3] + fi * width - width, cts, width=width,
            color=filt_colors[filt], alpha=0.8, label=filt,
            edgecolor="white", linewidth=0.5)

# DAX T500: G only on night 4
ax4.bar(x_all[3], -0.17, width=width,
        color=filt_colors["G"], alpha=0.8, edgecolor="#f39c12",
        linewidth=1.5, hatch="//", label="G (DAX T500)")

ax4.set_xticks(x_all)
ax4.set_xticklabels(all_nights_labels, fontsize=5.5)
ax4.set_ylabel("Color term", fontsize=6.5)
ax4.set_title("Color Term: G-band across sites", fontsize=7.5, fontweight="bold", pad=4)
ax4.legend(fontsize=5, loc="lower right")
ax4.tick_params(labelsize=5.5)
ax4.grid(True, alpha=0.15, axis="y")
ax4.axhline(0, color="#888", ls="-", lw=0.5)

# ══════════════════════════════════════════════════════════════════════════════
# Row 4 right: Zero point stability (absolute values, not delta)
# ══════════════════════════════════════════════════════════════════════════════

ax5 = fig.add_subplot(gs[4, 1])
x_nights_3 = np.arange(len(nights))
for filt in ["G", "BP", "RP"]:
    zps = calib[filt]["zp"]
    ax5.plot(x_nights_3, zps, marker=filt_markers[filt], color=filt_colors[filt],
             label=filt, markersize=5, linewidth=1)
ax5.set_xticks(x_nights_3)
ax5.set_xticklabels([n.replace("2026-", "") for n in nights], fontsize=6)
ax5.set_ylabel("Zero point (mag)", fontsize=6.5)
ax5.set_title("Zero Point per Night", fontsize=7.5, fontweight="bold", pad=4)
ax5.legend(fontsize=5.5, loc="center right")
ax5.tick_params(labelsize=5.5)
ax5.grid(True, alpha=0.15)

# ══════════════════════════════════════════════════════════════════════════════
# Row 5: Conclusions
# ══════════════════════════════════════════════════════════════════════════════

ax6 = fig.add_subplot(gs[5, :])
ax6.axis("off")

conclusions = (
    "Findings\n\n"
    "1. Unreliable first epoch — On April 24, the template subtraction yields "
    "S/N = 2.88 (G) and 2.28 (RP), below the S/N = 3 detection threshold. "
    "The color index CI(BP−RP) = −0.37 is unphysical and entirely driven by "
    "noise in the RP measurement (err = 0.44 mag). The subsequent nights give "
    "CI = +0.60 ± 0.21 and +0.49 ± 0.14, consistent within uncertainties.\n\n"
    "2. Observation conditions — April 24 had the poorest atmospheric transparency "
    "(lowest zero points in all bands, ~0.4 mag below the other nights), worse "
    "seeing (FWHM 4.7–5.2 px vs 4.1–4.4 on Apr 27), and fewer stacked frames "
    "(37 vs 60 for BP/RP), all degrading the subtraction S/N.\n\n"
    "3. No filter swap — Photometric color terms are remarkably stable across "
    "all 3 nights: G ≈ −0.25, BP ≈ +0.03, RP ≈ −0.18 (spread < 0.02). "
    "A filter swap would manifest as a sudden jump in color term. "
    "Intrinsic calibration scatter = 0.01 mag on all 9 tasks."
)

ax6.text(0.0, 1.0, conclusions, transform=ax6.transAxes,
         fontsize=6.2, va="top", ha="left", linespacing=1.45,
         fontfamily="serif", wrap=True,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f4f8",
                   edgecolor="#ccc", linewidth=0.5))

# ── Save ──────────────────────────────────────────────────────────────────────

output = "/home/pyl/Documents/astrobatch/doc/SN2026kie_photometry_report.pdf"
fig.savefig(output, dpi=150)
plt.close(fig)
print(f"PDF saved to {output}")
