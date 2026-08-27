import math
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np


R_CAL = 1.98720425864083
LINESTYLES = ["-", "--", "-.", ":", (0, (6, 2, 1, 2, 1, 2)), (0, (3, 1, 1, 1))]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "<", ">"]
COLORS = [
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
    "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


@dataclass
class Reaction:
    equation: str
    reactants: tuple
    products: tuple
    product_label: str
    plog: dict = field(default_factory=dict)


def is_number_token(token):
    try:
        float(token)
        return True
    except ValueError:
        return token.lower() in {"inf", "+inf", "-inf"}


def extract_equation(raw):
    clean = raw.split("!", 1)[0].strip()
    fields = clean.split()
    first_numeric = next((i for i, token in enumerate(fields) if is_number_token(token)), len(fields))
    return "".join(fields[:first_numeric])


def parse_side(side):
    return tuple(x.strip() for x in side.split("+") if x.strip())


def canonical_side(side):
    return tuple(sorted(side))


def split_equation(eq):
    if "=>" in eq:
        left, right = eq.split("=>", 1)
    elif "=" in eq:
        left, right = eq.split("=", 1)
    else:
        raise ValueError(eq)
    reactants = canonical_side(parse_side(left))
    products_ordered = parse_side(right)
    return reactants, canonical_side(products_ordered), "+".join(products_ordered)


def is_reaction_line(line):
    s = line.strip()
    if not s or s.startswith("!"):
        return False
    if s.upper().startswith(("PLOG", "LOW", "TROE", "SRI", "DUP", "DUPLICATE", "REV", "FORD")):
        return False
    return "=>" in s or "=" in s


def parse_kinetics(path):
    reactions = []
    current = None
    in_reactions = False
    with open(path, "r") as handle:
        for raw in handle:
            line = raw.strip()
            upper = line.upper()
            if upper.startswith("REACTIONS"):
                in_reactions = True
                continue
            if in_reactions and upper == "END":
                break
            if not in_reactions and is_reaction_line(raw):
                in_reactions = True
            if not in_reactions:
                continue
            if is_reaction_line(raw):
                equation = extract_equation(raw)
                try:
                    reactants, products, product_label = split_equation(equation)
                except ValueError:
                    current = None
                    continue
                current = Reaction(equation, reactants, products, product_label)
                reactions.append(current)
                continue
            if current and line.upper().startswith("PLOG") and not line.startswith("!"):
                clean = line.split("!", 1)[0]
                nums = re.findall(
                    r"[-+]?inf|[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?",
                    clean,
                    flags=re.I,
                )
                if len(nums) >= 4:
                    try:
                        pressure, A, beta, Ea = float(nums[0]), float(nums[1]), float(nums[2]), float(nums[3])
                    except ValueError:
                        continue
                    if math.isfinite(A) and A > 0:
                        current.plog.setdefault(pressure, []).append((A, beta, Ea))
    return reactions


def k_arrhenius(params, T):
    A, beta, Ea = params
    return A * (T ** beta) * math.exp(-Ea / (R_CAL * T))


def k_plog(params_list, T):
    return sum(k_arrhenius(params, T) for params in params_list)


def species_label(species):
    return "+".join(species)


def unit_for_reactants(reactants):
    return "cm3/mol/s" if len(reactants) == 2 else "1/s"


def pressure_style_map(pressures):
    ordered = sorted(pressures, reverse=True)
    return {p: LINESTYLES[i % len(LINESTYLES)] for i, p in enumerate(ordered)}


def pressure_marker_map(pressures):
    ordered = sorted(pressures, reverse=True)
    return {p: MARKERS[i % len(MARKERS)] for i, p in enumerate(ordered)}


def force_log_minor_ticks(ax, FixedLocator, NullFormatter):
    ymin, ymax = ax.get_ylim()
    if ymin <= 0 or ymax <= 0:
        return
    lo = math.floor(math.log10(ymin)) - 1
    hi = math.ceil(math.log10(ymax)) + 1
    ticks = []
    for exp in range(lo, hi + 1):
        decade = 10.0 ** exp
        for mult in range(2, 10):
            value = mult * decade
            if ymin < value < ymax:
                ticks.append(value)
    ax.yaxis.set_minor_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="y", which="minor", left=True, right=False, length=3.8, width=0.8)
    ax.tick_params(axis="y", which="major", left=True, right=False, length=6.0, width=1.0)


def set_bold_axes(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.set_title(title, fontweight="bold")
    for label in ax.get_xticklabels() + ax.get_yticklabels(minor=False) + ax.get_yticklabels(minor=True):
        label.set_fontweight("bold")


def parse_pressure_from_folder(folder_name):
    if not folder_name.endswith("atm"):
        return None
    try:
        return float(folder_name[:-3])
    except ValueError:
        return None


def read_rates_file(path):
    header = None
    data = []
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header = stripped.lstrip("#").split()
                continue
            parts = stripped.split()
            try:
                data.append([float(x) for x in parts])
            except ValueError:
                continue
    if header is None or len(header) < 2 or len(data) == 0:
        return None
    arr = np.array(data, dtype=float)
    if arr.shape[1] != len(header):
        ncols = min(arr.shape[1], len(header))
        arr = arr[:, :ncols]
        header = header[:ncols]
    return header, arr


def load_tabulated_rates(cwd):
    lumping_dir = os.path.join(cwd, "lumping")
    tabulated = defaultdict(lambda: defaultdict(dict))
    if not os.path.isdir(lumping_dir):
        return tabulated
    for reactant in sorted(os.listdir(lumping_dir)):
        reactant_dir = os.path.join(lumping_dir, reactant)
        if not os.path.isdir(reactant_dir):
            continue
        for pressure_folder in sorted(os.listdir(reactant_dir)):
            pressure = parse_pressure_from_folder(pressure_folder)
            if pressure is None:
                continue
            rates_path = os.path.join(reactant_dir, pressure_folder, "rates.txt")
            if not os.path.isfile(rates_path):
                continue
            parsed = read_rates_file(rates_path)
            if parsed is None:
                continue
            header, arr = parsed
            T = arr[:, 0]
            for idx, name in enumerate(header[1:], start=1):
                tabulated[reactant][pressure][name] = (T, arr[:, idx])
    return tabulated


def tabulated_for_reaction(tabulated, reaction, pressure):
    matching_reactants = [reactant for reactant in reaction.reactants if reactant in tabulated]
    if not matching_reactants:
        return None
    product_names = list(reaction.products)
    product_label = reaction.product_label
    for reactant in matching_reactants:
        pressure_data = tabulated[reactant].get(pressure)
        if not pressure_data:
            continue
        if product_label in pressure_data:
            return pressure_data[product_label]
        matching_products = [name for name in product_names if name in pressure_data]
        if matching_products:
            T = pressure_data[matching_products[0]][0]
            y = np.zeros_like(T)
            for name in matching_products:
                y = y + pressure_data[name][1]
            return T, y
    return None


def product_colors(reactions):
    labels = sorted({reaction.product_label for reaction in reactions if reaction.plog})
    return {label: COLORS[i % len(COLORS)] for i, label in enumerate(labels)}


def group_channels(reactions):
    groups = defaultdict(list)
    for reaction in reactions:
        if reaction.plog:
            groups[reaction.reactants].append(reaction)
    return {reactants: channels for reactants, channels in groups.items() if len(channels) >= 2}


def plot_all(cwd, tmin=500.0, tmax=2500.0, nT=250, min_bf=1e-3):
    kin_path = os.path.join(cwd, "lumpedmech", "kin.txt")
    if not os.path.isfile(kin_path):
        raise FileNotFoundError("plot_rates requires lumpedmech/kin.txt, not found: " + kin_path)

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FixedLocator, NullFormatter

    outdir = os.path.join(cwd, "lumpedmech_plots")
    os.makedirs(outdir, exist_ok=True)
    rates_pdf = os.path.join(outdir, "rates.pdf")
    bf_pdf = os.path.join(outdir, "branching_fractions.pdf")

    reactions = parse_kinetics(kin_path)
    reactions_plog = [reaction for reaction in reactions if reaction.plog]
    tabulated = load_tabulated_rates(cwd)
    temperatures = np.linspace(tmin, tmax, nT)
    xfit = 1000.0 / temperatures

    with PdfPages(rates_pdf) as pdf:
        for reaction in reactions_plog:
            pressures = sorted(reaction.plog)
            styles = pressure_style_map(pressures)
            markers = pressure_marker_map(pressures)
            fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
            pressure_handles = []
            marker_handles = []
            for pressure in sorted(pressures, reverse=True):
                yfit = [k_plog(reaction.plog[pressure], T) for T in temperatures]
                ax.plot(xfit, yfit, color="tab:blue", linestyle=styles[pressure], linewidth=2.0)
                pressure_handles.append(Line2D([0], [0], color="tab:blue", linestyle=styles[pressure],
                                               linewidth=2.0, label=f"{pressure:g} atm fit"))
                tab = tabulated_for_reaction(tabulated, reaction, pressure)
                if tab is not None:
                    Ttab, ytab = tab
                    mask = np.isfinite(Ttab) & np.isfinite(ytab) & (ytab > 0)
                    if np.any(mask):
                        ax.scatter(
                            1000.0 / Ttab[mask],
                            ytab[mask],
                            marker=markers[pressure],
                            facecolors="none",
                            edgecolors="black",
                            s=32,
                            linewidths=0.9,
                            zorder=4,
                        )
                        marker_handles.append(Line2D([0], [0], color="black", marker=markers[pressure],
                                                     linestyle="None", markerfacecolor="none",
                                                     label=f"{pressure:g} atm data"))
            ax.set_yscale("log")
            force_log_minor_ticks(ax, FixedLocator, NullFormatter)
            set_bold_axes(ax, "1000/T [1/K]", f"k [{unit_for_reactants(reaction.reactants)}]", reaction.equation)
            handles = pressure_handles + marker_handles
            labels_seen = set()
            handles_unique = []
            for handle in handles:
                label = handle.get_label()
                if label not in labels_seen:
                    handles_unique.append(handle)
                    labels_seen.add(label)
            legend = ax.legend(handles=handles_unique, fontsize=8.5, loc="best", framealpha=0.82)
            for text in legend.get_texts():
                text.set_fontweight("bold")
            ax.grid(True, which="major", alpha=0.35)
            ax.grid(True, which="minor", alpha=0.15)
            pdf.savefig(fig)
            plt.close(fig)

    colors = product_colors(reactions_plog)
    groups = group_channels(reactions_plog)
    with PdfPages(bf_pdf) as pdf:
        for reactants, channels in sorted(groups.items(), key=lambda item: species_label(item[0])):
            pressures = sorted(set.union(*(set(channel.plog) for channel in channels)))
            styles = pressure_style_map(pressures)
            fig, ax = plt.subplots(figsize=(10.8, 5.8), constrained_layout=True)
            visible_products = set()
            for pressure in sorted(pressures, reverse=True):
                rate_by_index = {}
                for idx, reaction in enumerate(channels):
                    if pressure in reaction.plog:
                        params = reaction.plog[pressure]
                    else:
                        nearest_pressure = min(reaction.plog, key=lambda p: abs(math.log(p / pressure)))
                        params = reaction.plog[nearest_pressure]
                    rate_by_index[idx] = np.array([k_plog(params, T) for T in temperatures])
                totals = np.sum(np.vstack(list(rate_by_index.values())), axis=0)
                for idx, reaction in enumerate(channels):
                    bf = np.divide(
                        rate_by_index[idx],
                        totals,
                        out=np.full_like(totals, min_bf),
                        where=totals > 0,
                    )
                    if np.nanmax(bf) < min_bf and reaction.product_label not in visible_products:
                        continue
                    ax.plot(
                        xfit,
                        bf,
                        color=colors[reaction.product_label],
                        linestyle=styles[pressure],
                        linewidth=1.8,
                    )
                    visible_products.add(reaction.product_label)
            if not visible_products:
                plt.close(fig)
                continue
            ax.set_yscale("log")
            ax.set_ylim(min_bf, 1.0)
            force_log_minor_ticks(ax, FixedLocator, NullFormatter)
            set_bold_axes(ax, "1000/T [1/K]", "BF [-]", f"{species_label(reactants)} branching")
            color_handles = [
                Line2D([0], [0], color=colors[label], linestyle="-", linewidth=2.4, label=label)
                for label in sorted(visible_products)
            ]
            style_handles = [
                Line2D([0], [0], color="black", linestyle=styles[p], linewidth=2.2, label=f"{p:g} atm")
                for p in sorted(pressures, reverse=True)
            ]
            leg1 = ax.legend(handles=color_handles, title="Product", fontsize=9.0, title_fontsize=10.0,
                             loc="upper right", framealpha=0.82)
            leg2 = ax.legend(handles=style_handles, title="Pressure", fontsize=9.0, title_fontsize=10.0,
                             loc="lower right", framealpha=0.82)
            ax.add_artist(leg1)
            for legend in (leg1, leg2):
                legend.get_title().set_fontweight("bold")
                for text in legend.get_texts():
                    text.set_fontweight("bold")
            ax.grid(True, which="major", alpha=0.35)
            ax.grid(True, which="minor", alpha=0.18)
            pdf.savefig(fig)
            plt.close(fig)

    return {"rates_pdf": rates_pdf, "bf_pdf": bf_pdf, "n_reactions": len(reactions_plog), "n_bf_groups": len(groups)}
