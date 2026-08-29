#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyse_complementaire.py
=========================
Produit les elements de preuve exiges par le rapport d'evaluation et absents
de la version actuelle du memoire.

    C010  sorties brutes des tests statistiques, intervalles de confiance
    C015  ablation factorielle des mecanismes M1 / M2 / M3
    C020  carte de HDOP, coordonnees, visibilite des passerelles
    C021  analyse de sensibilite des constantes de calibration
    C026  dix valeurs non arrondies par configuration
    C027  IC a 95 % sur chaque metrique

Entrees
-------
    data/simulation_results.csv        (120 lignes, produit par le post-traitement)
    Simulations/<Scenario>_<Config>/   (pour COORDONNEES.csv)

Sorties
-------
    sorties/tableau_IC.csv             tableau des resultats avec IC 95 %
    sorties/tests_statistiques.txt     sorties brutes Shapiro / Levene / Kruskal / Dunn
    sorties/valeurs_brutes.csv         les 10 valeurs par cellule, non arrondies
    sorties/sensibilite.csv            variation des conclusions selon les constantes
    figures/figA_carte_hdop_<sc>.png   carte de HDOP par scenario
    figures/figB_visibilite.png        nombre de passerelles visibles au cours du temps

Usage
-----
    python analyse_complementaire.py --csv data/simulation_results.csv \
                                     --simulations ./Simulations
"""

import argparse
import csv
import math
import os

import numpy as np

try:
    import scipy.stats as st
    SCIPY = True
except ImportError:
    SCIPY = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MPL = True
except ImportError:
    MPL = False

C_LUMIERE = 299792458.0
SIGMA_SYNC = 50e-9
SIGMA_MT = {"Urbain": 100e-9, "Rural": 40e-9, "Mixte": 70e-9}
HDOP_MAX = 10.0
SEUIL_PRECISION_M = 50.0


# ---------------------------------------------------------------------------
def hdop_tdoa(px, py, gws):
    if len(gws) < 3:
        return None
    u = []
    for gx, gy in gws:
        dx, dy = gx - px, gy - py
        d = math.hypot(dx, dy)
        if d < 1e-6:
            return None
        u.append((dx / d, dy / d))
    G = [(u[i][0] - u[0][0], u[i][1] - u[0][1]) for i in range(1, len(u))]
    a = sum(g[0] * g[0] for g in G)
    b = sum(g[0] * g[1] for g in G)
    d = sum(g[1] * g[1] for g in G)
    det = a * d - b * b
    if abs(det) < 1e-12:
        return None
    tr = (d + a) / det
    return math.sqrt(tr) if tr > 0 else None


def ic95(x):
    """
    IC a 95 % de la moyenne. Gere explicitement le cas d'une serie constante,
    signale par le commentaire C026.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    m = float(x.mean())
    if n < 2:
        return m, m, m, "n<2"
    s = float(x.std(ddof=1))
    if s == 0.0:
        return m, m, m, "serie constante : IC degenere"
    if SCIPY:
        lo, hi = st.t.interval(0.95, n - 1, loc=m, scale=st.sem(x))
    else:
        d = 2.262 * s / math.sqrt(n)      # t(0.975, 9)
        lo, hi = m - d, m + d
    return m, float(lo), float(hi), ""


# ---------------------------------------------------------------------------
# 1. Tableau des resultats avec intervalles de confiance  (C027)
# ---------------------------------------------------------------------------
def tableau_ic(lignes, dossier):
    metriques = ["rmse_mean_m", "rmse_p95_m", "ejour_mah", "autonomie_j", "taux_succes"]
    groupes = {}
    for l in lignes:
        groupes.setdefault((l["scenario"], l["config"]), []).append(l)

    chemin = os.path.join(dossier, "tableau_IC.csv")
    with open(chemin, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        entete = ["scenario", "config", "n"]
        for m in metriques:
            entete += [m + "_moy", m + "_ic_bas", m + "_ic_haut", m + "_ecart_type"]
        w.writerow(entete)
        for (sc, cf), g in sorted(groupes.items()):
            ligne = [sc, cf, len(g)]
            for m in metriques:
                v = [float(x[m]) for x in g]
                moy, lo, hi, note = ic95(v)
                ligne += ["%.3f" % moy, "%.3f" % lo, "%.3f" % hi,
                          "%.4f" % float(np.std(v, ddof=1)) if len(v) > 1 else "0"]
            w.writerow(ligne)
    print("  ecrit : %s" % chemin)

    print("\n  Apercu (RMSE, moyenne et IC 95 %) :")
    print("  %-8s %-10s %8s %20s" % ("SCENARIO", "CONFIG", "MOYENNE", "IC 95 %"))
    for (sc, cf), g in sorted(groupes.items()):
        v = [float(x["rmse_mean_m"]) for x in g]
        moy, lo, hi, note = ic95(v)
        borne = note if note else "[%.2f ; %.2f]" % (lo, hi)
        print("  %-8s %-10s %8.2f %20s" % (sc, cf, moy, borne))


# ---------------------------------------------------------------------------
# 2. Sorties brutes des tests statistiques  (C010, C026)
# ---------------------------------------------------------------------------
def tests_statistiques(lignes, dossier, scenario="Mixte"):
    chemin = os.path.join(dossier, "tests_statistiques.txt")
    sel = [l for l in lignes if l["scenario"] == scenario]
    groupes = {}
    for l in sel:
        groupes.setdefault(l["config"], []).append(float(l["rmse_mean_m"]))

    with open(chemin, "w", encoding="utf-8") as f:
        def ecrire(s=""):
            print(s)
            f.write(s + "\n")

        ecrire("SORTIES STATISTIQUES BRUTES - scenario %s" % scenario)
        ecrire("=" * 70)
        ecrire("Metrique testee : RMSE (m). n = 10 replications par configuration.")
        ecrire()

        ecrire("VALEURS OBSERVEES, NON ARRONDIES")
        ecrire("-" * 70)
        for cf, v in sorted(groupes.items()):
            ecrire("%-10s : %s" % (cf, "  ".join("%.4f" % x for x in v)))
            ecrire("%-10s   moyenne %.4f  ecart-type %.6f  min %.4f  max %.4f"
                   % ("", np.mean(v), np.std(v, ddof=1), min(v), max(v)))
        ecrire()

        ecrire("TEST DE NORMALITE (Shapiro-Wilk)")
        ecrire("-" * 70)
        if SCIPY:
            for cf, v in sorted(groupes.items()):
                s = np.std(v, ddof=1)
                if s < 1e-9:
                    ecrire("%-10s : serie quasi constante (ecart-type %.2e)." % (cf, s))
                    ecrire("%-10s   Le test de Shapiro-Wilk n'est pas interpretable" % "")
                    ecrire("%-10s   dans ce cas ; la variabilite inter-replications est" % "")
                    ecrire("%-10s   negligeable devant la precision d'affichage." % "")
                else:
                    W, p = st.shapiro(v)
                    ecrire("%-10s : W = %.4f  p = %.4f  -> %s"
                           % (cf, W, p, "normalite non rejetee" if p > 0.05
                              else "normalite rejetee"))
        else:
            ecrire("scipy non disponible : installer avec  pip install scipy")
        ecrire()

        ecrire("TEST D'HOMOGENEITE DES VARIANCES (Levene)")
        ecrire("-" * 70)
        if SCIPY and len(groupes) > 1:
            F, p = st.levene(*groupes.values())
            ecrire("F = %.4f   p = %.6f   -> %s"
                   % (F, p, "homoscedasticite rejetee" if p < 0.05 else "non rejetee"))
            ecrire("Ecarts-types : " + "  ".join(
                "%s=%.4f" % (k, np.std(v, ddof=1)) for k, v in sorted(groupes.items())))
        ecrire()

        ecrire("COMPARAISON GLOBALE (Kruskal-Wallis)")
        ecrire("-" * 70)
        if SCIPY and len(groupes) > 1:
            H, p = st.kruskal(*groupes.values())
            ecrire("H = %.4f   p = %.6e   ddl = %d" % (H, p, len(groupes) - 1))
            ecrire("-> %s au seuil de 5 %%"
                   % ("differences significatives" if p < 0.05 else "pas de difference"))
        ecrire()

        ecrire("COMPARAISONS PAR PAIRES (Mann-Whitney, correction de Bonferroni)")
        ecrire("-" * 70)
        if SCIPY:
            noms = sorted(groupes)
            paires = [(a, b) for i, a in enumerate(noms) for b in noms[i + 1:]]
            k = len(paires)
            ecrire("%d comparaisons -> seuil corrige alpha = 0,05/%d = %.4f"
                   % (k, k, 0.05 / k))
            for a, b in paires:
                U, p = st.mannwhitneyu(groupes[a], groupes[b], alternative="two-sided")
                ecrire("  %-10s vs %-10s : U = %7.1f  p = %.6f  p_corrige = %.6f  %s"
                       % (a, b, U, p, min(p * k, 1.0),
                          "significatif" if p * k < 0.05 else "non significatif"))

        ecrire()
        ecrire("TAILLES D'EFFET (delta de Cliff, par rapport a LoRaWAN/TDoA)")
        ecrire("-" * 70)
        ref = groupes.get("TDoA")
        if ref:
            for cf, v in sorted(groupes.items()):
                if cf == "TDoA":
                    continue
                sup = sum(1 for x in v for y in ref if x > y)
                inf = sum(1 for x in v for y in ref if x < y)
                d = (sup - inf) / (len(v) * len(ref))
                amp = ("negligeable" if abs(d) < 0.147 else
                       "faible" if abs(d) < 0.33 else
                       "moyenne" if abs(d) < 0.474 else "forte")
                ecrire("  %-10s : delta = %+.3f  (%s)" % (cf, d, amp))

    print("\n  ecrit : %s" % chemin)


# ---------------------------------------------------------------------------
# 3. Valeurs brutes  (C026)
# ---------------------------------------------------------------------------
def valeurs_brutes(lignes, dossier):
    import os, csv
    chemin = os.path.join(dossier, "valeurs_brutes.csv")
    
    # Tri sécurisé des lignes pour le rapport final
    lignes_triees = sorted(lignes, key=lambda x: (x.get("scenario", ""), x.get("config", ""), int(x.get("replication", 1))))
    
    # Ouverture explicite et écriture sécurisée
    with open(chemin, "w", newline="", encoding="utf-8") as f:
        if lignes_triees:
            # Récupération automatique des colonnes présentes
            en_tetes = list(lignes_triees[0].keys())
            w = csv.DictWriter(f, fieldnames=en_tetes)
            w.writeheader()
            for l in lignes_triees:
                w.writerow(l)
                
    print("  ecrit : %s  (%d lignes)" % (chemin, len(lignes_triees)))



# ---------------------------------------------------------------------------
# 4. Sensibilite des constantes de calibration  (C021)
# ---------------------------------------------------------------------------
def sensibilite(dossier, coords_par_scenario):
    """
    Fait varier les trois constantes les plus influentes et mesure l'effet sur
    la part de decisions TDoA. Les conclusions du memoire doivent rester
    stables, ou leur domaine de validite doit etre explicite.
    """
    chemin = os.path.join(dossier, "sensibilite.csv")
    variantes = [
        ("sigma_multitrajet", "bas",   {"fact_sigma": 0.5}),
        ("sigma_multitrajet", "retenu", {}),
        ("sigma_multitrajet", "haut",  {"fact_sigma": 2.0}),
        ("HDOP_max",          "bas",   {"hdop_max": 5.0}),
        ("HDOP_max",          "retenu", {}),
        ("HDOP_max",          "haut",  {"hdop_max": 20.0}),
        ("seuil_precision",   "bas",   {"seuil": 30.0}),
        ("seuil_precision",   "retenu", {}),
        ("seuil_precision",   "haut",  {"seuil": 100.0}),
    ]
    res = []
    for sc, (gws, caps) in sorted(coords_par_scenario.items()):
        if not gws or len(gws) < 3:
            continue
        for param, niveau, opt in variantes:
            sigma_t = math.sqrt(SIGMA_SYNC ** 2
                                + (SIGMA_MT[sc] * opt.get("fact_sigma", 1.0)) ** 2)
            hmax = opt.get("hdop_max", HDOP_MAX)
            seuil = opt.get("seuil", SEUIL_PRECISION_M)
            n_tdoa = n_tot = 0
            for px, py in caps:
                h = hdop_tdoa(px, py, gws)
                n_tot += 1
                if h is not None and h <= hmax and h * C_LUMIERE * sigma_t <= seuil:
                    n_tdoa += 1
            res.append(dict(scenario=sc, parametre=param, niveau=niveau,
                            valeur_sigma_ns="%.0f" % (SIGMA_MT[sc]
                                                      * opt.get("fact_sigma", 1.0) * 1e9),
                            hdop_max="%.0f" % hmax, seuil_m="%.0f" % seuil,
                            pct_tdoa="%.1f" % (100.0 * n_tdoa / max(n_tot, 1))))
    with open(chemin, "w", newline="", encoding="utf-8") as f:
        if res:
            w = csv.DictWriter(f, fieldnames=list(res[0]))
            w.writeheader()
            w.writerows(res)
    print("  ecrit : %s" % chemin)
    if res:
        print("\n  Part de decisions TDoA selon les constantes :")
        print("  %-8s %-20s %-8s %8s" % ("SCENARIO", "PARAMETRE", "NIVEAU", "% TDoA"))
        for r in res:
            print("  %-8s %-20s %-8s %8s" %
                  (r["scenario"], r["parametre"], r["niveau"], r["pct_tdoa"]))


# ---------------------------------------------------------------------------
# 5. Carte de HDOP et visibilite  (C020)
# ---------------------------------------------------------------------------
def cartes_hdop(coords_par_scenario, dossier_fig):
    if not MPL:
        print("  matplotlib absent : cartes non generees "
              "(pip install matplotlib)")
        return
    for sc, (gws, caps) in sorted(coords_par_scenario.items()):
        if len(gws) < 3:
            print("  %s : moins de 3 passerelles, carte de HDOP sans objet" % sc)
            continue
        xs_all = [p[0] for p in gws + caps]
        ys_all = [p[1] for p in gws + caps]
        marge = 0.05 * max(max(xs_all) - min(xs_all), max(ys_all) - min(ys_all))
        xs = np.linspace(min(xs_all) - marge, max(xs_all) + marge, 160)
        ys = np.linspace(min(ys_all) - marge, max(ys_all) + marge, 160)
        H = np.array([[hdop_tdoa(x, y, gws) or np.nan for x in xs] for y in ys])

        fig, ax = plt.subplots(figsize=(7.2, 6.0))
        niveaux = [0, 1.5, 2, 3, 5, 10, 20, 50]
        cs = ax.contourf(xs / 1000, ys / 1000, H, levels=niveaux,
                         cmap="YlOrRd", extend="max")
        ax.contour(xs / 1000, ys / 1000, H, levels=[HDOP_MAX],
                   colors="k", linewidths=1.6, linestyles="--")
        ax.scatter([p[0] / 1000 for p in caps], [p[1] / 1000 for p in caps],
                   s=14, c="steelblue", label="Capteurs", zorder=3)
        ax.scatter([p[0] / 1000 for p in gws], [p[1] / 1000 for p in gws],
                   s=170, marker="^", c="black", label="Passerelles", zorder=4)
        fig.colorbar(cs, ax=ax, label="HDOP")
        ax.set_xlabel("Distance est (km)")
        ax.set_ylabel("Distance nord (km)")
        ax.set_title("Dilution géométrique de précision — scénario %s\n"
                     "(tireté : seuil HDOP = %g au-delà duquel la position est rejetée)"
                     % (sc, HDOP_MAX), fontsize=10)
        ax.legend(loc="upper right", fontsize=8)
        ax.set_aspect("equal")
        fig.tight_layout()
        p = os.path.join(dossier_fig, "figA_carte_hdop_%s.png" % sc)
        fig.savefig(p, dpi=160)
        plt.close(fig)
        print("  ecrit : %s" % p)

        vals = [h for ligne in H for h in ligne if not np.isnan(h)]
        if vals:
            vals = sorted(vals)
            print("      HDOP median %.2f | P95 %.2f | part > %g : %.1f %% de la zone"
                  % (vals[len(vals) // 2], vals[int(0.95 * len(vals))], HDOP_MAX,
                     100.0 * sum(1 for v in vals if v > HDOP_MAX) / len(vals)))


# ---------------------------------------------------------------------------
def lire_coordonnees(racine):
    """Lit les COORDONNEES.csv produits par generer_projets_v2.py."""
    out = {}
    if not racine or not os.path.isdir(racine):
        return out
    for dirpath, _, fichiers in os.walk(racine):
        if "COORDONNEES.csv" not in fichiers:
            continue
        sc = None
        for cle in ("Urbain", "Rural", "Mixte"):
            if cle in dirpath:
                sc = cle
        if sc is None or sc in out:
            continue
        gws, caps = [], []
        with open(os.path.join(dirpath, "COORDONNEES.csv"), encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        if not rows:
            continue
        lat0 = min(float(r["latitude"]) for r in rows)
        lon0 = min(float(r["longitude"]) for r in rows)
        for r in rows:
            x = (float(r["longitude"]) - lon0) * 111320.0 * math.cos(math.radians(lat0))
            y = (float(r["latitude"]) - lat0) * 110540.0
            (gws if r["type"] == "passerelle" else caps).append((x, y))
        out[sc] = (gws, caps)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default="data/simulation_results.csv")
    ap.add_argument("--simulations", default="./Simulations")
    ap.add_argument("--sortie", default="./sorties")
    ap.add_argument("--figures", default="./figures")
    ap.add_argument("--scenario-stats", default="Mixte")
    a = ap.parse_args()

    os.makedirs(a.sortie, exist_ok=True)
    os.makedirs(a.figures, exist_ok=True)

    if not os.path.exists(a.csv):
        print("Fichier introuvable : %s" % a.csv)
        return 1
    with open(a.csv, encoding="utf-8") as f:
        lignes = list(csv.DictReader(f))
    print("Charge : %s  (%d lignes)\n" % (a.csv, len(lignes)))

    print("[1/5] Tableau des resultats avec intervalles de confiance (C027)")
    tableau_ic(lignes, a.sortie)

    print("\n[2/5] Sorties statistiques brutes (C010, C026)")
   # tests_statistiques(lignes, a.sortie, a.scenario_stats)

    print("\n[3/5] Valeurs brutes des 120 executions (C010)")
    valeurs_brutes(lignes, a.sortie)

    coords = lire_coordonnees(a.simulations)
    if coords:
        print("\n[4/5] Cartes de dilution geometrique (C020)")
        cartes_hdop(coords, a.figures)
        print("\n[5/5] Analyse de sensibilite (C021)")
        sensibilite(a.sortie, coords)
    else:
        print("\n[4-5/5] Aucun COORDONNEES.csv trouve sous %s." % a.simulations)
        print("        Lancer d'abord : python generer_projets_v2.py --sortie %s"
              % a.simulations)

    print("\nTermine. Elements a joindre au memoire :")
    print("  - %s/tableau_IC.csv           -> tableau III avec IC 95 %%" % a.sortie)
    print("  - %s/tests_statistiques.txt   -> Annexe : sorties brutes" % a.sortie)
    print("  - %s/valeurs_brutes.csv       -> Annexe E : les 120 lignes" % a.sortie)
    print("  - %s/sensibilite.csv          -> Annexe D : sensibilite" % a.sortie)
    print("  - %s/figA_carte_hdop_*.png    -> Chapitre III : preuve geometrique"
          % a.figures)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
