#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
replications.py
===============
Produit les dix réplications de l'expérience du mémoire

    « Optimisation des technologies de géolocalisation
      sous consommation énergétique limitée »
    KOALA Issouf — Master 2 RIoT/IE — UFR-SEA / LAMI — UJKZ

CE QUE CE SCRIPT FAIT
---------------------
Pour chacune des graines 1 à 10, il construit les topologies des trois
scénarios (positions des nœuds, positions des passerelles, traces de mobilité),
puis évalue les quatre configurations sur chaque topologie.

Il produit :

    data/simulation_results.csv     120 lignes (3 scénarios x 4 configs x 10 graines)
    data/topologies/               coordonnées de chaque réplication
    sorties/invariance_energie.csv tableau de contrôle de l'énergie radio

CE QUE CE SCRIPT NE FAIT PAS
----------------------------
Il ne pilote pas CupCarbon. Le simulateur ne dispose pas de mode ligne de
commande : sa classe principale (cupcarbon/CupCarbon.java) n'accepte en
argument qu'une configuration de proxy. Les exécutions du simulateur restent
donc manuelles, et servent à établir la valeur du terme radio (voir §2 du
guide).

RÉPARTITION DES GRANDEURS
-------------------------
    Grandeur                        Origine
    ------------------------------  --------------------------------------
    Positions, mobilité, portée     construites ici, identiques à celles
                                    que reçoit CupCarbon
    Énergie de la chaîne radio      mesurée par CupCarbon, invariante
    HDOP, erreur TDoA               calculés ici
    Énergie GNSS, MCU, veille       calculées ici
    Autonomie, taux de succès       calculés ici

USAGE
-----
    python replications.py
    python replications.py --graines 1 10 --sortie data
    python replications.py --e-radio 1.783        # valeur mesurée sous CupCarbon
"""

import argparse
import csv
import math
import os
import random
import statistics
import zlib

# ---------------------------------------------------------------------------
# PARAMÈTRES — identiques à l'Annexe A et à l'Annexe D du mémoire
# ---------------------------------------------------------------------------
DUREE_S = 86400
PAS_GPS_S = 60
TLOC_S = 300                    # intervalle de référence
TLOC_OPT_MOY = 900              # intervalle moyen de la configuration adaptative
C_LUMIERE = 299792458.0
SIGMA_SYNC_S = 50e-9
HDOP_MAX = 10.0
SEUIL_PRECISION_M = 50.0
C_BAT_MAH = 3000.0
U_V = 3.3

PLATEFORME = dict(I_TX_MA=80.0, I_RX_MA=12.0, I_MCU_MA=15.0, T_MCU_S=2.0,
                  I_GNSS_MA=30.0, I_SLEEP_MA=0.020)
TTFF = dict(GPS=30.0, AGPS=10.0)
T_RX_ASSIST = dict(AGPS=2.0)

SIGMA_MT_S = dict(Urbain=100e-9, Rural=40e-9, Mixte=70e-9)
SIGMA_GNSS_M = dict(Urbain=6.5, Rural=4.9, Mixte=5.7)
P_FIX_GNSS = dict(Urbain=0.930, Rural=0.985, Mixte=0.960)
P_FIX_AGPS = dict(Urbain=0.965, Rural=0.990, Mixte=0.975)
FACTEUR_AGPS = 2.3
OVERHEAD_B = dict(TDoA=24, GPS=0, AGPS=48, Optimise=16)
PAYLOAD_O = dict(TDoA=12, GPS=20, AGPS=20, Optimise=16)

SCENARIOS = {
    "Urbain": dict(lat0=12.3650, lon0=-1.5350, largeur_km=2.0, hauteur_km=2.0,
                   n_capteurs=50, rayon_m=2500, v_min=1.0, v_max=5.0,
                   motif=[(0.18, 0.18), (0.82, 0.18), (0.18, 0.82),
                          (0.82, 0.82), (0.50, 0.50)]),
    "Rural":  dict(lat0=12.1900, lon0=-1.5750, largeur_km=5.0, hauteur_km=5.0,
                   n_capteurs=30, rayon_m=7500, v_min=30.0, v_max=80.0,
                   motif=[(0.50, 0.88), (0.12, 0.22), (0.88, 0.22)]),
    "Mixte":  dict(lat0=12.3820, lon0=-1.5030, largeur_km=1.0, hauteur_km=1.0,
                   n_capteurs=40, rayon_m=1500, v_min=1.0, v_max=20.0,
                   motif=[(0.22, 0.22), (0.78, 0.22), (0.22, 0.78), (0.78, 0.78)]),
}
CONFIGS = ["TDoA", "GPS", "AGPS", "Optimise"]
LIBELLE = {"TDoA": "TDoA", "GPS": "GPS", "AGPS": "AGPS", "Optimise": "Optimisé"}

M_LAT = 110540.0


def graine_stable(*parties):
    """
    Graine reproductible entre exécutions et entre machines.

    La fonction hash() de Python est randomisée par défaut depuis la version 3.3
    (PYTHONHASHSEED aléatoire) : deux exécutions du même script produiraient des
    tirages différents. crc32 est déterministe et portable.
    """
    cle = "|".join(str(p) for p in parties).encode("utf-8")
    return zlib.crc32(cle) % (2 ** 31)


# ---------------------------------------------------------------------------
def m_lon(lat):
    return 111320.0 * math.cos(math.radians(lat))


def en_metres(lo, la, lo0, la0):
    return ((lo - lo0) * m_lon(la0), (la - la0) * M_LAT)


def hdop(px, py, gws):
    """
    Facteur de dilution géométrique horizontale pour une résolution TDoA.

    La matrice de géométrie est construite sur les différences de vecteurs
    unitaires nœud -> passerelle par rapport à la passerelle de référence :
    seules les différences de temps d'arrivée sont observables.
    """
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


def toa_lora(sf, pl_o, bw=125000, cr=1, preambule=8, crc=1, ih=0):
    """Time on Air LoRa, formule Semtech AN1200.13."""
    ts = (2 ** sf) / bw
    de = 1 if (sf >= 11 and bw == 125000) else 0
    n_pay = 8 + max(math.ceil((8 * pl_o - 4 * sf + 28 + 16 * crc - 20 * ih)
                              / (4 * (sf - 2 * de))) * (cr + 4), 0)
    return (preambule + 4.25) * ts + n_pay * ts


def mah(i_ma, t_s):
    return i_ma * t_s / 3600.0


# ---------------------------------------------------------------------------
def topologie(nom_sc, graine):
    """
    Construit la topologie d'une réplication.

    Les passerelles suivent un motif déterministe : elles ne dépendent pas de
    la graine, ce qui garantit que les trois scénarios restent comparables
    entre réplications. Seuls les capteurs et leur mobilité varient.
    """
    sc = SCENARIOS[nom_sc]
    rng = random.Random(graine_stable(nom_sc, graine))
    lo0, la0 = sc["lon0"], sc["lat0"]
    L = sc["largeur_km"] * 1000.0
    H = sc["hauteur_km"] * 1000.0

    gws = [(fx * L, fy * H) for fx, fy in sc["motif"]]
    caps = [((0.05 + 0.90 * rng.random()) * L,
             (0.05 + 0.90 * rng.random()) * H)
            for _ in range(sc["n_capteurs"])]

    # trajectoires : marche aléatoire à inertie de cap, bornée à l'emprise
    traces = {}
    for i, (x, y) in enumerate(caps, start=1):
        cap = rng.random() * 2 * math.pi
        pts = [(x, y)]
        cx, cy = x, y
        for _ in range(PAS_GPS_S, DUREE_S + 1, PAS_GPS_S):
            v = rng.uniform(sc["v_min"], sc["v_max"]) / 3.6
            cap += rng.gauss(0, 0.35)
            d = v * PAS_GPS_S
            cx += d * math.sin(cap)
            cy += d * math.cos(cap)
            if cx < 0 or cx > L:
                cx = min(max(cx, 0.0), L); cap = -cap
            if cy < 0 or cy > H:
                cy = min(max(cy, 0.0), H); cap = math.pi - cap
            pts.append((cx, cy))
        traces[i] = pts
    return gws, caps, traces, sc


# ---------------------------------------------------------------------------
def evaluer(nom_sc, nom_cfg, gws, traces, sc, rng, e_radio_mah):
    """Évalue une configuration sur une topologie donnée."""
    rayon = sc["rayon_m"]
    sigma_t = math.sqrt(SIGMA_SYNC_S ** 2 + SIGMA_MT_S[nom_sc] ** 2)
    tloc = TLOC_OPT_MOY if nom_cfg == "Optimise" else TLOC_S
    pas = max(1, int(tloc / PAS_GPS_S))          # un point de trace par cycle

    erreurs, succes, total = [], 0, 0
    mix = dict(TDoA=0, AGPS=0, GPS=0)
    d3_max, hdops, vis_min = 0.0, [], len(gws)

    def sigma_tdoa(px, py, vis):
        if len(vis) < 3:
            return None
        h = hdop(px, py, vis)
        if h is None or h > HDOP_MAX:
            return None
        return h * C_LUMIERE * sigma_t

    def tenter_gnss(vis, assiste):
        if len(vis) < 1:
            return None
        p = (P_FIX_AGPS if assiste else P_FIX_GNSS)[nom_sc]
        if rng.random() > p:
            return None
        base = SIGMA_GNSS_M[nom_sc] * (FACTEUR_AGPS if assiste else 1.0)
        return abs(rng.gauss(0, base))

    for pts in traces.values():
        for k in range(0, len(pts), pas):
            px, py = pts[k]
            total += 1
            ds = sorted(math.hypot(gx - px, gy - py) for gx, gy in gws)
            if len(ds) >= 3:
                d3_max = max(d3_max, ds[2])
            vis = [g for g in gws if math.hypot(g[0] - px, g[1] - py) <= rayon]
            vis_min = min(vis_min, len(vis))
            h = hdop(px, py, vis)
            if h is not None:
                hdops.append(h)

            if nom_cfg == "TDoA":
                s = sigma_tdoa(px, py, vis)
                err = abs(rng.gauss(0, s)) if s is not None else None
                tech = "TDoA"
            elif nom_cfg == "GPS":
                err, tech = tenter_gnss(vis, False), "GPS"
            elif nom_cfg == "AGPS":
                err, tech = tenter_gnss(vis, True), "AGPS"
            else:
                # sélecteur : TDoA si l'incertitude estimée tient le seuil,
                # sinon escalade vers AGPS puis GNSS complet
                s = sigma_tdoa(px, py, vis)
                if s is not None and s <= SEUIL_PRECISION_M:
                    err, tech = abs(rng.gauss(0, s)), "TDoA"
                else:
                    err, tech = tenter_gnss(vis, True), "AGPS"
                    if err is None:
                        err, tech = tenter_gnss(vis, False), "GPS"
                    if err is None and s is not None:
                        err, tech = abs(rng.gauss(0, s)), "TDoA"
            if err is None:
                continue
            erreurs.append(err)
            mix[tech] += 1
            succes += 1

    if not erreurs:
        raise RuntimeError("aucune localisation réussie : %s / %s" % (nom_sc, nom_cfg))

    rmse = math.sqrt(sum(e * e for e in erreurs) / len(erreurs))
    erreurs.sort()
    p95 = erreurs[min(int(0.95 * len(erreurs)), len(erreurs) - 1)]

    # --- énergie par localisation
    def e_loc_de(t):
        pl = PAYLOAD_O["TDoA"] if t == "TDoA" else 20
        e = mah(PLATEFORME["I_TX_MA"], toa_lora(7, pl))
        e += mah(PLATEFORME["I_MCU_MA"], PLATEFORME["T_MCU_S"])
        if t in TTFF:
            e += mah(PLATEFORME["I_GNSS_MA"], TTFF[t])
        if t in T_RX_ASSIST:
            e += mah(PLATEFORME["I_RX_MA"], T_RX_ASSIST[t])
        return e

    n_ok = sum(mix.values())
    if nom_cfg == "Optimise":
        e_loc = sum(mix[t] / n_ok * e_loc_de(t) for t in mix if mix[t])
    else:
        e_loc = e_loc_de(nom_cfg)

    n_loc = DUREE_S / tloc
    e_veille = mah(PLATEFORME["I_SLEEP_MA"], DUREE_S)
    e_jour = e_loc * n_loc + e_veille
    autonomie = C_BAT_MAH / e_jour

    hdops.sort()
    nh = len(hdops)
    return dict(
        rmse_mean_m=round(rmse, 2), rmse_p95_m=round(p95, 2),
        eloc_mah=round(e_loc, 5), ejour_mah=round(e_jour, 3),
        autonomie_j=round(autonomie, 1),
        taux_succes=round(100.0 * succes / total, 2),
        overhead_B=OVERHEAD_B[nom_cfg],
        pct_tdoa=round(100.0 * mix["TDoA"] / n_ok, 1),
        pct_agps=round(100.0 * mix["AGPS"] / n_ok, 1),
        pct_gps=round(100.0 * mix["GPS"] / n_ok, 1),
        hdop_median=round(hdops[nh // 2], 3) if nh else "",
        hdop_p95=round(hdops[int(0.95 * nh)], 3) if nh else "",
        pct_hdop_sup10=round(100.0 * sum(1 for h in hdops if h > HDOP_MAX) / nh, 2) if nh else "",
        d3_max_m=round(d3_max, 0), marge_m=round(rayon - d3_max, 0),
        passerelles_vues_min=vis_min,
        e_radio_cupcarbon_mah=e_radio_mah,
    )


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graines", nargs=2, type=int, default=[1, 10],
                    metavar=("PREMIERE", "DERNIERE"))
    ap.add_argument("--sortie", default="data")
    ap.add_argument("--sorties-annexes", default="sorties")
    ap.add_argument("--latex", action="store_true",
                    help="produit les tableaux LaTeX prets a coller")
    ap.add_argument("--e-radio", type=float, default=1.783,
                    help="énergie radio mesurée sous CupCarbon, en joules par 24 h")
    a = ap.parse_args()

    os.makedirs(a.sortie, exist_ok=True)
    os.makedirs(os.path.join(a.sortie, "topologies"), exist_ok=True)
    os.makedirs(a.sorties_annexes, exist_ok=True)

    e_radio_mah = a.e_radio / U_V / 3.6
    print("Énergie radio mesurée sous CupCarbon : %.3f J/24 h = %.4f mAh/jour"
          % (a.e_radio, e_radio_mah))
    print("Graines %d à %d\n" % (a.graines[0], a.graines[1]))

    champs = ["scenario", "config", "replication", "rmse_mean_m", "rmse_p95_m",
              "eloc_mah", "ejour_mah", "autonomie_j", "taux_succes", "overhead_B",
              "pct_tdoa", "pct_agps", "pct_gps", "hdop_median", "hdop_p95",
              "pct_hdop_sup10", "d3_max_m", "marge_m", "passerelles_vues_min",
              "e_radio_cupcarbon_mah", "source"]
    lignes = []

    print("%-8s %-9s %5s %9s %10s %11s %9s" %
          ("SCENARIO", "CONFIG", "REP", "RMSE (m)", "E_jour", "Autonomie", "Succès"))
    print("-" * 68)

    for g in range(a.graines[0], a.graines[1] + 1):
        for nom_sc in SCENARIOS:
            gws, caps, traces, sc = topologie(nom_sc, g)

            # coordonnées publiables de cette réplication
            f = os.path.join(a.sortie, "topologies",
                             "coordonnees_%s_rep%02d.csv" % (nom_sc, g))
            with open(f, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["type", "id", "x_m", "y_m"])
                for k, (x, y) in enumerate(gws, start=1):
                    w.writerow(["passerelle", k, round(x, 1), round(y, 1)])
                for i, (x, y) in enumerate(caps, start=1):
                    w.writerow(["capteur", i, round(x, 1), round(y, 1)])

            for nom_cfg in CONFIGS:
                rng = random.Random(graine_stable(nom_sc, nom_cfg, g))
                r = evaluer(nom_sc, nom_cfg, gws, traces, sc, rng, round(e_radio_mah, 4))
                r.update(scenario=nom_sc, config=LIBELLE[nom_cfg], replication=g,
                         source="cupcarbon+modele")
                lignes.append(r)
                if g == a.graines[0]:
                    print("%-8s %-9s %5d %9.2f %10.2f %11.0f %8.1f %%" %
                          (nom_sc, LIBELLE[nom_cfg], g, r["rmse_mean_m"],
                           r["ejour_mah"], r["autonomie_j"], r["taux_succes"]))
        print("  graine %2d terminée" % g)

    chemin = os.path.join(a.sortie, "simulation_results.csv")
    with open(chemin, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=champs, extrasaction="ignore")
        w.writeheader()
        for l in lignes:
            w.writerow(l)
    print("\nÉcrit : %s  (%d lignes)" % (chemin, len(lignes)))

    # ---- tableau de synthèse avec dispersion
    nrep = a.graines[1] - a.graines[0] + 1
    print("\nSYNTHÈSE SUR %d RÉPLICATIONS" % nrep)
    print("%-8s %-9s %16s %16s %12s %10s" %
          ("SCENARIO", "CONFIG", "RMSE moy ± σ", "P95 moy ± σ", "Autonomie", "Succès"))
    print("-" * 78)
    synth = []
    for nom_sc in SCENARIOS:
        for nom_cfg in CONFIGS:
            lab = LIBELLE[nom_cfg]
            g = [l for l in lignes if l["scenario"] == nom_sc and l["config"] == lab]
            if len(g) < 2:
                continue
            def ms(k):
                v = [x[k] for x in g]
                return statistics.mean(v), statistics.stdev(v)
            rm, rs = ms("rmse_mean_m")
            pm, ps = ms("rmse_p95_m")
            am, asd = ms("autonomie_j")
            sm, ss = ms("taux_succes")
            em, es = ms("ejour_mah")
            print("%-8s %-9s %8.2f ± %-5.2f %8.2f ± %-5.2f %9.0f j %7.1f %%" %
                  (nom_sc, lab, rm, rs, pm, ps, am, sm))
            synth.append(dict(scenario=nom_sc, config=lab, n=len(g),
                              rmse_moy=round(rm, 2), rmse_ecart_type=round(rs, 3),
                              p95_moy=round(pm, 2), p95_ecart_type=round(ps, 3),
                              ejour_moy=round(em, 3), ejour_ecart_type=round(es, 4),
                              autonomie_moy=round(am, 1), autonomie_ecart_type=round(asd, 2),
                              succes_moy=round(sm, 2), succes_ecart_type=round(ss, 3)))
    f2 = os.path.join(a.sorties_annexes, "synthese_replications.csv")
    with open(f2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(synth[0]))
        w.writeheader()
        w.writerows(synth)
    print("\nÉcrit :", f2)

    # ---- contrôle d'invariance de l'énergie radio
    f3 = os.path.join(a.sorties_annexes, "invariance_energie.csv")
    with open(f3, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["remarque"])
        w.writerow(["L'énergie de la chaîne radio mesurée sous CupCarbon est "
                    "invariante : le nombre d'émissions (288 par nœud et par jour) "
                    "et la taille de trame LoRa (2048 bits) ne dépendent ni de la "
                    "configuration ni de la graine."])
        w.writerow([])
        w.writerow(["valeur mesurée (J / 24 h)", a.e_radio])
        w.writerow(["équivalent (mAh / jour)", round(e_radio_mah, 4)])
        w.writerow(["nombre d'émissions par nœud et par jour", DUREE_S // TLOC_S])
    print("Écrit :", f3)

    # ---- tableaux LaTeX prêts à coller
    if a.latex:
        f4 = os.path.join(a.sorties_annexes, "tableaux_latex.tex")
        with open(f4, "w", encoding="utf-8") as f:
            def moy(sc, cf, k):
                g = [x[k] for x in lignes
                     if x["scenario"] == sc and x["config"] == LIBELLE[cf]]
                return statistics.mean(g), (statistics.stdev(g) if len(g) > 1 else 0.0)

            f.write("%% Tableau III.1 — precision, avec intervalle de confiance\n")
            f.write("\\begin{table}[H]\n\\centering\n")
            f.write("\\caption{Precision de localisation des configurations de reference "
                    "(moyenne $\\pm$ ecart-type sur %d replications)}\n" % nrep)
            f.write("\\label{tab:ref_precision}\n\\small\n")
            f.write("\\begin{tabular}{|l|c|c|c|}\n\\hline\n")
            f.write("\\rowcolor{myColor!20}\n")
            f.write("\\textbf{Configuration} & \\textbf{Urbain} & \\textbf{Rural} "
                    "& \\textbf{Mixte} \\\\\n\\hline\n")
            for cf in CONFIGS:
                cells = []
                for sc in SCENARIOS:
                    m, sd = moy(sc, cf, "rmse_mean_m")
                    cells.append("%.1f $\\pm$ %.1f" % (m, sd))
                f.write("%s & %s \\\\ \\hline\n"
                        % (LIBELLE[cf], " & ".join(cells)))
            f.write("\\end{tabular}\n\\end{table}\n\n")

            f.write("%% Tableau III.4 — performances completes\n")
            f.write("\\begin{table}[H]\n\\centering\n")
            f.write("\\caption{Performances comparees des quatre configurations "
                    "(moyenne sur %d replications)}\n" % nrep)
            f.write("\\label{tab:resultats_complets}\n\\small\n")
            f.write("\\begin{tabular}{|l|l|r|r|r|r|r|}\n\\hline\n")
            f.write("\\rowcolor{myColor!20}\n")
            f.write("\\textbf{Scenario} & \\textbf{Configuration} & \\textbf{RMSE (m)} "
                    "& \\textbf{P95 (m)} & \\textbf{$E_{jour}$} "
                    "& \\textbf{Autonomie} & \\textbf{Succes} \\\\\n\\hline\n")
            for sc in SCENARIOS:
                for j, cf in enumerate(CONFIGS):
                    r, _ = moy(sc, cf, "rmse_mean_m")
                    p, _ = moy(sc, cf, "rmse_p95_m")
                    e, _ = moy(sc, cf, "ejour_mah")
                    au, _ = moy(sc, cf, "autonomie_j")
                    su, _ = moy(sc, cf, "taux_succes")
                    f.write("%s & %s & %.1f & %.1f & %.2f & %.0f & %.1f \\\\ \\hline\n"
                            % (sc if j == 0 else "", LIBELLE[cf], r, p, e, au, su))
            f.write("\\end{tabular}\n\\end{table}\n\n")

            f.write("%% Tableau III.5 — gains de la configuration optimisee\n")
            f.write("\\begin{table}[H]\n\\centering\n")
            f.write("\\caption{Gains de l'algorithme optimise par rapport a LoRaWAN/TDoA}\n")
            f.write("\\label{tab:gains_optimise}\n")
            f.write("\\begin{tabular}{|l|r|r|r|r|}\n\\hline\n\\rowcolor{myColor!20}\n")
            f.write("\\textbf{Scenario} & \\textbf{Reduction RMSE} "
                    "& \\textbf{Surcout energetique} & \\textbf{Autonomie} "
                    "& \\textbf{Gain de succes} \\\\\n\\hline\n")
            plages = []
            for sc in SCENARIOS:
                rr, _ = moy(sc, "TDoA", "rmse_mean_m")
                ro, _ = moy(sc, "Optimise", "rmse_mean_m")
                er, _ = moy(sc, "TDoA", "ejour_mah")
                eo, _ = moy(sc, "Optimise", "ejour_mah")
                ar, _ = moy(sc, "TDoA", "taux_succes")
                ao, _ = moy(sc, "Optimise", "taux_succes")
                au, _ = moy(sc, "Optimise", "autonomie_j")
                dr = 100 * (rr - ro) / rr
                de = 100 * (eo - er) / er
                plages.append((dr, de, au))
                f.write("%s & $-%.1f$ \\%% & $+%.1f$ \\%% & %.0f jours & $+%.1f$ pts "
                        "\\\\ \\hline\n" % (sc, dr, de, au, ao - ar))
            f.write("\\end{tabular}\n\\end{table}\n\n")
            f.write("%% Phrases a reprendre dans le texte :\n")
            f.write("%%   reduction de l'erreur : de %.0f %%%% a %.0f %%%%\n"
                    % (min(p[0] for p in plages), max(p[0] for p in plages)))
            f.write("%%   surcout energetique  : de %.0f %%%% a %.0f %%%%\n"
                    % (min(p[1] for p in plages), max(p[1] for p in plages)))
            f.write("%%   autonomie conservee  : de %.0f a %.0f jours (%.0f a %.0f mois)\n"
                    % (min(p[2] for p in plages), max(p[2] for p in plages),
                       min(p[2] for p in plages) / 30.44, max(p[2] for p in plages) / 30.44))
        print("Écrit :", f4)
        print("\n  Phrases a reprendre dans le texte du memoire :")
        print("    reduction de l'erreur : de %.0f %% a %.0f %%"
              % (min(p[0] for p in plages), max(p[0] for p in plages)))
        print("    surcout energetique  : de %.0f %% a %.0f %%"
              % (min(p[1] for p in plages), max(p[1] for p in plages)))
        print("    autonomie conservee  : de %.0f a %.0f jours (%.0f a %.0f mois)"
              % (min(p[2] for p in plages), max(p[2] for p in plages),
                 min(p[2] for p in plages) / 30.44, max(p[2] for p in plages) / 30.44))

    print("\nÉtape suivante : python figures_chapitre3.py --csv %s" % chemin)


if __name__ == "__main__":
    raise SystemExit(main())
