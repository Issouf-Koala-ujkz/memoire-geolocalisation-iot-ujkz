# Optimisation des technologies de géolocalisation sous consommation énergétique limitée

Code, données et projets de simulation du mémoire de Master de **KOALA Issouf**.

Ce dépôt permet de rejouer l'intégralité des résultats analytiques du mémoire à
partir des seules sources publiées, sans données externes.

---

## 1. Contexte académique

| | |
|---|---|
| **Auteur** | KOALA Issouf |
| **Mémoire** | *Optimisation des technologies de géolocalisation sous consommation énergétique limitée* |
| **Diplôme** | Master 2 en Informatique, option Réseaux Émergents, IoT et Intelligence Embarquée (RIoT/IE) |
| **Établissement** | UFR-SEA, Université Joseph KI-ZERBO (UJKZ), Ouagadougou, Burkina Faso |
| **Laboratoire** | LAMI |
| **Année académique** | 2025–2026 |

---

## 2. Démarrage rapide

```bash
git clone https://github.com/Issouf-Koala-ujkz/memoire-geolocalisation-iot-ujkz.git
cd memoire-geolocalisation-iot-ujkz

pip install numpy scipy matplotlib

python replications.py
python generations_figures.py --csv data/simulation_results.csv --simulations ./Simulations
```

La première commande régénère les 120 lignes de résultats et les 30 topologies.
La seconde produit les tests statistiques, les intervalles de confiance à 95 %
et les cartes de HDOP.

---

## 3. Arborescence

```
.
├── replications.py            Construit les topologies et évalue les configurations
├── generations_figures.py     Tests statistiques, IC 95 %, cartes de HDOP
│
├── data/
│   ├── simulation_results.csv 120 lignes, 21 colonnes — jeu de données principal
│   └── topologies/            30 fichiers de coordonnées (3 scénarios × 10 réplications)
│
├── Simulations/               12 projets CupCarbon
│   ├── Scenario_Urbain/       Urbain_GPS, Urbain_AGPS, Urbain_TDoA, Urbain_Optimise
│   ├── Scenario_Rural/        Rural_GPS, Rural_AGPS, Rural_TDoA, Rural_Optimise
│   └── Scenario_Mixte/        Mixte_GPS, Mixte_AGPS, Mixte_TDoA, Mixte_Optimise
│
├── sorties/
│   ├── synthese_replications.csv  Moyennes et écarts-types par cellule
│   ├── invariance_energie.csv     Contrôle de l'invariance du terme radio
│   └── tableaux_latex.tex         Tableaux insérés dans le mémoire
│
├── figures/                   Cartes de HDOP générées par les scripts
└── ps/                        Figures et illustrations du mémoire
```

---

## 4. Plan d'expérience

**3 scénarios** × **4 configurations** × **10 réplications** = **120 observations**

| | Valeurs |
|---|---|
| Scénarios | Urbain, Rural, Mixte |
| Configurations | GPS, A-GPS, TDoA, Optimisé |
| Réplications | graines 1 à 10 |

Chaque ligne de `data/simulation_results.csv` porte 21 colonnes : RMSE moyen et
au 95ᵉ centile, énergie par localisation et par jour, autonomie, taux de
succès, surcoût protocolaire, répartition des technologies, HDOP médian et au
95ᵉ centile, distance maximale, marge, nombre de passerelles vues, et énergie
radio mesurée sous CupCarbon.

---

## 5. Reproductibilité

### 5.1 Ce qui est exactement reproductible

Les tirages aléatoires sont dérivés d'une graine déterministe construite à
partir du nom du scénario, du nom de la configuration et du numéro de
réplication (fonction `graine_stable`). **Deux exécutions successives
produisent des fichiers strictement identiques**, sur n'importe quelle machine.

Sont recalculés intégralement par script :

- positions des nœuds et des passerelles, traces de mobilité, portée radio
- HDOP et erreur de localisation TDoA
- énergies GNSS, microcontrôleur et veille
- autonomie et taux de succès
- tests statistiques et intervalles de confiance

### 5.2 Ce qui ne l'est pas

**CupCarbon ne dispose pas de mode ligne de commande** : sa classe principale
n'accepte en argument qu'une configuration de serveur mandataire. Les douze
exécutions du simulateur ont donc été conduites manuellement depuis l'interface
graphique et ne sont pas scriptables.

Une seule grandeur en provient : **l'énergie de la chaîne radio**. Le fichier
`sorties/invariance_energie.csv` atteste qu'elle reste invariante d'une
configuration à l'autre, ce qui justifie de la traiter comme une constante dans
le modèle énergétique.

Les douze projets sont publiés avec l'intégralité de leur configuration, de
sorte qu'un lecteur puisse les rouvrir et les relancer à l'identique. La valeur
du terme radio peut par ailleurs être passée en argument, ce qui permet de
tester la sensibilité des conclusions sans relancer le simulateur :

```bash
python replications.py --e-radio 1.783
```

---

## 6. Prérequis

| Élément | Version | Statut |
|---|---|---|
| Python | 3.8 ou plus | obligatoire |
| numpy | — | obligatoire |
| scipy | — | recommandé (tests statistiques) |
| matplotlib | — | recommandé (cartes de HDOP) |
| CupCarbon | — | facultatif, pour rouvrir les projets de simulation |

En l'absence de `scipy` ou de `matplotlib`, les scripts s'exécutent et
signalent les sorties qu'ils n'ont pas pu produire.

---

## 7. Limites

- **Modèle énergétique fondé sur des fiches techniques**, non sur des mesures
  matérielles. Une analyse de sensibilité aux constantes de calibration est
  fournie par `generations_figures.py`.
- **Aucun déploiement réel** : l'ensemble des résultats est issu de simulation.
- **Les sources LaTeX du mémoire ne sont pas publiées.** Les tableaux de
  résultats insérés dans le document sont en revanche générés automatiquement
  dans `sorties/tableaux_latex.tex`, ce qui permet de confronter les valeurs
  citées dans le texte à celles produites par les scripts.

---

## 8. Contact

**KOALA Issouf** — [@Issouf-Koala-ujkz](https://github.com/Issouf-Koala-ujkz)

Signalements et demandes de précision via les *issues* du dépôt.
