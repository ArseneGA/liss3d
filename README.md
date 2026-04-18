# Lissajous 3D Tube Generator

Voisinage tubulaire R d'une courbe de Lissajous 3D, extrait en mesh solide manifold via SDF + Dual Contouring. Imprimable STL.

## Quick start

```bash
npm install
npm run dev
```

Ouvrir http://localhost:5173.

## Étape 1 — test manuel

- Les sliders `p`, `q`, `r` (fréquences, réels dans [0, 10]) modifient la courbe visible en temps réel.
- `φx, φy, φz` font tourner/déphaser chaque axe.
- La courbe apparaît en bleu ; axes RGB, grille au sol.

## Étape 2 — test manuel

- Le panneau **Validation** affiche `κ_max`, `1/κmax`, `N pts`, `h`.
- Quand `R` approche ou dépasse `0.8 / κ_max`, un warning/erreur apparaît.
- `N` augmente automatiquement avec `p, q, r` (max |γ'| = √(p²+q²+r²)) et avec `κ_max`.
- En dev, ouvrir la console : self-tests dans `[selftest] lissajous` — écarts γ'/γ'' vs différence finie (< 1e-4 et < 1e-2), et `estimateMaxCurvature` à deux résolutions (écart relatif faible).

## Étape 3 — test manuel

Dans la console de dev (groupes `[selftest] sdf …`) :

- **single segment** : pour chaque voxel dans la bande `d ≤ R+3h`, compare `field[idx]` à `dist(centre_voxel, segment) - R` calculé analytiquement. `maxErrInBand` doit être ≈ 0 (< 1e-6 en pratique, flottants float32).
- **sphere** : idem, pour un point (segment dégénéré). `maxErrInBand` ≈ 0.
- **smooth mode** :
  - `singleSegment_AvsB_maxDiff` ≈ 0 (avec un seul segment, smin = min).
  - `unionDeltaFromMin ≤ 0` (l'union rasterisée n'est pas plus grande que chaque SDF individuelle).
  - `smoothVsSharpAtCenter < 0` (mode B creuse au croisement, fillet vers l'intérieur).

## Étape 4 — test manuel

Bouton **Compute mesh** (vert) pour générer le maillage du tube à partir des paramètres courants. Le panneau affiche `V / F`, `χ (if closed)`, et les temps (sample, sdf, dc, total).

Dans la console :

- **dc sphere** : SDF analytique `|p| − R` sur grille 48³. `eulerIfClosed` doit valoir **2.0**. `maxRadiusErr` ≲ `h` (0.04).
- **dc cube** : SDF analytique d'un cube `[-0.3, 0.3]³`. `eulerIfClosed` doit valoir 2.0 et les arêtes doivent visuellement apparaître vives (QEF → coins nets).

## Étape 5, 6, 7 — test manuel

- **Auto rebuild** : bouger les sliders déclenche un rebuild 200 ms après la dernière modification ; le mesh change tout seul.
- **Bouton Compute mesh** : force un rebuild immédiat (bypass debounce).
- **Panneau Validation** : affiche `V, E, F`, `χ`, flags `watertight / manifold / orientable ✓/✗`, et `degenerate`.
- **Bouton Export STL** : activé seulement si `watertight && manifold`. Génère un `.stl` binaire recentré et scalé pour que la dimension max vaille `scaleMm` mm. Nom de fichier : `lissajous_a-b-c_Rxx_res_modeX.stl`.
- **UI non-bloquante** : tout le pipeline tourne dans un Web Worker. Les sliders restent réactifs pendant le calcul.
- **Self-tests console** :
  - `validate sphere` : tous flags à `true`, `χ = 2`.
  - `validate broken` : tétraèdre amputé → `watertight = false`.
