# Rapport de diagnostic — Crashs NR du modèle `NGammaTiTeNeutralEuler`

**Date :** 12 juin 2026
**Données analysées :** `/home/obrun/Euler_data/Solutions/*.h5` (runs TCABR splined 3K, P8, DPe = 0.933 et DPe = 20)
**Code audité :** `src/Models/NGammaTiTeNeutralEuler/physics.f90`, `src/HDG/hdg_ComputeJacobian.f90` (diff vs HEAD), `src/HDG/hdg_BC.f90` (diff vs HEAD), `test/param.txt`

---

## 1. Résumé exécutif

Le crash n'est **pas** dû au réglage des diffusions perpendiculaires (les deux runs DPe = 0.933 et DPe = 20 divergent de manière identique), mais à une combinaison de bugs structurels dans l'implémentation du sous-système Euler des neutres. Le bug principal est une **inversion de signe des termes sources de quantité de mouvement des neutres** dans `assemblyNeutral` (`hdg_ComputeJacobian.f90`) : le terme de friction ionisation/CX, qui devrait amortir Γₙ, agit comme une **anti-friction** (croissance exponentielle), et son signe dans la jacobienne rend la matrice locale quasi singulière dès la première itération NR (le taux de perte n·⟨σv⟩ dépasse 1/Δt avec Δt = 100). Les itérés NR produisent alors des Γₙ ~ O(1) parasites, puis le **clamp ρₙ → max(ρₙ, 10⁻⁷)** transforme ces momenta en vitesses uₙ ~ 10⁵–10⁷, qui contaminent la jacobienne de convection (~uₙ²), la stabilisation τ (~|uₙ|) et la source d'énergie ionique fEiN (~|Γₙ|²/ρₙ), provoquant l'explosion complète à NR2.

---

## 2. Constat factuel d'après les fichiers `.h5`

Évolution des champs (run DPe = 0.933, identique en structure pour DPe = 20) :

| Itération | État de la solution |
|---|---|
| `0000` (init) | ρ ∈ [0.044, 1], Γ = 0, ρₙ = 10⁻⁸ uniforme, Γₙₓ = Γₙᵧ = 0, q ≈ 0. État sain. |
| `NR0001` | ρₙ ∈ [−7×10⁻⁴, 1.6×10⁻²] avec **27 844 points négatifs (18 %)** ; Γₙₓ jusqu'à 1.36, Γₙᵧ jusqu'à 3.12 (en unités n₀·c_s0, soit un flux de neutres comparable au flux sonique plasma !) ; vitesse de neutres reconstruite max ≈ 3×10⁵ (adim.). Plasma encore quasi intact. |
| `NR0002` | Explosion totale : ρ ∈ [−5×10⁸, 2×10⁹], nEi jusqu'à −3×10¹⁰, Γₙ ~ 10¹⁰, Ti < 0 sur 15 235 points, max\|q\| ~ 10¹¹. |

**Localisation :** la divergence de NR1 est concentrée autour de (R, Z) ≈ (0.61, −0.16) — région X-point/jambe de divertor, à ~6 cm du bord (donc un mécanisme **volumique**, pas une CL) — c'est-à-dire là où n·(σ_iz + σ_cx) et le couplage source plasma↔neutres sont maximaux. Le puff (face en (0.525, −0.245)) remplit ρₙ jusqu'à ~10⁻² sans poser problème en soi.

Paramètres du run : `stab = 5`, `tau = 1`, `dt0 = 100` (≈ 1.4×10⁻⁵ s), `nrp = 3` (3 itérations NR max), `shockcp = 0`, `limrho = 0`, `thresh = 0` (aucun limiteur/capture de choc actif), Mref = 12.49.

---

## 3. Bugs identifiés

### Bug n°1 (critique) — Signes inversés des sources de quantité de mouvement des neutres

Fichier : `src/HDG/hdg_ComputeJacobian.f90`, `assemblyNeutral`, blocs `#ifdef NEUTRALEULER`.

La convention du code (vérifiée sur les équations 1–4 et sur le modèle `NEUTRALGAMMA`) est :
- `Sn(i,:)` = **−∂S/∂U** (la jacobienne de la source, signe moins, ajoutée au membre de gauche via `Auu += Sn·NNi`) ;
- `Sn0(i)` = **−S(U₀) + S′(U₀)·U₀** (ajouté au RHS via `rhs −= Sn0`).

Pour `NEUTRALGAMMA`, la cohérence est garantie par construction : `Sn(ign,:) = −Sn(2,:)` et `Sn0(ign) = −Sn0(2)`.

Pour `NEUTRALEULER`, les termes ont été réécrits à la main **avec le signe de +∂S/∂U au lieu de −∂S/∂U** :

```4205:4226:src/HDG/hdg_ComputeJacobian.f90
#ifdef NEUTRALEULER
      Sn(ignx, :) = (dfGammacx_dU(:)*sigmavcx + dfGammarec_dU(:)*sigmavrec)*b(1)
      Sn(igny, :) = (dfGammacx_dU(:)*sigmavcx + dfGammarec_dU(:)*sigmavrec)*b(2)
      ! ...
        Sn(ignx, :) = Sn(ignx, :) - dloss_rate_dU(:)*U(ignx)
        Sn(ignx, ignx) = Sn(ignx, ignx) - loss_rate
```

La source physique attendue pour Γₙₓ est S = +b₁·(fΓcx·σcx + fΓrec·σrec) − Γₙₓ·n·(σiz + σcx). Avec la convention du code, on devrait avoir `Sn(ignx,:) = −b₁·d(gain)/dU + d(perte)/dU`, soit exactement **l'opposé** de ce qui est écrit. Conséquences :

1. **Terme de perte (friction) anti-amorti :** sur la diagonale du système, le terme `−loss_rate` (au lieu de `+loss_rate`) s'ajoute à la matrice de masse M/Δt. Estimation au point chaud : σ_adim = σ_phys·n₀·t₀ ≈ 4×10⁻², n ≈ 0.15–0.4 ⇒ loss ≈ 10⁻²–4×10⁻², tandis que 1/Δt = 10⁻². **Le terme déstabilisant est du même ordre ou supérieur à l'inertie** : l'opérateur implicite local pour Γₙ devient indéfini/quasi singulier précisément dans la zone X-point/divertor — c'est ce qui produit les Γₙ ~ O(1) parasites de NR1 au bon endroit.
2. **Résidu incohérent avec la jacobienne :** en combinant `Sn` et `Sn0` tels qu'écrits, la source effective au point de linéarisation vaut ≈ **−3×(gain CX/rec) + 2×(perte)·Γₙ** au lieu de +gain − perte·Γₙ. Même si l'algèbre linéaire convergeait, la dynamique discrète serait celle d'un système anti-amorti (croissance exponentielle de Γₙ au taux ~2·n·(σiz+σcx)).
3. Incohérence secondaire : dans `Sn0(ignx)`, le terme `fGammacx·dsigmavcx_dU·U·b(1)` est absent (seul le terme `fGammarec·dsigmavrec_dU` y figure), alors que son homologue est présent dans `Sn(ignx,:)`. La linéarisation de Newton est donc incomplète même à signes corrigés.

### Bug n°2 (critique) — Reconstruction de vitesse uₙ = Γₙ/max(ρₙ, 10⁻⁷) non bornée

Présent dans les trois copies du bloc `Ax_neut/Ay_neut` (`hdg_ComputeJacobian.f90`, volume + faces int. + faces ext.), dans `cons2phys` (`physics.f90`) et dans `computeTauGaussPoints` :

```2192:2197:src/HDG/hdg_ComputeJacobian.f90
      nn = MAX(uf(5), tol)        ! tol = 1.d-7
      unx = uf(6)/nn
      uny = uf(7)/nn
```

- Le clamp ne protège que contre la division par zéro, **pas contre l'incohérence physique** : dès que ρₙ < 10⁻⁷ (ou négatif — `MAX` sans `ABS` ici), uₙ = Γₙ/10⁻⁷ explose. À NR1 on mesure uₙ ~ 3×10⁵ là où ρₙ > 0, et ~10⁷ là où ρₙ < 0.
- Ces vitesses entrent **au carré** dans `Ax_neut(6,5) = −unx² + Ti` (termes ~10¹⁰–10¹⁴ dans la jacobienne de convection), **linéairement** dans τ (`tau_euler = |un·n| + cs_n` ~ 10⁵–10⁷, stabilisation devenue elle-même pathologique), et au carré dans `fEiN = ½·n·|Γₙ|²/ρₙ` ~ 10⁷ qui est injecté comme source dans l'équation d'énergie ionique — d'où le nEi ~ −3×10¹⁰ observé à NR2.
- **Aggravation :** l'état initial ρₙ = 10⁻⁸ (`analytical.f90`) est **inférieur au clamp 10⁻⁷**. Tout le domaine démarre donc dans le régime saturé du clamp, où la jacobienne (écrite pour la branche non clampée : termes `−unx²`, `2·unx`, etc.) ne correspond pas au flux réellement évalué — l'homogénéité F = A·U, sur laquelle repose l'assemblage HDG de la convection, est brisée dès t = 0.

### Bug n°3 (majeur) — Aucun mécanisme de positivité/régularisation pour le sous-système neutre

`setLocalDiff` (`physics.f90`, l. 1323–1329) met **explicitement à zéro** toute diffusion pour les équations 5–7 :

```1323:1329:src/Models/NGammaTiTeNeutralEuler/physics.f90
#ifdef NEUTRALEULER
    d_iso(phys%idx_rhon_eq, phys%idx_rhon_eq, :) = 0.d0
    d_iso(phys%idx_gammanx_eq, phys%idx_gammanx_eq, :) = 0.d0
    d_iso(phys%idx_gammany_eq, phys%idx_gammany_eq, :) = 0.d0
```

Le sous-système neutre est donc purement hyperbolique, discrétisé en P8, **sans diffusion physique, sans capture de choc (`shockcp = 0`), sans limiteur de positivité (`limrho = 0`, `thresh = 0`)**. Les oscillations de Gibbs au front du puff et dans la zone source rendent ρₙ < 0 inévitable (18 % des points dès NR1) ; combiné au Bug n°2, c'est l'amplificateur de l'instabilité. Pour un système d'Euler raide, l'absence totale de viscosité artificielle ou de limiteur est intenable numériquement.

### Bug n°4 (majeur) — Termes jacobiens orphelins de la diffusion Dₙₙ supprimée

La diffusion de base `Dₙₙ∇ρₙ` a été retirée (d_iso = 0), mais les **termes de correction de Newton** de cette diffusion non linéaire sont toujours assemblés sans garde `#ifndef NEUTRALEULER`, en volume :

```2755:2766:src/HDG/hdg_ComputeJacobian.f90
          ELSEIF (i == inn) THEN
                 DO j = 1,Neq
              z = i+(j-1)*Neq
                    DO k = 1,Ndim
                Auu(:,:,z) =Auu(:,:,z) + (NxyzNi(:,:,k)*Dnn_dU(j)*Qpr(k,i))
```

et aux faces (l. 3444, 3973). Le résidu stationnaire n'est pas modifié (les deux corrections se compensent au point de linéarisation), mais **la jacobienne contient des couplages massifs qui ne correspondent à aucun terme du résidu** : Dₙₙ est calculé via `compute_Dnn` avec saturation à `diff_nn = 7.6×10⁶` (adim.), et ces termes sont multipliés par Qpr(:,inn) = ∇ρₙ qui vaut ~30 dès NR1. Direction de Newton faussée garantie dès que ∇ρₙ ≠ 0.

### Bug n°5 (physique) — Facteur Mref manquant dans la pression des neutres

Dans `Ax_neut/Ay_neut`, la pression neutre est implémentée comme pₙ = ρₙ·Ti (hypothèse Tₙ = Ti). Or dans l'adimensionnement du code, la pression en unités de flux de quantité de mouvement vaut p = Mref·ρ·T (cf. plasma : pi = 2/3·(U₃ − ½U₂²/U₁) = Mref·ρ·Ti). Avec **Mref = 12.49**, la force de pression des neutres et la vitesse du son neutre (c_n = √Ti au lieu de √(Mref·Ti)) sont sous-estimées d'un facteur ~12.5 (~3.5 sur c_n). Ce bug est « stabilisant » (pression trop faible) mais fausse la physique et rend τ_euler sous-évalué par rapport au vrai rayon spectral. Les échanges de quantité de mouvement par CX/recombinaison étant, eux, dans les bonnes unités, le système couplé est dimensionnellement incohérent.

### Bug n°6 (mineur mais latent) — Tableau `numer%tau(1:5)` pour 7 équations

`types.f90` déclare `tau(1:5)` ; le modèle a Neq = 7. Avec `stab = 1`, la boucle `DO i = 1,Neq : tau(i,i) = numer%tau(i)` lit **hors bornes** (`numer%tau(6:7)`). Le run analysé utilise `stab = 5` (non affecté car le bloc NEUTRALEULER de `computeTauGaussPoints` n'utilise que `numer%tau(5)`), mais c'est une bombe à retardement.

### Observations secondaires

- **Terme géométrique axisymétrique manquant :** l'assemblage axisymétrique (dvolu ∝ R) fait que la divergence du flux de quantité de mouvement radial calcule effectivement 1/R·∂_R(R·(ρₙuₙₓ² + pₙ)) ; le terme de compensation (« hoop stress » +pₙ/R, +ρₙu_θ²/R) n'est pas ajouté en source. Erreur de physique en géométrie torique, sans rôle dans le crash.
- **CL de Bohm pour Γₙ :** la trace est imposée à û_Γₙ = δ·R_rec·(−u_i‖·b)·ρₙ via la pénalisation τ. Quand δ = 0 (incidence quasi tangente ou écoulement supersonique), la CL dégénère en û_Γₙ = 0 strict. Aucun bilan de flux (pression + convection) n'est écrit pour les lignes 6–7 au mur — c'est un choix de type Dirichlet acceptable, mais à documenter ; le couplage recyclage utilise ρₙ·u_i (convention héritée de NeutralGamma) et non Γ_i.
- `jacobianMatricesN` (physics.f90) calcule Ax/Ay puis les **écrase par zéro** (lignes 1192–1193) — code mort hérité, sans effet mais trompeur.
- Indices 5/6/7 codés en dur dans les trois blocs `Ax_neut` et les branches `if (i >= 5 .and. i <= 7)` : incompatible avec `KEQUATION` et fragile vis-à-vis de `phys%idx_*`.
- Paramétrage du run aggravant : Δt initial = 100 (≈ 14 µs) pour démarrer un sous-système raide depuis ρₙ = 10⁻⁸, avec seulement `nrp = 3` itérations NR — même sans les bugs ci-dessus, la convergence NR du transitoire de remplissage serait très improbable.

---

## 4. Scénario de défaillance reconstitué

1. **t = 0 :** ρₙ = 10⁻⁸ (sous le clamp 10⁻⁷), Γₙ = 0, Γ = 0, Δt = 100.
2. **NR1 :** la jacobienne des équations Γₙₓ/Γₙᵧ porte le terme de friction au **mauvais signe** (−n(σiz+σcx) sur la diagonale, comparable ou supérieur à 1/Δt) → opérateur local mal conditionné/indéfini dans la zone X-point/divertor → la résolution produit Γₙ ~ O(1) parasites, maximaux en (0.61, −0.16). En parallèle, le puff remplit ρₙ, et l'absence de toute diffusion/limiteur sur un champ P8 crée des oscillations : ρₙ < 0 sur 18 % des points.
3. **NR1 → NR2 :** uₙ = Γₙ/max(ρₙ, 10⁻⁷) atteint 10⁵–10⁷ ; τ_euler ~ |uₙ| et Ax_neut ~ uₙ² explosent ; fEiN ~ n|Γₙ|²/10⁻⁷ ~ 10⁷ injecte une source démesurée dans nEi ; la source anti-amortie +2·loss·Γₙ amplifie encore Γₙ ; les termes orphelins Dₙₙ_dU·∇ρₙ (∇ρₙ ~ 30) polluent la direction de Newton de l'équation ρₙ.
4. **NR2 :** solution à 10⁹–10¹⁰, Ti < 0 sur 73 000 points, ρ < 0 → toutes les fermetures (taux EIRENE, c_s, α_i/α_e) sont évaluées hors domaine → NR3 inutilisable, crash (matrices singulières/overflow).

L'insensibilité au DPe (0.933 vs 20) confirme que le mécanisme est porté par le nouveau sous-système neutre et ses sources, pas par le transport perpendiculaire du plasma.

---

## 5. Recommandations (par ordre de priorité, sans code)

1. **Corriger les signes** des blocs `NEUTRALEULER` dans `assemblyNeutral` (Sn(ignx/igny,:) et Sn0(ignx/igny)) pour respecter la convention Sn = −∂S/∂U, Sn0 = −S + S′·U₀, et compléter le terme `fGammacx·dsigmavcx_dU` manquant dans Sn0. Vérifier par différences finies (jacobienne vs résidu) sur un élément.
2. **Borner la reconstruction de uₙ de manière cohérente** : plancher de ρₙ nettement supérieur (et identique partout : assemblage volume/faces, cons2phys, τ, fEiN), usage d'ABS cohérent, et surtout jacobienne cohérente avec la branche clampée (dérivées gelées dans le régime saturé). Initialiser ρₙ au-dessus du plancher.
3. **Réintroduire une régularisation du sous-système neutre** : viscosité artificielle/diffusion résiduelle sur les équations 5–7, capture de choc activée, et idéalement un traitement de positivité pour ρₙ. Un système d'Euler P8 sans aucune dissipation n'est pas viable.
4. **Supprimer (ou garder cohérents) les termes Dₙₙ_dU** sous `NEUTRALEULER` en volume et aux faces.
5. **Restaurer le facteur Mref** dans pₙ et dans la vitesse du son neutre de τ_euler.
6. Dimensionner `numer%tau` à Neq (ou ≥ 7) pour éliminer l'accès hors bornes en `stab = 1`.
7. Pour la mise en route : Δt initial réduit de plusieurs ordres de grandeur (le temps caractéristique d'ionisation adim. est ~25), `nrp` ≫ 3, éventuellement rampe sur le puff ; ajouter le terme géométrique axisymétrique de l'équation de quantité de mouvement radiale lors d'une passe « physique ».

---

## 6. Annexe — éléments de preuve quantitatifs

- Taux adimensionnels : σ·n₀·t₀ ≈ 4×10⁻² (σ ≈ 3×10⁻¹⁴ m³/s, n₀ = 10¹⁹ m⁻³, t₀ = 1.374×10⁻⁷ s) ; avec n ≈ 0.15–0.4 dans la zone critique, loss ≈ 0.6–4×10⁻² ≥ 1/Δt = 10⁻².
- NR1 : max|Γₙᵧ| = 3.12 en (0.606, −0.163) où ρₙ = 10⁻⁵ ⇒ uₙᵧ ≈ 3×10⁵ ; min ρₙ = −7×10⁻⁴ en (0.611, −0.161).
- NR2 : max|nEi| = 3.06×10¹⁰ en (0.607, −0.156) — colocalisé avec le défaut de NR1, signature de la chaîne fEiN/Ax_neut.
- fEiN à NR1 (zone critique) : ½·0.4·(1.36² + 3.12²)/10⁻⁷ ≈ 2×10⁷ (adim.) injecté dans l'équation nEi.
