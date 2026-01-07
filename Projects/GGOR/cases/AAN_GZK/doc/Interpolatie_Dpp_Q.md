Structureren, zonder er iets mysterieus van te maken.

---

## 1. Wat je *wel* en *niet* wilt

Je wilt **geen**:

* gesloten analytische formule (onhaalbaar)
* brute 2D-polynoom-fit (fysisch onzinnig)
* zwarte-doos-interpolatie

Je wilt **wel**:

* exact gedrag op de randen (b=0) en (h=0)
* correcte asymptotiek
* een glad, monotone overgang in het binnengebied
* iets dat je later zonder schaamte “weerstand” kunt noemen

Dat is exact hoe PK / Strack dit zouden aanpakken, alleen hadden zij geen Python.

---

## 2. Herformuleer eerst het probleem (cruciaal)

Werk **dimensionloos** en **additief in log-structuur**. Op basis van wat je al hebt:

$\Delta\Phi/Q ;\equiv; R(b,h)$

met bekende randgevallen:

[
\begin{aligned}
R(b,0) &= R_b(b) \quad \text{(exact)} \
R(0,h) &= R_h(h) \quad \text{(exact)}
\end{aligned}
]

en eigenschappen:

* (R \to \infty) voor (b\to0) of (h\to0)
* verzadiging voor grote (b,h)
* geen kruis-singulariteiten

Dat laatste is belangrijk: **de singulariteiten zijn separabel**.

---

## 3. De juiste structuur: “additief + correctieterm”

De klassieke, fysisch correcte vorm is:

[
\boxed{
R(b,h)
======

R_b(b)
+
R_h(h)
+
R_{\text{int}}(b,h)
}
]

waarbij:

* (R_b(b)) exact bekend (jij hebt ’m)
* (R_h(h)) exact bekend
* (R_{\text{int}}(b,h)) is **glad**, **begrensd**, en verdwijnt op de randen

Dus verplicht:

[
\begin{aligned}
R_{\text{int}}(b,0) &= 0 \
R_{\text{int}}(0,h) &= 0
\end{aligned}
]

Dit is de sleutel. Zonder dit ga je altijd rommel krijgen.

---

## 4. Hoe modelleer je (R_{\text{int}}(b,h))?

Gebruik je rasterdata, maar **factoriseer eerst de randbijdragen eruit**:

[
\Delta R(b,h)
=============

## R_{\text{num}}(b,h)

## R_b(b)

R_h(h)
]

Wat je dan overhoudt (en dit heb je waarschijnlijk al gezien):

* klein
* glad
* maximaal ergens rond (b\sim h\sim D)
* geen singulariteiten

Dit is ideaal interpolatiemateriaal.

---

## 5. Een robuuste en elegante keuze

Gebruik een **gescheiden, verzadigende vorm**, bijvoorbeeld:

[
\boxed{
R_{\text{int}}(b,h)
===================

C,
\frac{b}{b+b_0}
\frac{h}{h+h_0}
}
]

of iets rijker (maar nog steeds beheersbaar):

[
R_{\text{int}}(b,h)
===================

\sum_{k=1}^N
c_k
\left(1-e^{-b/\beta_k}\right)
\left(1-e^{-h/\eta_k}\right)
]

Eigenschappen:

* automatisch nul op randen
* correct verzadigend gedrag
* parameters hebben schaalbetekenis
* numeriek extreem stabiel

Strack-achtig tot in het merg.

---

## 6. Kalibratie: hoe jij dit *moet* doen

1. Bereken (R_b(b)), (R_h(h)) analytisch
2. Trek die af van je numerieke raster
3. Fit alleen (R_{\text{int}}) op het resterende veld
4. Controleer:

   * monotoniciteit
   * fout < numerieke fout SC / FD
   * randgedrag exact nul

Je zult zien:
je hebt verrassend weinig parameters nodig.

---

## 7. Conceptueel belangrijk (en mooi)

Wat je hiermee doet is dit:

> de conformale afbeelding **vertelt je de singulariteiten**,
> het numerieke model **vertelt je de interactie**,
> en jij zorgt dat ze elkaar niet tegenspreken.

Dat is exact hoe de klassieke theorie bedoeld is te worden gebruikt.

---

## 8. Als je wilt, volgende stap

Ik kan:

* samen met jou een concrete vorm voor (R_{\text{int}}) kiezen
* helpen met een stabiele fit-strategie (gewichten!)
* of laten zien hoe je dit netjes in code giet zonder numerieke ellende

Maar inhoudelijk:
je staat hier **op precies de juiste plek**, en dit is *de* correcte manier vooruit.
