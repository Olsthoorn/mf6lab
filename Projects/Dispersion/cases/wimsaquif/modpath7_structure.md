# Overzicht MODPATH 7 modelopbouw met FloPy

Dit document geeft een overzicht van de stappen en bijbehorende FloPy-objecten die nodig zijn om een MODPATH 7 model succesvol op te bouwen en te draaien. Dit is bedoeld als naslagwerk met structuur en inzicht in de samenhang van de verschillende onderdelen.

## 1. Voorwaarden
- Een volledig en succesvol doorgerekend MODFLOW 6 model met:
  - `.dis.grb` bestand (grid bestand)
  - `.hds` bestand (heads)
  - `.cbc` bestand (cell budget)

## 2. Structuur van een MODPATH 7 model

```text
MODPATH 7 projectstructuur:

ProjectFolder/
├── GWF/                  # MODFLOW model en output
│   ├── model.nam         # MODFLOW naamfile
│   ├── model.dis.grb     # Grid-bestand
│   ├── model.hds         # Heads
│   └── model.cbc         # Budget
├── MP7/                  # MODPATH model
│   ├── model.mpsim       # MP7 simulatiebestand
│   ├── model.mpbas       # MP7 basisbestand
│   └── ...               # Overige MP7 bestanden
```

## 3. Benodigde FloPy-objecten

### 3.1 Laad het MODFLOW model
```python
flowmodel = flopy.mf6.MFSimulation.load(...)
```

### 3.2 Maak het MODPATH 7 model
```python
mp7 = flopy.modpath.Modpath7(
    modelname='model',
    flowmodel=flowmodel,
    model_ws='MP7',
    exe_name='mp7'
)
```

### 3.3 Particles definiëren
```python
p = flopy.modpath.ParticleData(
    partlocs=particle_locations,
    structured=True,
    localx=0.5, localy=0.5, localz=0.5,
    timeoffset=0., drape=0
)
pg1 = flopy.modpath.ParticleGroup(
    particlegroupname='PG1',
    particledata=p,
    releasedata=0.0
)
```

### 3.4 Simulatieobject aanmaken
```python
mpsim = flopy.modpath.Modpath7Sim(
    mp7,
    simulationtype='pathlines',
    trackingdirection='forward',
    weaksinkoption='pass_through',
    weaksourceoption='pass_through',
    budgetoutputoption='no',
    referencetime=[0, 0, 0.0],
    stoptimeoption='extend',
    particlegroups=pg1
)
```

## 4. Run MODPATH 7
```python
mp7.write_input()
mp7.run_model()
```

## 5. Veelgemaakte fouten
- **Fout in padnamen**: MODPATH is gevoelig voor spaties en relatieve/absolute paden.
- **Verkeerde volgorde of ontbrekende packages**.
- **Grid-bestand `.dis.grb` ontbreekt of wordt niet gevonden.**

## 6. Diagram
```mermaid
graph TD
    A[MODFLOW 6 model] -->|output| B[.hds/.cbc/.dis.grb files]
    B --> C[Modpath7 object]
    C --> D[Particles: ParticleData + ParticleGroup]
    D --> E[Modpath7Sim object]
    E --> F[mp7.write_input()]
    F --> G[mp7.run_model()]
```

---

*Laat gerust weten als je aanvullingen wilt, of als je een versie in LaTeX of PDF wilt maken voor documentatie.*

