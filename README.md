# TP 3: árboles y alineamiento múltiple

**El informe está en `informe.pdf`.** Es autocontenido, con las figuras adentro. También está el mismo texto en `informe.md`, que es la fuente del pdf (`python3 build_informe.py` lo regenera, junto con `informe.html`).

## Cómo correrlo

Hace falta Python 3 con `biopython`, `numpy`, `matplotlib` y `markdown`, y dos programas externos:

```
brew install blast mafft
```

Cada ejercicio se corre por separado, con el nombre que tiene en la consigna:

```
python3 main.py 1-i      # 10 secuencias de una misma familia (beta-globina)
python3 main.py 1-ii     # all-vs-all con BLAST: las cuatro matrices
python3 main.py 1-sim    # ¿son simétricas?
python3 main.py 1-cmp    # comparación entre bit score, identidad y E-value
python3 main.py 1-mat    # nota: BLOSUM45 / 62 / 80
python3 main.py 2a       # árbol guía de cada matriz, y su figura
python3 main.py 2b       # comparación de las topologías
python3 main.py 3a       # alineamiento múltiple, mostrado en pantalla
python3 main.py 3b       # comparación con MAFFT
```

Las secuencias, la salida de BLAST y los alineamientos de MAFFT quedan cacheados en `data/`, así que con `data/` completo todo corre sin red. Si se borra, la primera corrida vuelve a descargarlo y a correr BLAST y MAFFT.

## Dónde está cada cosa

| Archivo | Qué tiene |
|---|---|
| `sequences.py` | descarga de las 10 secuencias, con las entradas fijas y validadas, y el all-vs-all con `blastp` |
| `matrices.py` | matrices de puntaje, simetría, conversión a distancias y árbol guía |
| `msa.py` | del árbol y los alineamientos de a pares al alineamiento múltiple |
| `plots.py` | figuras (se guardan en `figures/`) |
| `exercises/` | un archivo por ejercicio |
| `build_informe.py` | genera `informe.html` e `informe.pdf` a partir de `informe.md` |
