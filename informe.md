# Módulo 2B: árboles y alineamiento múltiple

Cada ejercicio se corre por separado desde la raíz del proyecto con el comando `python3 main.py <ejercicio>`, por ejemplo `python3 main.py 1-ii`. El código está en `exercises/`, las funciones compartidas en `ncbi.py` (descarga), `blast.py` (comparación todos contra todos), `matrices.py` (matrices, simetría y distancias) y `plots.py` (gráficos). Las secuencias descargadas y la salida de BLAST quedan cacheadas en `data/`; si se borran, la primera corrida vuelve a generarlas.

A diferencia del TP anterior, hay una dependencia externa: **BLAST+ instalado localmente** (`brew install blast`). La alternativa era mandar las comparaciones al servidor del NCBI con `NCBIWWW.qblast`, pero son 90 pares ordenados encolados de a uno, lo que da corridas de más de una hora y falla seguido por timeouts. Con BLAST local el all-vs-all completo es un solo comando y tarda segundos.

---

## Datos

### La familia

Usamos **beta-globina**, una entrada por especie, diez especies. La elegimos por dos motivos que tienen que ver con lo que el TP pide después:

- **La filogenia real es conocida.** El árbol guía del ejercicio 2 tiene entonces contra qué contrastarse: si no agrupa humano con chimpancé y no deja al pez cebra afuera, el problema está en la matriz de distancias, no en la biología. Un conjunto de secuencias sin filogenia conocida no permite distinguir "mi algoritmo funciona" de "mi algoritmo devolvió algo".
- **La divergencia es gradual.** La otra manera obvia de conseguir diez secuencias emparentadas es tomar los diez primeros hits de un BLAST, pero esos vienen ordenados por puntaje, así que tienden a ser los diez parientes *más cercanos* y la matriz sale casi uniforme. Con una matriz uniforme cada decisión del clustering es un empate resuelto por el orden de las filas.

El árbol de especies que estas diez deberían recuperar es:

```
(((((humano, chimpancé), macaco), ((ratón, rata), conejo)), vaca), pollo, rana, pez cebra)
```

primates y roedores como clados compactos, el pez en la raíz.

### El primer intento y por qué falló

La primera versión de `ncbi.py` no fijaba las entradas: buscaba, por especie, con `esearch` sobre la base `protein` del NCBI restringida a Swiss-Prot, y se quedaba con el primer resultado. La consulta estricta era

```
"hemoglobin subunit beta"[Protein Name] AND "<especie>"[Organism] AND swissprot[filter]
```

y, cuando esa no devolvía nada, caía a una laxa:

```
hemoglobin beta[All Fields] AND "<especie>"[Organism] AND swissprot[filter]
```

**El fallback laxo arruinó cuatro de las diez.** La consulta estricta falla en varias especies porque `[Protein Name]` es un campo de frase exacta y Swiss-Prot no las nombra igual en todos lados: en ratón la entrada se llama "Hemoglobin subunit **beta-1**". Al caer al `[All Fields]`, que busca en el texto completo de la entrada, bajamos:

| id | qué bajó realmente | largo |
|---|---|---|
| `Mus_musculus` | `NFE2_MOUSE`, factor de transcripción NF-E2 | 373 aa |
| `Rattus_norvegicus` | `NFE2_RAT`, factor de transcripción NF-E2 | 373 aa |
| `Oryctolagus_cuniculus` | `HPT_RABIT`, haptoglobina | 347 aa |
| `Danio_rerio` | `ADN1A_DANRE`, proteína neuroprotectora dependiente de actividad | 969 aa |

Las tres primeras tienen alguna relación con la hemoglobina —NF-E2 regula el locus de la beta-globina, la haptoglobina une hemoglobina libre en plasma— y por eso la mencionan en su anotación. La cuarta ni eso: le pegó por alguna referencia bibliográfica.

Nada se quejó. El error apareció recién tres pasos después, en la salida de BLAST, y no como un error sino como resultados raros:

```
Oryctolagus_cuniculus   Gallus_gallus    80.000    5      <- 80% de identidad sobre 5 residuos
Oryctolagus_cuniculus   Danio_rerio      71.429    7
Mus_musculus            (ningún hit contra ninguna globina)
```

Esa primera línea es, por accidente, la mejor ilustración del problema de `pident` que discutimos más abajo: haptoglobina contra beta-globina de pollo da 80% de identidad, y con una matriz armada sobre `pident` ese par se uniría en el árbol antes que humano-vaca.

Sacamos dos conclusiones que están implementadas en la versión final:

1. **Una búsqueda automática no puede elegir por nosotros**, y no por un problema de nomenclatura sino de biología: la familia es paralóga (ver abajo).
2. **Un dato malo tiene que hacer ruido en el momento en que entra**, no tres etapas después. `ncbi.py` ahora valida cada descarga y corta la corrida si algo no cierra.

### La elección final

Listamos los candidatos con la API REST de UniProt, que devuelve tabla plana:

```
https://rest.uniprot.org/uniprotkb/search
  ?query=(gene:hbb OR protein_name:"hemoglobin subunit beta")
         AND reviewed:true AND (organism_id:9606 OR organism_id:9598 OR ...)
  &fields=accession,id,protein_name,length,organism_name
  &format=tsv
```

y elegimos a mano. Las diez entradas quedaron fijas en el código (`FAMILY`, en `ncbi.py`):

| Especie | Accession | Entry name | Largo | Por qué esta y no sus parálogos |
|---|---|---|---|---|
| *Homo sapiens* | P68871 | `HBB_HUMAN` | 147 | única beta adulta |
| *Pan troglodytes* | P68873 | `HBB_PANTR` | 147 | única beta adulta |
| *Macaca mulatta* | P02026 | `HBB_MACMU` | 146 | única beta adulta |
| *Mus musculus* | P02088 | `HBB1_MOUSE` | 147 | beta-major, adulta; sobre beta-2 / H0 / H1 / epsilon-Y2 |
| *Rattus norvegicus* | P02091 | `HBB1_RAT` | 147 | major-form, adulta; sobre beta-2 minor |
| *Oryctolagus cuniculus* | P02057 | `HBB_RABIT` | 147 | beta-1/2, única beta adulta revisada |
| *Bos taurus* | P02070 | `HBB_BOVIN` | 145 | única beta adulta |
| *Gallus gallus* | P02112 | `HBB_CHICK` | 147 | beta adulta; sobre las embrionarias |
| *Xenopus laevis* | P02132 | `HBB1_XENLA` | 146 | beta-major, adulta; sobre beta-2 y la larval |
| *Danio rerio* | Q90486 | `HBB1_DANRE` | 148 | beta-A1, adulta; sobre beta-A2 |

Fijar accessions en vez de fijar una consulta tiene una ventaja concreta más allá de haber sido lo que arregló el problema: una búsqueda re-corrida dentro de seis meses puede devolver otra cosa, diez accessions no. Es la diferencia entre un método descrito y un método reproducible.

### Paralogía: por qué había que elegir a mano

El punto que hizo imposible automatizar la elección no es que las entradas se llamen distinto, es que **hay varias beta-globinas por especie y no son intercambiables**. Solo el ratón tiene cinco revisadas:

| Accession | Nombre |
|---|---|
| P02088 | beta-1 (beta-major) |
| P02089 | beta-2 (beta-minor) |
| P04443 | beta-H0 |
| P04444 | beta-H1 |
| P02104 | epsilon-Y2 |

Rata, *Xenopus* y pez cebra tienen dos cada una. El cluster de beta-globina se duplicó de manera **independiente** en roedores, en anfibios y en peces, así que "la beta-globina del ratón" no es una pregunta que una query pueda responder: es una decisión. La convención que adoptamos es la **cadena adulta mayoritaria**, que UniProt anota explícitamente ("beta-major chain", "major-form").

Esto deja una salvedad que conviene tener presente al leer el árbol del ejercicio 2: **el árbol génico no es exactamente el árbol de especies**. Como las duplicaciones son posteriores a varias de las separaciones entre linajes, la beta-1 del ratón y la beta-1 del pez cebra no son ortólogos estrictos 1:1, sino "la que quedó como adulta principal en cada linaje". En los nodos superficiales (primates, roedores) esto no afecta nada; en los profundos (pez contra tetrápodos) puede desviar el árbol respecto de la filogenia de especies. Si el árbol final falla, es el primer lugar donde mirar antes de sospechar del algoritmo.

### Validación de la descarga

`ncbi.py` valida cada entrada al bajarla y **corta la corrida** si alguna no pasa:

- la descripción tiene que contener "hemoglobin" o "globin";
- el largo tiene que caer entre 120 y 200 aa.

La regla de largo sola habría descartado las cuatro entradas equivocadas del primer intento (373, 373, 347 y 969 aa). Son dos chequeos porque los dos modos de falla son distintos: un accession equivocado da la proteína equivocada con el largo equivocado, mientras que un accession correcto sobre un precursor da la proteína correcta con un largo de más.

---

## Ejercicio 1

> **Ejercicio 1)** Utilice lo aprendido en los trabajos prácticos previos (uso de biopython para obtener secuencias, BLAST y el "parser del mismo") para: i) obtener 10 secuencias pertenecientes a una misma familia y ii) alinearlas y obtener la correspondiente matriz de puntaje / identidad. Analice si la matriz es (o no) simétrica. Compare matrices utilizando diferentes puntajes que le otorga BLAST (Bit Score, Identidad, E-Value).

### i) Las diez secuencias

`python3 main.py 1-i`

Ver la sección **Datos**. La salida del comando imprime, por especie, accession, entry name, largo y el motivo de la elección, que es la misma información de la tabla de arriba pero generada desde `data/globins.json`, o sea desde lo que realmente se descargó y no desde lo que quisimos descargar.

**[PENDIENTE: pegar la salida de `1-i` después de rebajar los datos.]**

### ii) La matriz de puntaje / identidad

`python3 main.py 1-ii`

El all-vs-all es un solo comando:

```bash
blastp -query data/globins.fasta -subject data/globins.fasta \
       -outfmt "6 qseqid sseqid pident length ... nident qlen slen" \
       -evalue 10000 -max_hsps 1
```

Cuatro decisiones sobre ese comando, todas con consecuencias en la matriz:

**Un `-subject` en vez de `makeblastdb`.** Con un archivo de 10 como query y el mismo archivo como subject, BLAST hace las 100 comparaciones ordenadas en una sola llamada, incluyendo los self-hits.

**Los self-hits se conservan.** Son la diagonal contra la que hay que normalizar el bit score: un bit score de 200 no significa nada hasta saber que la secuencia contra sí misma puntúa 210.

**`-evalue 10000` en vez del default 10.** El corte por defecto descarta en silencio los pares más distantes. Un par faltante no es "distancia infinito", es un agujero en la matriz que el clustering del ejercicio 2 tendría que rellenar inventando un número; subir el corte mueve esa invención del código a la estadística de BLAST. `1-ii` imprime explícitamente los pares que aun así no alinearon.

**`-max_hsps 1`.** Sin esto, un par con dos alineamientos locales aporta dos filas y el parser tiene que decidir cuál es *la* comparación. En la corrida fallida del primer intento esto se vio bien: `Oryctolagus × Danio` ocupaba cuatro filas y `Gallus × Oryctolagus` dos. Vale la pena notar que la *forma* del resultado ya delata la falta de homología: dos proteínas no emparentadas que comparten pedacitos por azar dan varios HSPs cortos y dispersos, en lugar de uno solo largo.

#### Qué reporta BLAST realmente: el HSP

BLAST es un alineador **local** (una heurística de Smith-Waterman), así que nunca devuelve "el alineamiento de estas dos secuencias". Devuelve uno o más **HSP** (*High-scoring Segment Pair*): pares de tramos, uno de cada secuencia, tales que el puntaje no mejora ni extendiendo ni recortando. Los bordes del HSP no son un parámetro: emergen de la extensión con X-drop, que estira desde una semilla mientras el puntaje acumulado no caiga más de X por debajo del máximo alcanzado, y después **retrocede al máximo**. Por eso un HSP puede empezar y terminar en el medio de la secuencia.

Esto importa para el TP porque el HSP es una decisión que BLAST toma con criterio de *búsqueda* ("¿vale la pena reportar esto?") y no con criterio *filogenético* ("¿cuánto divergieron estas dos secuencias?"). Son objetivos distintos, y la matriz de distancias tiene que corregir por esa diferencia.

#### `pident` no es la identidad que queremos

La columna `pident` de BLAST es

$$\text{pident} = \frac{\text{nident}}{\text{length del HSP}}$$

o sea *de la parte que alineó, cuánto matcheó*. Lo que quedó fuera del HSP no aparece en ninguno de los dos términos. Dos secuencias que comparten un motivo conservado de 30 residuos y nada más reportan 95%.

Agregamos entonces una segunda identidad, que en el código llamamos `gident`:

$$\text{gident} = \frac{\text{nident}}{\min(\text{qlen}, \text{slen})}$$

que mete la cobertura adentro: lo que no alineó ahora cuenta como no-match. Coinciden exactamente cuando el alineamiento es de punta a punta. No hay re-alineamiento: `gident` sale de la misma fila de BLAST, y por eso el `-outfmt 6` pide `nident`, `qlen` y `slen`, que no vienen en el conjunto por defecto.

Dividir por la más corta es una convención entre varias; dividir por la más larga es más estricto y castiga que una de las dos tenga cola extra. Elegimos `min` por ser la más indulgente de las que igual incorporan cobertura: no penaliza que un ortólogo tenga un péptido señal de más.

Las dos identidades viajan hasta el final a propósito. La consigna pide comparar puntajes, y comparar dos cosas que *se llaman igual* y dan árboles distintos es más interesante que comparar Bit Score contra E-value.

**[PENDIENTE: pegar las cuatro matrices crudas de `1-ii`.]**

### ¿Es simétrica la matriz?

`python3 main.py 1-sym`

La respuesta corta es **depende de la métrica**, y la razón no es numérica sino estructural.

**Identidad y bit score son propiedades del alineamiento.** El alineamiento de (i, j) es esencialmente el de (j, i). Lo poco que puede quedar de asimetría viene de las *composition-based statistics*, que reescalan la matriz de sustitución usando la composición de aminoácidos del **query**, y de la heurística de semillas, que indexa el subject y recorre el query.

**El E-value es distinto en especie.** Vale

$$E = K \cdot m \cdot n \cdot e^{-\lambda S}$$

con $m$ el largo del query y $n$ el tamaño de la base. Intercambiar query y subject intercambia $m$ y $n$, así que salvo que las dos secuencias midan exactamente lo mismo el E-value **tiene** que diferir, y la diferencia crece con la razón de largos. No es ruido, es la definición. Un E-value es una afirmación sobre una *búsqueda*, no sobre un *par*, y ese es el motivo de fondo por el que es mal candidato a distancia.

Para simetrizar promediamos las dos direcciones, porque ninguna de las dos es la correcta: "la distancia entre humano y pollo" no tiene query ni subject. El E-value se promedia **en escala logarítmica**, ya que es una probabilidad que abarca doscientos órdenes de magnitud y $(10^{-100} + 10^{-3})/2$ es $10^{-3}$ con pasos de más.

**[PENDIENTE: pegar la tabla de simetría de `1-sym` y los pares con mayor discrepancia.]**

### Comparación entre los distintos puntajes

`python3 main.py 1-cmp`

Cada métrica necesita su propia conversión a distancia, y equivocarse acá es la manera más silenciosa de terminar con un árbol hermoso y al revés:

| Métrica | Conversión | Por qué |
|---|---|---|
| Identidad | $d = 1 - \text{id}$ | La p-distance de manual. Satura: pasado ~75% de divergencia todos los pares parecen igual de lejanos, y por eso la filogenia real le aplica encima una corrección de Poisson o de Kimura. |
| Bit score | $d = 1 - S_{ij} / \min(S_{ii}, S_{jj})$ | El bit score no es una similitud en [0,1]: crece con el largo, así que un par largo y distante puede superar a uno corto y casi idéntico. El self-score es lo máximo que el par podría haber puntuado. |
| E-value | $-\log_{10} E$ reescalado a [0,1] | Monótono en la similitud, así que sirve como *ranking*, pero no es una distancia: depende del tamaño de la base, con lo cual los mismos dos pares cambian de distancia cuando se agrega una secuencia al conjunto, y satura en 0.0 para todo lo cercano, colapsando los pares próximos en un mismo valor. |

La comparación se reporta en tres formas, de menos a más consecuente:

1. **Correlación de rangos (Spearman)** entre las cuatro matrices de distancia, sobre los 45 pares. De rangos y no de Pearson porque el árbol nunca mira los valores, solo cuál par es el más chico: dos matrices que difieren en todos los números pero coinciden en el orden construyen el mismo árbol.
2. **Tabla de vecino más cercano** por secuencia. Es la consecuencia práctica: dos métricas pueden correlacionar en 0,99 y aun así arrancar el árbol distinto.
3. **El primer par que uniría el árbol guía** según cada métrica, que es literalmente el paso i) del procedimiento del ejercicio 2.

Los gráficos que genera son `exercises/ex1_distance_matrices.png` (las cuatro matrices en escala común, para que se vea que el mismo conjunto se ve más o menos disperso según el puntaje) y `exercises/ex1_metric_agreement.png` (cada métrica contra `gident`, un punto por par; lo que importa no es la dispersión sino si el *orden* se preserva).

**[PENDIENTE: pegar la matriz de Spearman, la tabla de vecinos más cercanos y las dos figuras.]**

### Nota: la sensibilidad de BLAST y la matriz de sustitución

Corrimos también el all-vs-all con `-matrix BLOSUM45` para ver cuánto de la matriz de distancias es propiedad de las secuencias y cuánto del sistema de puntaje. BLOSUM45 tolera más divergencia y produce HSPs más largos; BLOSUM80 corta antes y devuelve solo el núcleo conservado.

En la corrida con los datos malos ya se vio un efecto claro y en la dirección esperada: BLOSUM45 armó un HSP de **115 residuos al 27% de identidad** entre la beta-globina humana y una proteína neuroprotectora sin ningún parentesco. Ése es el costo de una matriz pensada para la zona crepuscular: gana sensibilidad y fabrica alineamientos de ruido.

**[PENDIENTE: rehacer la comparación BLOSUM62 / BLOSUM45 con los datos corregidos, comparando la columna `length`. En los pares cercanos no debería moverse; entre humano y pez cebra probablemente sí.]**

Vale la pena dejar anotado el límite de fondo, porque explica los pares sin hit: BLAST no examina todos los alineamientos posibles. Solo extiende desde donde encontró una **semilla** —dos palabras de 3 residuos, en la misma diagonal, que puntúen ≥ 11 en la matriz de sustitución—, así que una región homóloga sin semilla no se descarta: no se mira nunca. Es el precio de no correr Smith-Waterman completo, y en la práctica casi no se paga, porque los homólogos reales conservan residuos idénticos salpicados. Se paga por debajo de ~25% de identidad, que es justo donde uno más querría una respuesta, y es el motivo por el que existen PSI-BLAST y los métodos basados en HMM.

Entonces, si `1-ii` reporta algún par sin hit —el candidato natural es pez cebra contra pollo o contra rana—, eso **no** significa distancia infinita: significa que BLAST no encontró por dónde empezar. Es un agujero de método, no un dato biológico. Las dos salidas honestas son bajar el umbral de semilla con `-threshold 9` (vecindario más grande, más lento) o completar ese par con un Smith-Waterman de verdad vía `PairwiseAligner` de Biopython, dejándolo anotado.

---

## Ejercicio 2

> **Ejercicio 2)** Programe un algoritmo que dada una matriz de distancias construya un árbol guía. Utilice el mismo para trabajar con las diferentes matrices obtenidas a partir de los puntajes de BLAST y analice los resultados de manera comparativa.

`python3 main.py 2a`

### El algoritmo

Está en `tree.py`, y es el procedimiento del enunciado tal cual: buscar el par más cercano, unirlo en un nodo, recalcular la distancia de todos los demás al nodo nuevo como el promedio de las distancias a sus dos hijos, repetir. Nueve fusiones para diez secuencias.

Una observación sobre la implementación: **las opciones a) y b) del enunciado no son dos ramas del código**. Un nodo ya creado *es* un elemento de la matriz nueva, con la misma jerarquía que una secuencia suelta, así que "buscar el par más cercano" cubre los dos casos sin distinguirlos. La distinción del enunciado es sobre cómo se dibuja el árbol resultante, no sobre cómo se decide en cada paso.

El promedio simple que pide el enunciado es **WPGMA**, no UPGMA: pondera igual a los dos hijos sin importar cuántas secuencias tenga cada uno. La diferencia aparece recién cuando se fusionan dos nodos de tamaños distintos.

Las alturas siguen la convención ultramétrica: un nodo que une dos elementos a distancia $d$ queda a altura $d/2$, y la rama hacia cada hijo es lo que falta desde la altura del hijo hasta la del padre. Las hojas están a altura 0.

Dos detalles que no son features sino correctitud:

- **Los empates se desempatan por nombre.** Humano y chimpancé tienen la beta-globina idéntica, o sea distancia exactamente 0, así que el empate existe de verdad; sin un criterio fijo, dos corridas con la misma matriz podrían dar árboles distintos según el orden de iteración.
- **Las ramas negativas se recortan a 0.** Pueden aparecer si las distancias no son ultramétricas, es decir si un hijo queda más alto que su padre. Es síntoma de que los datos no cumplen el reloj molecular que WPGMA asume, no de un error de cálculo.

### La figura: árbol y matriz juntos

Cada corrida produce una figura de dos paneles: el árbol a la izquierda, dibujado con `Phylo.draw` con las ramas a escala, y la matriz de distancias a la derecha con los valores en cada celda.

**Las filas de la matriz están reordenadas según el orden de las hojas del árbol**, no según el orden original de `ids`. Ése es el punto de ponerlas al lado: cada fila queda a la altura de su hoja, y entonces las dos mitades de la figura son la misma información en dos lenguajes. Los bloques claros sobre la diagonal son los clados.

### Cómo leer el heatmap como si fuera el árbol

Con las filas ordenadas por el árbol, la **subdiagonal** —la distancia entre cada par de hojas adyacentes— cuenta la topología. En la corrida con `gident`:

```
Xenopus–Danio       .575
Danio–Gallus        .435
Gallus–Mus          .347
Mus–Rattus          .082   <- pozo
Rattus–Bos          .248
Bos–Oryctolagus     .179
Oryctolagus–Macaca  .110
Macaca–Homo         .055
Homo–Pan            .000
```

Decreciente en todo salvo en un lugar: `.347 → .082 → .248`.

Si el árbol fuera una **escalera pura** —cada secuencia colgándose del grupo que ya viene creciendo, la opción b) del enunciado en todos los pasos— la subdiagonal saldría monótona decreciente, porque cada hoja adyacente se uniría a una altura menor que la anterior. Un valor que baja y vuelve a subir es la firma de una **cereza**: dos hojas vecinas que se emparejan entre sí a una altura mucho más superficial que la estructura que las rodea. Eso es exactamente la opción a): dos elementos más cercanos entre sí que de cualquier nodo ya creado.

Y el `.248` que viene después no es un retroceso: es la altura a la que ese par ya formado se engancha con el resto, por eso salta.

Hay una segunda cereza, `Homo–Pan` en `.000`, que no se ve como pozo porque cae al final de la secuencia y no queda nada después con qué contrastar.

Vale aclarar que esto **no es un atajo para construir el árbol**, es un diagnóstico visual del árbol ya construido: la subdiagonal se lee así justamente porque las filas ya se ordenaron con el resultado. En el orden original de `ids` la adyacencia sería arbitraria y los pozos no significarían nada.

### Las cuatro métricas

El algoritmo se corre sobre las cuatro matrices de distancia del ejercicio 1 (`pident`, `gident`, `bitscore`, `evalue`), produciendo `exercises/ex2_tree_<métrica>.png` y el Newick de cada una.

**[PENDIENTE: pegar los cuatro Newick, las cuatro figuras, y la comparación entre topologías: qué clados recupera cada una y cuáles coinciden con el árbol de especies esperado.]**

---

## Ejercicio 3

> **Ejercicio 3)** Construya un algoritmo que dados los AP entre las secuencias del objetivo 1 y el árbol guía genere el AM y lo muestre en pantalla. Compare sus resultados con los obtenidos utilizando algún programa on-line.

**[PENDIENTE]**
