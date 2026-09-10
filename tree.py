"""De la matriz de distancias al arbol guia: el procedimiento del enunciado.

Repetidamente:
  i)  buscar los dos elementos mas cercanos y unirlos en un nodo;
  ii) recalcular la distancia de los "otros" elementos al nodo nuevo como el
      promedio de las distancias a sus dos hijos, y volver a buscar el par
      mas cercano.

Los dos casos que el enunciado distingue -- opcion a) los dos mas cercanos son
elementos sueltos, opcion b) son un elemento y un nodo ya creado -- no
necesitan ramas separadas en el codigo: un nodo ya creado ES un elemento de la
matriz nueva, asi que "buscar el par mas cercano" los cubre a los dos.
"""


class Cluster:
    """Un elemento de la matriz: al principio una secuencia, despues un nodo.

    `height` es la altura del nodo, la mitad de la distancia a la que se
    fusionaron sus dos hijos. Las hojas estan a altura 0.
    """

    def __init__(self, newick, height=0.0, leaves=None):
        self.newick = newick
        self.height = height
        self.leaves = leaves or {newick}


def build(D, ids):
    """El arbol guia en formato Newick, a partir de una matriz simetrica."""
    n = len(ids)
    active = list(range(n))
    clusters = {i: Cluster(ids[i]) for i in range(n)}
    dist = {frozenset((i, j)): float(D[i][j])
            for i in range(n) for j in range(i + 1, n)}
    next_id = n

    while len(active) > 1:
        # i) el par mas cercano. Los empates se desempatan por nombre para que
        # dos corridas con la misma matriz den el mismo arbol: con secuencias
        # identicas (humano y chimpance) la distancia es 0 y el empate existe.
        pair = min(
            (frozenset((a, b)) for k, a in enumerate(active) for b in active[k + 1:]),
            key=lambda p: (dist[p], sorted(n for k in p for n in clusters[k].leaves)),
        )
        a, b = sorted(pair, key=lambda x: sorted(clusters[x].leaves))
        height = dist[pair] / 2

        # Las ramas hacia cada hijo son lo que falta desde la altura del hijo
        # hasta la del padre. El max(..., 0) cubre el caso de distancias no
        # ultrametricas, donde un hijo puede quedar mas alto que el padre.
        la = max(height - clusters[a].height, 0.0)
        lb = max(height - clusters[b].height, 0.0)

        clusters[next_id] = Cluster(
            newick=f"({clusters[a].newick}:{la:.4f},{clusters[b].newick}:{lb:.4f})",
            height=height,
            leaves=clusters[a].leaves | clusters[b].leaves,
        )

        # ii) recalcular la distancia del nodo nuevo al resto.
        active.remove(a)
        active.remove(b)
        for other in active:
            dist[frozenset((next_id, other))] = (
                dist[frozenset((a, other))] + dist[frozenset((b, other))]) / 2
        active.append(next_id)
        next_id += 1

    return clusters[active[0]].newick + ";"
