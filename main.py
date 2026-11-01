"""Entry point for the TP. Runs one exercise at a time, by name:

    python3 main.py 1-i
    python3 main.py 1-sim
    python3 main.py          # lists the available exercises

The names match the numbering in the TP PDF, so the command you type, the
section in the report and the file on disk all agree.

Needs BLAST+ on the PATH for anything past 1-i (`brew install blast`), and
MAFFT for the comparison of 3b (`brew install mafft`).
"""

import sys

from exercises import ex1, ex2, ex3

EXERCISES = {
    "1-i":   ex1.i,          # 10 secuencias de una misma familia
    "1-ii":  ex1.ii,         # alinearlas y obtener la matriz de puntaje
    "1-sim": ex1.symmetry,   # analice si la matriz es (o no) simetrica
    "1-cmp": ex1.compare,    # compare Bit Score, Identidad y E-value
    "1-mat": ex1.matrix_sensitivity,   # nota: BLOSUM62 vs BLOSUM45 vs BLOSUM80

    "2a":    ex2.a,          # construir el arbol guia y graficarlo
    "2b":    ex2.b,          # comparar las topologias entre metricas y con las especies

    "3a":    ex3.a,          # del arbol y los AP al AM, mostrado en pantalla
    "3b":    ex3.b,          # compararlo con un programa externo (MAFFT)
}


def usage():
    print("usage: python3 main.py <exercise>")
    print("available exercises:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' is not a known exercise.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
