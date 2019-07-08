# Examen Final CC471-2019-1 
# Problema N° 1b 
# Autor: Lazaro Camasca Edson
# Nombre del Programa: crear_Arbol.py

#========================================================================
# 1.a) Crear el arbol filogenetico

from Bio import AlignIO

def leer_SecAli(input_clustal = 'data_gen/Sec_Alineadas.clustal'):

    """ Recibe la direccion del archivo clutal con la secuencias aliendas 
        Retorna un objeto alineamiento"""
    
    aln = open(input_clustal, "r")
    # AlignIO para leer el archivo de alineamiento en formato 'clustal' format
    alignment = AlignIO.read(aln, "clustal")
    return alignment

aln = leer_SecAli()
print("\nLas secuencias alineadas son\n",aln)

from Bio.Phylo.TreeConstruction import DistanceCalculator  # crear la matriz de distancias
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

calculator = DistanceCalculator('blosum62')
# añade la matriz de  distancias al objeto calculator y lo retorna
dm = calculator.get_distance(aln)
print(dm)


def create_Tree(alignment, tipo = 'upgma'):
    """Se recibe de entrada el alimenamiento y el tipo de arbol nj o upgma 
    Genera ell arbol filogenetico 
    """
    print('Creando el arbol filogenetico ......')
    # 3. Creamos la matriz de distancias
    calculator = DistanceCalculator('identity')
    # añade la matriz de  distancias al objeto calculator y lo retorna
    dm = calculator.get_distance(alignment)

    #initialize a DistanceTreeConstructor object based on our distance calculator object
    constructor = DistanceTreeConstructor(calculator)

    #De acuerdo a las opciones elegimos el tipo de arbol
    #build the tree
    if tipo == 'upgma':
        tree = constructor.upgma(dm)
    if tipo == 'nj':
        tree = constructor.nj(dm)

    return tree

tree = create_Tree(aln)

from Bio import Phylo
import pylab

#Mostrar el arbol
def draw_Tree(tree):
    """ Muestra el arbol en pantalla"""
    Phylo.draw(tree)

draw_Tree(tree)


tree = create_Tree(aln, 'nj')

draw_Tree(tree)