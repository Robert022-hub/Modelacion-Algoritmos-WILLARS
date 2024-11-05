from fastaReader import fastaReader
import random
import numpy 
import copy
from evaluadorBlosum import evaluadorBlosum

class bacteria():

    def __init__(self, path):
        self.matrix = fastaReader(path)
        self.blosumScore = 0
        self.fitness = 0
        self.interaction = 0
        self.NFE = 0
        self.attraction = 1.0  # Nueva variable para atracción adaptativa
        self.repulsion = 1.0   # Nueva variable para repulsión adaptativa

    def showGenome(self):
        for seq in self.matrix.seqs:
            print(seq)

    def clonar(self, path):
        newBacteria = bacteria(path)
        newBacteria.matrix.seqs = numpy.array(copy.deepcopy(self.matrix.seqs))
        return newBacteria    

    def tumboNado(self, numGaps, fitness_threshold):

        self.cuadra()
        matrixCopy = copy.deepcopy(self.matrix.seqs)
        matrixCopy = matrixCopy.tolist()
        
        # Ajuste adaptativo: si el fitness es alto, se reduce la cantidad de gaps
        if self.fitness > fitness_threshold:
            gapRandomNumber = random.randint(0, max(1, numGaps//2))
            self.attraction += 0.1  # Incrementa atracción si el fitness mejora
            self.repulsion -= 0.1   # Reduce repulsión
        else:
            gapRandomNumber = random.randint(0, numGaps)
            self.attraction -= 0.1  # Reduce atracción si el fitness no mejora
            self.repulsion += 0.1   # Incrementa repulsión

        # Inserta gaps adaptativamente
        for i in range(gapRandomNumber):
            seqnum = random.randint(0, len(matrixCopy)-1)
            pos = random.randint(0, len(matrixCopy[0]))
            part1 = matrixCopy[seqnum][:pos]
            part2 = matrixCopy[seqnum][pos:]
            temp = "-".join([part1, part2])     
            matrixCopy[seqnum] = temp
        
        self.matrix.seqs = numpy.array(matrixCopy)
        self.cuadra()
        self.limpiaColumnas()


        # Evaluación adaptativa
        self.autoEvalua()
    def cuadra(self):
        seq = self.matrix.seqs
        maxLen = len(max(seq, key=len))
        for i in range(len(seq)):
            if len(seq[i]) < maxLen:
                seq[i] = seq[i] + "-" * (maxLen - len(seq[i]))
        self.matrix.seqs = numpy.array(seq)

    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True

    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteColumn(i)
            else:
                i += 1

    def deleteColumn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]

    def getColumn(self, col):
        column = []
        for i in range(len(self.matrix.seqs)):
            column.append(self.matrix.seqs[i][col])
        return column

    def autoEvalua(self):   
        evaluador = evaluadorBlosum()
        score = 0
        for i in range(len(self.matrix.seqs[0])):
            column = self.getColumn(i)
            gapCount = column.count("-")
            column = [x for x in column if x != "-"]
            pares = self.obtener_pares_unicos(column)
            for par in pares:
                score += evaluador.getScore(par[0], par[1])
            score -= gapCount * 2
        self.blosumScore = score
        self.NFE += 1

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)