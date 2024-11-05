import math
import random
from bacteria import bacteria

class chemiotaxis():
    def __init__(self):
       parcialNFE = 0  
    
    def compute_cell_interaction(self, bacteria, poblacion, d, w):
      total = 0.0
      for other in poblacion:
        diff = 0.0
        diff += (bacteria.blosumScore - other.blosumScore) ** 2.0
        total += d * math.exp(w * diff)
      return total

    def attract_repel(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
      attract = self.compute_cell_interaction(bacteria, poblacion, -d_attr, -w_attr)
      repel = self.compute_cell_interaction(bacteria, poblacion, h_rep, -w_rep)
      return attract + repel   #interaction
    
    def chemio(self, bacteria, poblacion, d_attr, w_attr, h_rep, w_rep):
      bacteria.interaction = self.attract_repel(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
      bacteria.fitness = bacteria.blosumScore + bacteria.interaction
      
    def doChemioTaxis(self, poblacion, d_attr, w_attr, h_rep, w_rep):
      self.parcialNFE = 0
      for bacteria in poblacion:
        self.chemio(bacteria, poblacion, d_attr, w_attr, h_rep, w_rep)
        self.parcialNFE += bacteria.NFE
        bacteria.NFE = 0
        
    def eliminarClonar(self, path, poblacion):
      """Elimina el 50% de las bacterias con menos fitness."""
      poblacion.sort(key=lambda x: x.fitness)
      for i in range(int(len(poblacion) / 2)):
          del poblacion[0]
          
      clones = self.clonacion(path, poblacion)
      for clone in clones:
          poblacion.append(clone)
    
    def clonacion(self, path, poblacion):
       poblacionClones = []
       best = max(poblacion, key=lambda x: x.fitness)
       fitness_threshold = 0.5  # Ajusta este valor si es necesario
       for bacteria in poblacion:
         newBacteria = bacteria.clonar(path)
         mutacion = int((best.fitness - bacteria.fitness) / 10)  # Mutación en función de la diferencia de fitness
         newBacteria.tumboNado(mutacion, fitness_threshold)
         newBacteria.autoEvalua()
         poblacionClones.append(newBacteria)
       return poblacionClones

    def randomBacteria(self, path):
       bact = bacteria(path)
       fitness_threshold = 0.5  # Ajusta este valor si es necesario
       bact.tumboNado(random.randint(1, 10), fitness_threshold)
       return bact 
   
    def insertRamdomBacterias(self, path, num, poblacion):
      for i in range(num):
         poblacion.append(self.randomBacteria(path))
         # Eliminar la bacteria con menos fitness
         poblacion.sort(key=lambda x: x.fitness)
         del poblacion[0]
