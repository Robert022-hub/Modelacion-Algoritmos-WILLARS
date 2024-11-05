from bacteria import bacteria
from chemiotaxis import chemiotaxis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parámetros de configuración
path = "multifasta.fasta"
numeroDeBacterias = 5
numRandomBacteria = 1
iteraciones = 30
tumbo = 1  # Número de gaps a insertar
nado = 3
chemio = chemiotaxis()
veryBest = bacteria(path)  # Mejor bacteria
globalNFE = 0  # Número de evaluaciones de la función objetivo

# Parámetros de quimiotaxis
dAttr = 0.1
wAttr = 0.2
hRep = dAttr
wRep = 10

# Inicializar población
poblacion = [bacteria(path) for _ in range(numeroDeBacterias)]

# Variables para almacenar los datos de cada iteración
resultados = []

# Función para clonar la mejor bacteria
def clonaBest(veryBest, best):
    veryBest.matrix.seqs = np.array(best.matrix.seqs)
    veryBest.blosumScore = best.blosumScore
    veryBest.fitness = best.fitness
    veryBest.interaction = best.interaction

# Proceso de iteración
for iteracion in range(iteraciones):
    for bacterium in poblacion:
        fitness_threshold = 0.5
        bacterium.tumboNado(tumbo, fitness_threshold)
        bacterium.autoEvalua()
    
    chemio.doChemioTaxis(poblacion, dAttr, wAttr, hRep, wRep)
    globalNFE += chemio.parcialNFE
    best = max(poblacion, key=lambda x: x.fitness)
    
    if (veryBest is None) or (best.fitness > veryBest.fitness):
        clonaBest(veryBest, best)
    
    # Guardar resultados de esta iteración
    resultados.append({
        "Iteración": iteracion + 1,
        "Mejor Fitness": veryBest.fitness,
        "Blosum Score": veryBest.blosumScore,
        "Interacción": veryBest.interaction,
        "NFE": globalNFE,
        "Tamaño de Población": len(poblacion)
    })

    chemio.eliminarClonar(path, poblacion)
    chemio.insertRamdomBacterias(path, numRandomBacteria, poblacion)

# Mostrar el genoma de la mejor bacteria final
veryBest.showGenome()

# Convertir resultados a DataFrame
df = pd.DataFrame(resultados)

# Mostrar la tabla de resultados
print(df)

# Graficar los resultados
plt.figure(figsize=(12, 6))
plt.plot(df["Iteración"], df["Mejor Fitness"], label="Mejor Fitness", marker='o')
plt.plot(df["Iteración"], df["Blosum Score"], label="Blosum Score", marker='x')
plt.plot(df["Iteración"], df["Interacción"], label="Interacción", marker='s')
plt.xlabel("Iteración")
plt.ylabel("Valor")
plt.title("Evolución del Fitness y la Interacción por Iteración")
plt.legend()
plt.grid(True)
plt.show()
