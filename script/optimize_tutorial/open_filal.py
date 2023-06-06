import pickle
# You may want to import your own packages if the pickle file contains custom objects

with open("final.p", "rb") as f:
    data = pickle.load(f)
# ``data`` is now a dictionary containing all results, including the final container, all solutions, the algorithm parameters, etc.

grid = data['container']
print(grid.best)
print(grid.best.fitness)
print(grid.best.features)
