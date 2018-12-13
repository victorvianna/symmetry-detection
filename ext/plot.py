from matplotlib import pyplot

x_id = 3
y_id = 4
for filename in ["complete_transf_space.txt", "pruned_transf_space.txt", "final_transf_space.txt"]:
    file = open(filename)
    if file:
        x = []
        y = []
        for line in file:
            transf = [float(c) for c in line.split(" ")[:-1]]
            x.append(transf[x_id])
            y.append(transf[y_id])
        pyplot.scatter(x,y)
pyplot.savefig('transformations.png')
