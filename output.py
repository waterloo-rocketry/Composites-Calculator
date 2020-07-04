import matplotlib.pyplot as plt
import csv


def write_matrix(file, m, title):
    with open (file, mode='w') as csv_file:
        writer = csv.writer(csv_file, delimiter = ',', lineterminator = '\n')
        writer.writerows([[title]])
        writer.writerows(m)

def plotdata(xlabel, ylabel, title, data, file_type):
    plt.plot(data)
    plt.gca().invert_yaxis()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(axis="both")
    plt.title(title)
    plt.savefig(title+"."+file_type, format = file_type)
    plt.clf()
    

