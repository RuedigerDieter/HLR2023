
import matplotlib.pyplot as plt
import numpy as np
import os

# get all .dat files in current directory
files = [f for f in os.listdir('.') if os.path.isfile(f) and f.endswith(".dat")]
#for every dat file, create array
communication_a_gs = []
communication_a_ja = []
strong_scaling_gs = []
strong_scaling_ja = []
weak_scaling_gs = []
weak_scaling_ja = []
names = [
    "COMMUNICATION_A_GS",
    "COMMUNICATION_A_JA",
    "STRONG_SCALING_GS",
    "STRONG_SCALING_JA",
    "WEAK_SCALING_GS",
    "WEAK_SCALING_JA"
]

for f in files:
    with open(f, 'r') as file:
        # read all lines
        lines = file.readlines()
        for l in lines:
            if l.startswith("#"):
                continue
            # split line at whitespace
            split = l.split()
            procs = int(split[0])
            nodes = int(split[1])
            lines = int(split[2])
            times = float(split[3])
            if f.startswith(names[0]):
                communication_a_gs.append([procs, nodes, lines, times])
            elif f.startswith(names[1]):
                communication_a_ja.append([procs, nodes, lines, times])
            elif f.startswith(names[2]):
                strong_scaling_gs.append([procs, nodes, lines, times])
            elif f.startswith(names[3]):
                strong_scaling_ja.append([procs, nodes, lines, times])
            elif f.startswith(names[4]):
                weak_scaling_gs.append([procs, nodes, lines, times])
            elif f.startswith(names[5]):
                weak_scaling_ja.append([procs, nodes, lines, times])
            else:
                print("Error: Unknown file name " + f)
                exit(1)

overall = [communication_a_gs, communication_a_ja, strong_scaling_gs, strong_scaling_ja, weak_scaling_gs, weak_scaling_ja]

for i in range(len(overall)):
    overall[i] = np.array(overall[i])

fig = plt.figure(figsize=(10, 10))  # Create a new figure with a custom size (you can adjust the size as needed)

for j in range(len(overall)):
    ax = fig.add_subplot(3, 2, j+1)  # Create a new subplot. Adjust the numbers (3, 2) as needed to change the layout of the subplots.
    ax.set_title(names[j])
    if j == 4 or j == 5: # weak scaling
        ax.set_xlabel("Lines per Process")
        ax.set_ylabel("Time in s")
        xvals = overall[j][:,2] / overall[j][:,0]
        yvals = overall[j][:,3]
    elif j == 0 or j == 1:
        ax.set_xlabel("Number of Nodes")
        ax.set_ylabel("Time in s")
        xvals = overall[j][:,1]
        yvals = overall[j][:,3]
    else:
        ax.set_xlabel("Number of Processes")
        ax.set_ylabel("Time in s")
        xvals = overall[j][:,0]
        yvals = overall[j][:,3]

    std_dev = np.std(yvals)
    plt.axhline(y=np.mean(yvals), color='r', linestyle='--', label='Mittelwert')  # Mittelwert als gestrichelte Linie
    plt.axhline(y=np.mean(yvals) + std_dev, color='g', linestyle='--', label='Standardabweichung')  # Mittelwert + Std Dev
    plt.axhline(y=np.mean(yvals) - std_dev, color='g', linestyle='--')  # Mittelwert - Std Dev
    ax.grid(True)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    ax.plot(xvals, yvals)
    ax.legend(fontsize=5)


plt.tight_layout()  # Adjust the layout so that the subplots do not overlap
plt.savefig("plots.png")  # Save the figure with all subplots
plt.show()  # Show the figure (remove this if you don't want to see it)

#reset plot for second plot
fig = plt.figure(figsize=(10, 10))  # Create a new figure with a custom size (you can adjust the size as needed)

for j in range(len(overall)):
    ax = fig.add_subplot(3, 2, j+1) 
    xlabel = "None"
    ylabel = "Time in s"
    if j == 4 or j == 5: # weak scaling
        xlabel = "Lines per Process"
        xvals = overall[j][:,2] / overall[j][:,0]
        yvals = overall[j][:,3]
    elif j == 0 or j == 1:
        xlabel = "Number of Nodes"
        xvals = overall[j][:,1]
        yvals = overall[j][:,3]
    else:
        xlabel = "Number of Processes"
        xvals = overall[j][:,0]
        yvals = overall[j][:,3]
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    data = list(zip(xvals, yvals))
    table = ax.table(cellText=data, loc='upper left', colLabels=[xlabel, ylabel])
    #set space not so big
    table.scale(1, 1.5)
    table.auto_set_font_size(True)
    ax.set_title(names[j])
    #lower margin
    plt.subplots_adjust(bottom=0.2)
    ax.axis('off')

plt.savefig("tables.png")  # Save the figure with all subplots
plt.show()  # Show the figure (remove this if you don't want to see it)