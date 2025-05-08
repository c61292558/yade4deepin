import matplotlib.pyplot as plt

def read_data(filename):
    times = []
    velocities = []
    positions = []

    with open(filename, 'r') as file:
        # Skip header line if there is one
        header = file.readline()

        for line in file:
            # Split the line into values
            values = line.strip().split()
            if len(values) == 3:
                time, velocity, position = map(float, values)
                times.append(time)
                velocities.append(velocity)
                positions.append(position)

    return times, velocities, positions


################ Analytical solution ############
densityF=1000
nu=0.0001
mu=1000*nu
densityParticle=1500
Radius=0.0005
g=-9.81

TermVel_analytical=(2*Radius*Radius*(densityParticle-densityF)*g)/(9*mu)
Reyonolds=abs(densityF*2*Radius*TermVel_analytical)/mu
print ("Reynolds number: ",Reyonolds)
print ("Terminal velocity: ",TermVel_analytical)


################ Plots ############

filename = 'data.txt'
times, velocities, positions = read_data(filename)

timeAnalytical=[times[0],times[-1]]
velAnalytical=[TermVel_analytical,TermVel_analytical]


fig, axs = plt.subplots(2, 1, figsize=(10, 8))

# Plot Time vs Velocity
axs[0].plot(times, velocities, marker='None', linestyle='--', color='b')
axs[0].plot(timeAnalytical, velAnalytical, marker='None', linestyle='-', color='k')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Velocity')

# Plot Time vs Position
axs[1].plot(times, positions, marker='None', linestyle='--', color='b')
axs[1].set_ylabel('Position')

# Adjust layout
plt.tight_layout()

# Show plots
plt.savefig("terminalVel.png")
plt.show()
