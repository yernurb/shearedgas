import numpy as np
import matplotlib.pyplot as plt


def maxwell_2D(velocity, T):
    return velocity*np.exp(-velocity*velocity/(2*T)) / T

def maxwell_3D(velocity, T):
    return 4*np.pi*velocity*velocity*np.exp(-velocity*velocity/(2*T)) / np.exp(1.5*np.log(2*np.pi*T))

class DSMC3D:
    def __init__(self):
        self.N = 10000  # number of particles
        self.Nbox = 10  # number of boxes per dimension
        self.R = 1       # radii of particles   
        self.density = 0.25   # number density
        self.L = np.exp(np.log(self.N / self.density)/3)    # total size of the simulation area
        self.eps = 0.8   # coefficient of restitution
        self.tend = 100 # simulation time length
        self.T0 = 0 # initial temperature

        # theoretical parameters
        self.C1 = 4*np.pi*self.R*self.R*self.Nbox*self.Nbox*self.Nbox / (self.L*self.L*self.L)
        self.a2 = 16*(1-self.eps)*(1-2*self.eps*self.eps) / (81-17*self.eps+30*self.eps*self.eps*(1-self.eps))
        self.eta = 4*np.pi*self.R*self.R*self.R*self.density / 3
        self.g2 = (2-self.eta) / 2 / (1-self.eta) / (1-self.eta)
        self.ctime = 0
        self.th = 0


        self.lattice = []    # simulation boxes that contain particle indices
        self.error = []      # roundoff errors

        self.vmax = 0    # maximal velocity to cut the DF
        self.energy = 0  # total energy of the system
        self.dt = 0      # timestep of the simulation

        self.x  = np.zeros(self.N)    # x coordinate of particles      
        self.y  = np.zeros(self.N)    # y coordinate of particles 
        self.z  = np.zeros(self.N)    # z coordinate of particles
        self.vx = np.zeros(self.N)    # x velocity component of particles
        self.vy = np.zeros(self.N)    # y velocity component of particles
        self.vz = np.zeros(self.N)    # z velocity component of particles
        self.init()


    def vmax_estimate(self):
        return 20 * np.sqrt(self.energy / self.N)

    def init(self):
        print("Initializing system...")
        print("...")
        print("Number of simulated particles N =", self.N)
        print("Linear size of the system L =", self.L)
        print("Particle radii R =", self.R)
        print("Coefficient of restitution eps =", self.eps)
        print("Number of boxes in single dimension Nbox =", self.Nbox)
        for i in range(self.Nbox*self.Nbox*self.Nbox):
            self.lattice.append([])
            self.error.append(0)
        print("Total lattice cells:", i+1)
        print("Number density:", self.N/self.L/self.L/self.L)
        print("...")
        self.x = self.L * np.random.rand(self.N)
        self.y = self.L * np.random.rand(self.N)
        self.z = self.L * np.random.rand(self.N)
        self.vx = 3*np.random.randn(self.N)
        self.vy = 3*np.random.randn(self.N)
        self.vz = 3*np.random.randn(self.N)
        self.vx = self.vx - self.vx.mean()
        self.vy = self.vy - self.vy.mean()
        self.vz = self.vz - self.vz.mean()
        print("Coordinates and velocities are distributed")
        print("...")
        self.energy = 0.5 * (self.vx*self.vx + self.vy*self.vy + self.vz*self.vz).sum()
        self.T0 = 2 * self.energy / (3 * self.N) 
        print("Total initial energy of the system E =", self.energy)
        self.ctime = 16*np.pi*self.g2*self.density*self.R*self.R*np.sqrt(self.T0)
        self.th = (1-self.eps*self.eps)*(1+3*self.a2/16)*self.ctime / 8
        print("...")
        print("Theoretical parameters:")
        print("T0:", self.T0)
        print("eta:", self.eta)
        print("g2:", self.g2)
        print("a2:", self.a2)
        print("tc-1:", self.ctime)
        print("th-1", self.th)
        print("...")



    def propagate(self):
        self.x = (self.x + self.vx*self.dt + self.L) % self.L
        self.y = (self.y + self.vy*self.dt + self.L) % self.L
        self.z = (self.z + self.vz*self.dt + self.L) % self.L

    def sort(self):
        self.lattice = []
        for i in range(self.Nbox*self.Nbox*self.Nbox):
            self.lattice.append([])
        ix = (self.Nbox * self.x / self.L).astype(int)
        iy = (self.Nbox * self.y / self.L).astype(int)
        iz = (self.Nbox * self.z / self.L).astype(int)
        for i in range(self.N):
            index = ix[i] * self.Nbox*self.Nbox + iy[i] * self.Nbox + iz[i]
            self.lattice[index].append(i)

    def collision(self):
        ncols = 0
        for index in range(self.Nbox*self.Nbox*self.Nbox):
            npart = len(self.lattice[index])
            if npart > 1:
                dcoll = self.C1 * npart * (npart - 1) * self.vmax * self.dt + self.error[index]
                ncoll = int(dcoll)
                self.error[index] = dcoll - ncoll
                for icol in range(ncoll):
                    pindex1 = int(npart * np.random.rand())
                    pindex2 = (int(pindex1 + (npart-1)*np.random.rand()) + 1) % npart
                    p1 = self.lattice[index][pindex1]
                    p2 = self.lattice[index][pindex2]

                    ncols += self.paircollision(p1, p2)
        
        return ncols

    def paircollision(self, i, j):
        dvx = self.vx[i] - self.vx[j]
        dvy = self.vy[i] - self.vy[j]
        dvz = self.vz[i] - self.vz[j]
        phi = 2 * np.pi * np.random.rand()
        costheta = 2 * np.random.rand() - 1
        sintheta = np.sqrt(1 - costheta*costheta)
        ndx = np.cos(phi)*sintheta
        ndy = np.sin(phi)*sintheta
        ndz = costheta
        vnorm = dvx * ndx + dvy * ndy + dvz * ndz

        if np.fabs(vnorm) < np.random.rand() * self.vmax:
            return 0
        else:
            h = (1 + self.eps) * vnorm / 2
            self.vx[i] -= h * ndx
            self.vy[i] -= h * ndy
            self.vz[i] -= h * ndz
            self.vx[j] += h * ndx
            self.vy[j] += h * ndy
            self.vz[j] += h * ndz
            self.energy = self.energy - (1 - self.eps * self.eps) * vnorm * vnorm / 4
            return 1

    def make_step(self):    
        #self.energy = 0.5 * (self.vx*self.vx + self.vy*self.vy + self.vz*self.vz).sum()
        self.vmax = self.vmax_estimate()
        self.dt = 0.2 * self.L / self.Nbox / self.vmax
        self.propagate()
        self.sort()
        ncols = self.collision()
        return (ncols, self.dt)

    def get_distribution(self):
        vel = np.sqrt(self.vx*self.vx + self.vy*self.vy + self.vz*self.vz)
        df, bins = np.histogram(vel, bins=100, density=True)
        bins = bins[:-1] + np.diff(bins) / 2
        return bins, df


class DSMC2D:
    def __init__(self):
        self.N = 10000  # number of particles
        self.Nbox = 10  # number of boxes per dimension
        self.R = 1       # radii of particles   
        self.density = 0.25 # number density of the system
        self.L = np.sqrt(self.N / self.density)    # total size of the simulation area
        self.eps = 0.8   # coefficient of restitution
        self.tend = 100  # simulation time length
        self.T0 = 0 # initial temperature

        # theoretical parameters
        self.C1 = 2*np.pi*self.R*self.Nbox*self.Nbox / (self.L*self.L)
        self.a2 = 16*(1-self.eps)*(1-2*self.eps*self.eps) / (57-25*self.eps+30*self.eps*self.eps*(1-self.eps))
        self.eta = np.pi*self.R*self.R*self.density
        self.g2 = (1-7*self.eta/16) / (1-self.eta) / (1-self.eta)
        self.ctime = 0
        self.th = 0

        self.lattice = []    # simulation boxes that contain particle indices
        self.error = []      # roundoff errors

        self.vmax = 0    # maximal velocity to cut the DF
        self.energy = 0  # total energy of the system
        self.dt = 0      # timestep of the simulation

        self.x  = np.zeros(self.N)    # x coordinate of particles      
        self.y  = np.zeros(self.N)    # y coordinate of particles 
        self.vx = np.zeros(self.N)    # x velocity component of particles
        self.vy = np.zeros(self.N)    # y velocity component of particles
        self.init()


    def vmax_estimate(self):
        return 20 * np.sqrt(self.energy / self.N)

    def init(self):
        print("Initializing system...")
        print("...")
        print("Number of simulated particles N =", self.N)
        print("Linear size of the system L =", self.L)
        print("Particle radii R =", self.R)
        print("Coefficient of restitution eps =", self.eps)
        print("Number of boxes in single dimension Nbox =", self.Nbox)
        for i in range(self.Nbox*self.Nbox*self.Nbox):
            self.lattice.append([])
            self.error.append(0)
        print("Total lattice cells:", i+1)
        print("Number density:", self.N/self.L/self.L)
        print("...")
        self.x = self.L * np.random.rand(self.N)
        self.y = self.L * np.random.rand(self.N)
        self.vx = 10*np.random.randn(self.N)
        self.vy = 10*np.random.randn(self.N)
        self.vx = self.vx - self.vx.mean()
        self.vy = self.vy - self.vy.mean()
        print("Coordinates and velocities are distributed")
        print("...")
        self.energy = 0.5 * (self.vx*self.vx + self.vy*self.vy).sum()
        self.T0 = self.energy / self.N 
        print("Total initial energy of the system E =", self.energy)
        self.ctime = 4*self.g2*self.density*self.R*np.sqrt(np.pi*self.T0)
        self.th = (1-self.eps*self.eps)*(1+3*self.a2/16)*self.ctime / 4
        print("...")
        print("Theoretical parameters:")
        print("T0:", self.T0)
        print("eta:", self.eta)
        print("g2:", self.g2)
        print("a2:", self.a2)
        print("tc-1:", self.ctime)
        print("th-1", self.th)
        print("...")


    def propagate(self):
        self.x = (self.x + self.vx*self.dt + self.L) % self.L
        self.y = (self.y + self.vy*self.dt + self.L) % self.L

    def sort(self):
        self.lattice = []
        for i in range(self.Nbox*self.Nbox):
            self.lattice.append([])
        ix = (self.Nbox * self.x / self.L).astype(int)
        iy = (self.Nbox * self.y / self.L).astype(int)
        for i in range(self.N):
            index = ix[i] * self.Nbox + iy[i]
            self.lattice[index].append(i)

    def collision(self):
        ncols = 0
        for index in range(self.Nbox*self.Nbox):
            npart = len(self.lattice[index])
            if npart > 1:
                dcoll = self.C1 * npart * (npart - 1) * self.vmax * self.dt + self.error[index]
                ncoll = int(dcoll)
                self.error[index] = dcoll - ncoll
                for icol in range(ncoll):
                    pindex1 = int(npart * np.random.rand())
                    pindex2 = (int(pindex1 + (npart-1)*np.random.rand()) + 1) % npart
                    p1 = self.lattice[index][pindex1]
                    p2 = self.lattice[index][pindex2]

                    ncols += self.paircollision(p1, p2)
        
        return ncols

    def paircollision(self, i, j):
        dvx = self.vx[i] - self.vx[j]
        dvy = self.vy[i] - self.vy[j]
        phi = 2 * np.pi * np.random.rand()
        ndx = np.cos(phi)
        ndy = np.sin(phi)
        vnorm = dvx * ndx + dvy * ndy

        if np.fabs(vnorm) < np.random.rand() * self.vmax:
            return 0
        else:
            h = (1 + self.eps) * vnorm / 2
            self.vx[i] -= h * ndx
            self.vy[i] -= h * ndy
            self.vx[j] += h * ndx
            self.vy[j] += h * ndy
            self.energy = self.energy - (1 - self.eps * self.eps) * vnorm * vnorm / 4
            return 1

    def make_step(self):    
        #self.energy = 0.5 * (self.vx*self.vx + self.vy*self.vy).sum()
        self.vmax = self.vmax_estimate()
        self.dt = 0.2 * self.L / self.Nbox / self.vmax
        self.propagate()
        self.sort()
        ncols = self.collision()
        return ncols, self.dt
    
    def get_distribution(self):
        vel = np.sqrt(self.vx*self.vx + self.vy*self.vy)
        df, bins = np.histogram(vel, bins=100, density=True)
        bins = bins[:-1] + np.diff(bins) / 2
        return bins, df




homo = DSMC2D()
T0 = homo.T0
time = 0
it = 0
ncols = 0
energy = []
timel = []

while time < homo.tend:
    energy.append(homo.energy)
    timel.append(time)
    scols, dt = homo.make_step()
    ncols += scols
    time += dt
    it += 1
    print("Time:", time)
    print("Time step:", dt)
    print("Number of collisions:", scols)
    print("Total collisions:", ncols)
    print("Energy:", homo.energy)
    print("Vmax:", homo.vmax)
    print("....")

'''
time = np.asarray(timel)
energy = 2*np.asarray(energy) / (2*homo.N)

#homo.th = 8 * (1-homo.eps*homo.eps)*homo.R*homo.R*homo.density*np.sqrt(np.pi*homo.T0) / 3
#homo.th = 4 * np.pi * (1-homo.eps*homo.eps)*homo.R*homo.density*np.sqrt(homo.T0)


energy_th = homo.T0 / (1 + time*homo.th) / (1 + time*homo.th)

print(time[0])

plt.plot(np.log(time), np.log(energy))
plt.plot(np.log(time), np.log(energy_th))

plt.show()
'''

vel, df = homo.get_distribution()
plt.plot(vel, df)

Temp = homo.energy / (homo.N)

maxwell = maxwell_2D(vel, Temp)
print((df*vel).sum())
print(maxwell.sum())
plt.plot(vel, maxwell)
print(it)

#plt.scatter(homo.vx, homo.vy)
plt.show()


