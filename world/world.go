package world

import "github.com/yernurb/shearedgas/particle"

type World struct {
	Particles      []particle.Particle // Array of particles in the world
	X_size, Y_size float64             // X and Y sizes of the world
	Omega          float64             // Orbital speed
	cell_x, cell_y float64             // Sizes of each cell in the world (private)
	Cell           [][][]int           // Cell structure of the physical world. Contains particle indices
}
