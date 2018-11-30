package particle

import "github.com/yernurb/shearedgas/vector"

type Particle struct {
	Pos                              vector.Vec2 // Position vector of the particle
	Drdt, Drdt2, Drdt3, Drdt4, Drdt5 vector.Vec2 // Time derivatives of the position vector
	M, R, E, Y, N, G, A              float64     // mass, radius, restitution, Young modulus, Poisson ratio, adhesion, damping parameters
	Force                            vector.Vec2 // Resulting force acting on the particle
}

func (p *Particle) AddForce(f vector.Vec2) {
	p.Force.Add(f)
}

func (p *Particle) ResetForce() {
	p.Force.Null()
}
