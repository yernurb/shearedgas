package particle

import "github.com/yernurb/shearedgas/vector"

// Particle structure contains all necessary information about particles
type Particle struct {
	Pos                               vector.Vector // Position vector of the particle
	Drdt0, Drdt1, Drdt2, Drdt3, Drdt4 vector.Vector // Time derivatives of the position vector
	M, R, E, Y, N, G, A               float64       // mass, radius, restitution, Young modulus, Poisson ratio, adhesion, damping parameters
	Force                             vector.Vector // Resulting force acting on the particle
}

// AddForce - adds force to the particle
func (p *Particle) AddForce(f vector.Vector) {
	vector.Add(p.Force, f)
}

// ResetForce - set acting force to zero
func (p *Particle) ResetForce() {
	p.Force.Null()
}

// Predict - predictor of the Gear's integrating scheme
func (p *Particle) Predict(dt float64) {
	a1 := dt
	a2 := a1 * dt / 2
	a3 := a2 * dt / 3
	a4 := a3 * dt / 4

	rdt0 := vector.Add(vector.Add(vector.Add(vector.Mul(p.Drdt1, a1), vector.Mul(p.Drdt2, a2)), vector.Mul(p.Drdt3, a3)), vector.Mul(p.Drdt4, a4))
	rdt1 := vector.Add(vector.Add(vector.Mul(p.Drdt2, a1), vector.Mul(p.Drdt3, a2)), vector.Mul(p.Drdt4, a3))
	rdt2 := vector.Add(vector.Mul(p.Drdt3, a1), vector.Mul(p.Drdt4, a2))
	rdt3 := vector.Mul(p.Drdt4, a1)

	p.Drdt0 = vector.Add(p.Drdt0, rdt0)
	p.Drdt1 = vector.Add(p.Drdt1, rdt1)
	p.Drdt2 = vector.Add(p.Drdt2, rdt2)
	p.Drdt3 = vector.Add(p.Drdt3, rdt3)
}

// Correct - corrector of the Gear's integrating scheme
func (p *Particle) Correct(dt float64) {
	dtrez := float64(1) / dt

	coeff0 := (float64(19) / float64(90)) * (dt * dt / float64(2))
	coeff1 := (float64(3) / float64(4)) * (dt / float64(2))
	coeff3 := (float64(1) / float64(2)) * (float64(3) * dtrez)
	coeff4 := (float64(1) / float64(12)) * (float64(12) * (dtrez * dtrez))

	accel := vector.Mul(p.Force, float64(1)/p.M)
	corr := vector.Sub(accel, p.Drdt2)

	p.Drdt0 = vector.Add(p.Drdt0, vector.Mul(corr, coeff0))
	p.Drdt1 = vector.Add(p.Drdt1, vector.Mul(corr, coeff1))
	p.Drdt2 = accel
	p.Drdt3 = vector.Add(p.Drdt3, vector.Mul(corr, coeff3))
	p.Drdt4 = vector.Add(p.Drdt4, vector.Mul(corr, coeff4))
}
