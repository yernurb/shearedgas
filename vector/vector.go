package vector

import "math"

// Vector is a 2D structure
type Vector struct {
	X, Y float64
}

// ***** Methods *****

// Magnitude - returns the magnitude of the vector
func (p *Vector) Magnitude() float64 {
	return math.Sqrt(p.X*p.X + p.Y*p.Y)
}

// Null - sets vector to zero
func (p *Vector) Null() {
	p.X = 0.0
	p.Y = 0.0
}

// Scale - scales the vector
func (p *Vector) Scale(s float64) {
	p.X *= s
	p.Y *= s
}

// ***** Applications *****

// Add - returns sum of two vectors
func Add(v1 Vector, v2 Vector) Vector {
	var res Vector
	res.X = v1.X + v2.X
	res.Y = v1.Y + v2.Y

	return res
}

// Sub - returns difference of two vectors
func Sub(v1 Vector, v2 Vector) Vector {
	var res Vector
	res.X = v1.X - v2.X
	res.Y = v1.Y - v2.Y

	return res
}

// Dot - returns the dot product of two vectors
func Dot(v1 Vector, v2 Vector) float64 {
	return v1.X*v2.X + v1.Y*v2.Y
}

// Cross - returns cross product of two vectors
func Cross(v1 Vector, v2 Vector) float64 {
	return v1.X*v2.Y - v1.Y*v2.X
}

// Mul - returns scaled copy of the vector
func Mul(v Vector, s float64) Vector {
	var res Vector
	res.X = s * v.X
	res.Y = s * v.Y

	return res
}
