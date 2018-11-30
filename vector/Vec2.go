package vector

import "math"

type Vec2 struct {
	X, Y float64
}

func (v *Vec2) Add(p Vec2) {
	v.X += p.X
	v.Y += p.Y
}

func (v *Vec2) Sub(p Vec2) {
	v.X -= p.X
	v.Y -= p.Y
}

func (v *Vec2) Mul(k float64) {
	v.X *= k
	v.Y *= k
}

func (v Vec2) Dot(p Vec2) float64 {
	return v.X*p.X + v.Y*p.Y
}

func (v Vec2) Cross(p Vec2) float64 {
	return v.X*p.Y - v.Y*p.X
}

func (v Vec2) Norm2() float64 {
	return v.X*v.X + v.Y*v.Y
}

func (v Vec2) Norm() float64 {
	return math.Sqrt(v.Norm2())
}
