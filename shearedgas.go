package main

import (
	"fmt"

	"github.com/yernurb/shearedgas/vector"
)

func main() {
	fmt.Println("Hello, World")

	a1, a2, a3 := 0.1, 0.2, 0.15

	v1 := vector.Vector{X: 0.3, Y: 0.2}
	v2 := vector.Vector{X: 0.15, Y: 0.1}
	v3 := vector.Vector{X: 1.2, Y: -0.4}

	res := vector.Add(vector.Add(vector.Mul(v1, a1), vector.Mul(v2, a2)), vector.Mul(v3, a3))

	fmt.Println(res)
	fmt.Println(v1, v2, v3)
	v1 = v2
	fmt.Println(v1, v2, v3)
}
