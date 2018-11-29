package main

import (
	"fmt"
	"github.com/yernurb/shearedgas/vector"
)

func main() {
	v1 := vector.Vec2{1.2, 2.3}
	v2 := vector.Vec2{-4.3, 2.1}
	fmt.Println("Hello, World")
	fmt.Println(v1.X, v1.Y)
	fmt.Println(v2.X, v2.Y)
	v1.X = 4.0
	v1.Y = 3.0
	fmt.Println(v1.Norm2())
}
