package main

import (
	"fmt"

	gmcl "github.com/alinush/go-mcl"
)

func main() {
	fmt.Println("Hello, World!")
	gmcl.InitFromString("bls12-381")
}
