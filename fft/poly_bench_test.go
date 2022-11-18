package fft

import (
	"fmt"
	"testing"

	"github.com/accumulators-agg/go-poly/ff"
	gmcl "github.com/alinush/go-mcl"
)

func BenchmarkPolyMul(b *testing.B) {

	for scale := uint8(4); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_ = PolyMul(A, B)
			}
		})
	}
}

func BenchmarkPolyTree(b *testing.B) {

	for scale := uint8(4); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_ = PolyTree(A)
			}
		})
	}
}

func BenchmarkPolyDiv(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_, _ = PolyDiv(A, B)
			}
		})
	}
}

func BenchmarkPolyDifferentiate(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		data := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			data[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_ = PolyDifferentiate(data)
			}
		})
	}
}

func BenchmarkPolyXGCD1Balanced(b *testing.B) {

	for scale := uint8(10); scale < 13; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_, _, _ = xGCD1(A, B)
			}
		})
	}
}

func BenchmarkPolyXGCD2Balanced(b *testing.B) {

	for scale := uint8(10); scale < 13; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
			B[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_, _, _ = xGCD2(A, B)
			}
		})
	}
}
func BenchmarkPolyXGCD1Lopsided(b *testing.B) {

	for scale := uint8(10); scale < 13; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_, _, _ = xGCD1(A, B)
			}
		})
	}
}

func BenchmarkPolyXGCD2Lopsided(b *testing.B) {

	for scale := uint8(10); scale < 13; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		B := make([]gmcl.Fr, 2, 2)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		B[0] = *(ff.RandomFr())
		B[1] = *(ff.RandomFr())
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_, _, _ = xGCD2(A, B)
			}
		})
	}
}

func BenchmarkPolySubProductTree(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_ = SubProductTree(A)
			}
		})
	}
}

func BenchmarkPolyMultiEvaluate(b *testing.B) {

	for scale := uint8(10); scale < 15; scale++ {
		n := uint64(1) << scale
		A := make([]gmcl.Fr, n, n)
		for i := uint64(0); i < n; i++ {
			A[i] = *(ff.RandomFr())
		}
		b.Run(fmt.Sprintf("scale_%d", scale), func(t *testing.B) {
			M := SubProductTree(A)
			t.ResetTimer()
			for i := 0; i < t.N; i++ {
				_ = PolyMultiEvaluate(A, M)
			}
		})
	}
}
