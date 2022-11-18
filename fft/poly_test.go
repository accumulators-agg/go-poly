package fft

import (
	"fmt"
	"testing"

	"github.com/accumulators-agg/go-poly/debug"
	"github.com/accumulators-agg/go-poly/ff"
	gmcl "github.com/alinush/go-mcl"
)

func CheckEqualVec(a []gmcl.Fr, b []gmcl.Fr) bool {
	n := len(a)
	if n == len(b) && n > 0 {
		flag := true
		for i := 0; i < n; i++ {
			flag = flag && a[i].IsEqual(&b[i])
		}
		return flag
	}
	return false
}

func printSubTree(M [][][]gmcl.Fr, printMsg string) {
	for i := 0; i < len(M); i++ {
		for j := 0; j < len(M[i]); j++ {
			msg := fmt.Sprintf("%s [%d, %d]", printMsg, i, j)
			debug.DebugFrs(msg, M[i][j])
		}
	}
}

func TestPolyIsPolyZero(t *testing.T) {
	var tests = []struct {
		a    []int64
		want bool
	}{
		{[]int64{1, 2, 3, 4, 7, 8}, false},
		{[]int64{0}, true},
		{[]int64{1, 2, 3, 4, 7, 8, 0, 0}, false},
		{[]int64{0, 0, 0, 0, 0}, true},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			ans := IsPolyZero(aFr)

			if ans != tt.want {
				t.Errorf("IsPolyZero: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyCondense(t *testing.T) {
	var tests = []struct {
		a, want []int64
	}{
		{[]int64{1, 2, 3, 4, 7, 8}, []int64{1, 2, 3, 4, 7, 8}},
		{[]int64{0}, []int64{0}},
		{[]int64{0, 0, 0, 0, 0}, []int64{0}},
		{[]int64{1, 2, 3, 4, 7, 8, 0, 0, 0, 0, 0}, []int64{1, 2, 3, 4, 7, 8}},
		{[]int64{1, 2, 3, 4, 7, -8, 0, 0, 0, 0, 0}, []int64{1, 2, 3, 4, 7, -8}},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolyCondense(aFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolyCondense: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyAdd(t *testing.T) {
	var tests = []struct {
		a, b, want []int64
	}{
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{5, 1, 3},
			[]int64{6, 3, 6, 4, 7, 8},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8}, []int64{0}, []int64{1, 2, 3, 4, 7, 8},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{-1, -2, -3, -4, -7, -8},
			[]int64{0},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{1, 2, 3, 4, 7, -8},
			[]int64{2, 4, 6, 8, 14},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{0, 0, 0, 0, 0, -8},
			[]int64{1, 2, 3, 4, 7},
		},
		{
			[]int64{0, 0, 0, 0},
			[]int64{1, 2, 3, 4, 7, -8},
			[]int64{1, 2, 3, 4, 7, -8},
		},
		{
			[]int64{1, 2, 3, 4, 7, -8},
			[]int64{0, 0, 0, 0},
			[]int64{1, 2, 3, 4, 7, -8},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolyAdd(aFr, bFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs("", ansFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolyAdd: Answer did not match with expected.")
			}
		})
	}
}

func TestPolySub(t *testing.T) {
	var tests = []struct {
		a, b, want []int64
	}{
		{
			[]int64{0, 0, 0, 0},
			[]int64{1, 2, 3, 4, 7, -8},
			[]int64{-1, -2, -3, -4, -7, 8},
		},
		{
			[]int64{1, 2, 3, 4, 7, -8},
			[]int64{0, 0, 0, 0},
			[]int64{1, 2, 3, 4, 7, -8},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{5, 1, 3},
			[]int64{-4, 1, 0, 4, 7, 8},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8}, []int64{0}, []int64{1, 2, 3, 4, 7, 8},
		},
		{
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{1, 2, 3, 4, 7, 8},
			[]int64{0},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolySub(aFr, bFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs("", ansFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolySub: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyMul(t *testing.T) {
	var tests = []struct {
		a, b, want []int64
	}{
		{[]int64{1, 2, 3, 4, 7, 8}, []int64{5, 1, 3}, []int64{5, 11, 20, 29, 48, 59, 29, 24}},
		{[]int64{5, 1, 3}, []int64{1, 2, 3, 4, 7}, []int64{5, 11, 20, 29, 48, 19, 21}},
		{[]int64{1, 2, 3, 4, 7, 8}, []int64{0}, []int64{0}},
		{[]int64{0}, []int64{1, 2, 3, 4, 7, 8}, []int64{0}},
		{[]int64{1}, []int64{0}, []int64{0}},
		{[]int64{112}, []int64{2}, []int64{224}},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolyMul(aFr, bFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs("", ansFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolyMul: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyLongDiv(t *testing.T) {

	var tests = []struct {
		a, b, want []int64
	}{
		{[]int64{1, 2, 3, 4}, []int64{5, 1}, []int64{87, -17, 4}},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolyLongDiv(aFr, bFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolyLongDiv: Quotient did not match with expected")
			}
		})
	}
}

func TestPolyDiv(t *testing.T) {

	var tests = []struct {
		a, b, qwant, rwant []int64
	}{
		{[]int64{1, 2, 3, 4}, []int64{5, 1}, []int64{87, -17, 4}, []int64{-434}},
		{[]int64{8, 10, -5, 3}, []int64{-3, 2, 1}, []int64{-11, 3}, []int64{-25, 41}},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			qwantFr := ff.FromInt64Vec(tt.qwant)
			rwantFr := ff.FromInt64Vec(tt.rwant)
			qFr, rFr := PolyDiv(aFr, bFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs("", ansFr)

			var flag bool

			flag = CheckEqualVec(qwantFr, qFr)
			if flag == false {
				t.Errorf("PolyDiv: Quotient did not match with expected.")
			}
			flag = CheckEqualVec(rwantFr, rFr)
			if flag == false {
				t.Errorf("PolyDiv: Remainder did not match with expected.")
			}
		})
	}
}

func TestPolyDerivative(t *testing.T) {
	var tests = []struct {
		a, want []int64
	}{
		{[]int64{1, 2, 3, 4, 7, 8}, []int64{2, 6, 12, 28, 40}},
		{[]int64{2, 6, 12, 28, 40}, []int64{6, 24, 84, 160}},
		{[]int64{0}, []int64{0}},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			wantFr := ff.FromInt64Vec(tt.want)
			ansFr := PolyDifferentiate(aFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs("", ansFr)

			flag := CheckEqualVec(wantFr, ansFr)
			if flag == false {
				t.Errorf("PolyDifferentiate: Derivate did not match with expected.")
			}
		})
	}
}

func TestPolyExtGCD(t *testing.T) {

	var tests = []struct {
		a []int64
		b []int64
		g []int64
		u []int64
		v []int64
	}{
		{
			[]int64{1, 1, 1, 1},
			[]int64{1, 0, 0, 1},
			[]int64{1, 1},
			[]int64{1, -1},
			[]int64{0, 1},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			bFr := ff.FromInt64Vec(tt.b)
			gFr := ff.FromInt64Vec(tt.g)
			uFr := ff.FromInt64Vec(tt.u)
			vFr := ff.FromInt64Vec(tt.v)

			var g, u, v []gmcl.Fr
			var flag bool
			flag = true

			g, u, v = xGCD1(aFr, bFr)
			flag = flag && CheckEqualVec(gFr, g)
			flag = flag && CheckEqualVec(uFr, u)
			flag = flag && CheckEqualVec(vFr, v)

			if flag == false {
				t.Errorf("xGCD1: Answer did not match with expected.")
			}
			g, u, v = xGCD2(aFr, bFr)
			flag = flag && CheckEqualVec(gFr, g)
			flag = flag && CheckEqualVec(uFr, u)
			flag = flag && CheckEqualVec(vFr, v)

			if flag == false {
				t.Errorf("xGCD2: Answer did not match with expected.")
			}
			// debug.DebugFrs("g", g)
			// debug.DebugFrs("u", u)
			// debug.DebugFrs("v", v)
			// debug.DebugFrs("gFr", gFr)
			// debug.DebugFrs("uFr", uFr)
			// debug.DebugFrs("vFr", vFr)
		})
	}
}

func TestPolySubProdTree(t *testing.T) {
	var tests = []struct {
		a    []int64
		want [][][]int64
	}{
		{
			[]int64{1, 2, 3, 4},
			[][][]int64{
				{
					{-1, 1},
					{-2, 1},
					{-3, 1},
					{-4, 1},
				},
				{
					{2, -3, 1},
					{12, -7, 1},
				},
				{
					{24, -50, 35, -10, 1},
				},
			},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			var wantFr [][][]gmcl.Fr

			wantFr = make([][][]gmcl.Fr, len(tt.want))
			for i := 0; i < len(wantFr); i++ {
				wantFr[i] = make([][]gmcl.Fr, len(tt.want[i]))
				for j := 0; j < len(tt.want[i]); j++ {
					wantFr[i][j] = ff.FromInt64Vec(tt.want[i][j])
				}
			}

			ansFr := SubProductTree(aFr)

			// debug.DebugFrs("", wantFr)
			// debug.DebugFrs(fmt.Sprintf("%d", len(ansFr)), wantFr)

			// printSubTree(wantFr, "Want")

			flag := true
			for i := 0; i < len(ansFr); i++ {
				for j := 0; j < len(ansFr[i]); j++ {
					flag = flag && CheckEqualVec(wantFr[i][j], ansFr[i][j])
					// msg := fmt.Sprintf("Ans [%d, %d]", i, j)
					// debug.DebugFrs(msg, ansFr[i][j])
				}
			}
			if flag == false {
				t.Errorf("SubProductTree: Answer did not match with expected.")
			}
		})
	}

	LIMIT := 10
	for l := 2; l < LIMIT; l++ {
		testname := fmt.Sprintf("scale-%d", l)
		t.Run(testname, func(t *testing.T) {
			n := 1 << l
			aFr := make([]gmcl.Fr, n, n)
			for i := 0; i < n; i++ {
				aFr[i] = *ff.RandomFr()
			}
			var temp gmcl.Fr
			result := []gmcl.Fr{ff.ONE}
			tempPoly := make([]gmcl.Fr, 2)

			for i := 0; i < n; i++ {
				gmcl.FrNeg(&temp, &aFr[i])
				tempPoly[0] = temp
				tempPoly[1] = ff.ONE
				result = PolyMul(result, tempPoly)
			}
			M := SubProductTree(aFr)
			// debug.DebugFrs("SubTree", M[len(M)-1][0])
			// debug.DebugFrs("Result", result)
			var flag bool
			flag = CheckEqualVec(M[len(M)-1][0], result)
			if flag == false {
				t.Errorf("SubProdTree 2: Answer did not match with expected.")
			}
			m := PolyTree(aFr)
			flag = CheckEqualVec(result, m)
			if flag == false {
				t.Errorf("PolyTree: Answer did not match with expected.")
			}
		})
	}
}

func TestPolySplitSubProdTree(t *testing.T) {
	var tests = []struct {
		a     []int64
		wantL [][][]int64
		wantR [][][]int64
	}{
		{
			[]int64{1, 2, 3, 4},
			[][][]int64{
				{
					{-1, 1},
					{-2, 1},
				},
				{
					{2, -3, 1},
				},
			},
			[][][]int64{
				{
					{-3, 1},
					{-4, 1},
				},
				{
					{12, -7, 1},
				},
			},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			aFr := ff.FromInt64Vec(tt.a)
			var wantLFr, wantRFr [][][]gmcl.Fr

			wantLFr = make([][][]gmcl.Fr, len(tt.wantL))
			wantRFr = make([][][]gmcl.Fr, len(tt.wantR))
			for i := 0; i < len(wantLFr); i++ {
				wantLFr[i] = make([][]gmcl.Fr, len(tt.wantL[i]))
				wantRFr[i] = make([][]gmcl.Fr, len(tt.wantR[i]))
				for j := 0; j < len(tt.wantL[i]); j++ {
					wantLFr[i][j] = ff.FromInt64Vec(tt.wantL[i][j])
					wantRFr[i][j] = ff.FromInt64Vec(tt.wantR[i][j])
				}
			}

			ansFr := SubProductTree(aFr)
			L, R := splitSubProdTree(ansFr)

			// printSubTree(L, "L")
			// printSubTree(R, "R")

			flag := true
			for i := 0; i < len(L); i++ {
				for j := 0; j < len(L[i]); j++ {
					flag = flag && CheckEqualVec(wantLFr[i][j], L[i][j])
					flag = flag && CheckEqualVec(wantRFr[i][j], R[i][j])
					// msg := fmt.Sprintf("Ans [%d, %d]", i, j)
					// debug.DebugFrs(msg, ansFr[i][j])
				}
			}

			if flag == false {
				t.Errorf("splitSubProdTree: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyMultiPointEval(t *testing.T) {

	var tests = []struct {
		polynomial  []int64
		evalPoints  []int64
		evaluations []int64
	}{
		{
			[]int64{1, 2, 3, 4},
			[]int64{1, 2, 3, 4},
			[]int64{10, 49, 142, 313},
		},
		{
			[]int64{12, 3, -17, 78, -91, 49, -82, 61},
			[]int64{10, 12, 11, 23, 8, 28, 13, 1},
			[]int64{532066342, 1951327776, 1050110403, 195846264759, 107702244, 784342707664, 3447624049, 13},
		},
	}

	for counter, tt := range tests {
		testname := fmt.Sprintf("%d", counter+1)
		t.Run(testname, func(t *testing.T) {

			polynomialFr := ff.FromInt64Vec(tt.polynomial)
			evalPointsFr := ff.FromInt64Vec(tt.evalPoints)
			evaluationsFr := ff.FromInt64Vec(tt.evaluations)

			subproductTree := SubProductTree(evalPointsFr)
			ansFr := PolyMultiEvaluate(polynomialFr, subproductTree)

			flag := CheckEqualVec(ansFr, evaluationsFr)

			if flag == false {
				t.Errorf("PolyMultiEvaluate: Answer did not match with expected.")
			}
		})
	}

	LIMIT := 10
	for l := 3; l < LIMIT; l++ {
		testname := fmt.Sprintf("MultiPointEval/scale-%d", l)
		t.Run(testname, func(t *testing.T) {
			n := 1 << l

			aFr := make([]gmcl.Fr, n-1, n-1) // Not n!
			evalPointsFr := make([]gmcl.Fr, n, n)
			evaluations := make([]gmcl.Fr, n, n)

			for i := 0; i < n; i++ {
				if i < n-1 {
					aFr[i] = *ff.RandomFr()
				}
				evalPointsFr[i] = *ff.RandomFr()
			}

			polynomialFr := PolyTree(aFr)
			M := SubProductTree(evalPointsFr)

			for i := 0; i < n; i++ {
				gmcl.FrEvaluatePolynomial(&evaluations[i], polynomialFr, &evalPointsFr[i])
			}
			ansFr := PolyMultiEvaluate(polynomialFr, M)

			flag := CheckEqualVec(ansFr, evaluations)
			if flag == false {
				t.Errorf("PolyMultiPointEval 2: Answer did not match with expected.")
			}
		})
	}
}

func TestPolyTreeVec(t *testing.T) {

	aLen := uint64(25)
	aFr := make([]gmcl.Fr, aLen)
	for i := uint64(0); i < aLen; i++ {
		aFr[i].Random()
	}
	N := nextPowOf2(uint64(aLen))
	M := make([][]gmcl.Fr, N)
	for j := uint64(0); j < N; j++ {
		if j < aLen {
			M[j] = make([]gmcl.Fr, 2)
			gmcl.FrNeg(&M[j][0], &aFr[j])
			M[j][1].SetInt64(1)
		} else {
			M[j] = make([]gmcl.Fr, 1)
			M[j][0].SetInt64(1)
		}
	}

	m := PolyTree(aFr)
	result := PolyTreeVec(M)
	flag := CheckEqualVec(result, m)
	if flag == false {
		t.Errorf("PolyTreeVec: Answer did not match with expected.")
	}
}
