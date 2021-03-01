// Original: https://github.com/ethereum/research/blob/master/mimc_stark/recovery.py

package fft

import (
	"fmt"

	"github.com/sshravan/go-poly/ff"
)

// Calculates modular inverses [1/values[0], 1/values[1] ...]
func multiInv(values []ff.Fr) []ff.Fr {
	partials := make([]ff.Fr, len(values)+1, len(values)+1)
	partials[0] = values[0]
	for i := 0; i < len(values); i++ {
		ff.MulModFr(&partials[i+1], &partials[i], &values[i])
	}
	var inv ff.Fr
	var tmp ff.Fr
	ff.InvModFr(&inv, &partials[len(partials)-1])
	outputs := make([]ff.Fr, len(values), len(values))
	for i := len(values); i > 0; i-- {
		ff.MulModFr(&outputs[i-1], &partials[i-1], &inv)
		ff.CopyFr(&tmp, &inv)
		ff.MulModFr(&inv, &tmp, &values[i-1])
	}
	return outputs
}

// Generates q(x) = poly(k * x)
func pOfKX(poly []ff.Fr, k *ff.Fr) []ff.Fr {
	out := make([]ff.Fr, len(poly), len(poly))
	powerOfK := ff.ONE
	var tmp ff.Fr
	for i := range poly {
		ff.MulModFr(&out[i], &poly[i], &powerOfK)
		ff.CopyFr(&tmp, &powerOfK)
		ff.MulModFr(&powerOfK, &tmp, k)
	}
	return out
}

func inefficientOddEvenDiv2(positions []uint64) (even []uint64, odd []uint64) { // TODO optimize away
	for _, p := range positions {
		if p&1 == 0 {
			even = append(even, p>>1)
		} else {
			odd = append(odd, p>>1)
		}
	}
	return
}

// Return (x - root**positions[0]) * (x - root**positions[1]) * ...
// possibly with a constant factor offset
func (fs *FFTSettings) _zPoly(positions []uint64, rootsOfUnityStride uint64) []ff.Fr {
	// If there are not more than 4 positions, use the naive
	// O(n^2) algorithm as it is faster
	if len(positions) <= 4 {
		/*
		   root = [1]
		   for pos in positions:
		       x = roots_of_unity[pos]
		       root.insert(0, 0)
		       for j in range(len(root)-1):
		           root[j] -= root[j+1] * x
		   return [x % modulus for x in root]
		*/
		root := make([]ff.Fr, len(positions)+1, len(positions)+1)
		root[0] = ff.ONE
		i := 1
		var v ff.Fr
		var tmp ff.Fr
		for _, pos := range positions {
			x := &fs.ExpandedRootsOfUnity[pos*rootsOfUnityStride]
			root[i] = ff.ZERO
			for j := i; j >= 1; j-- {
				ff.MulModFr(&v, &root[j-1], x)
				ff.CopyFr(&tmp, &root[j])
				ff.SubModFr(&root[j], &tmp, &v)
			}
			i++
		}
		// We did the reverse representation of 'root' as the python code, to not insert at the start all the time.
		// Now turn it back around.
		for i, j := 0, len(root)-1; i < j; i, j = i+1, j-1 {
			root[i], root[j] = root[j], root[i]
		}
		return root
	}
	// Recursively find the zpoly for even indices and odd
	// indices, operating over a half-size subgroup in each case
	evenPositions, oddPositions := inefficientOddEvenDiv2(positions)
	left := fs._zPoly(evenPositions, rootsOfUnityStride<<1)
	right := fs._zPoly(oddPositions, rootsOfUnityStride<<1)
	invRoot := &fs.ReverseRootsOfUnity[rootsOfUnityStride]
	// Offset the result for the odd indices, and combine the two
	out := fs.mulPolysWithFFT(left, pOfKX(right, invRoot), rootsOfUnityStride)
	// Deal with the special case where mul_polys returns zero
	// when it should return x ^ (2 ** k) - 1
	isZero := true
	for i := range out {
		if !ff.EqualZero(&out[i]) {
			isZero = false
			break
		}
	}
	if isZero {
		// TODO: it's [1] + [0] * (len(out) - 1) + [modulus - 1] in python, but strange it's 1 larger than out
		out[0] = ff.ONE
		for i := 1; i < len(out); i++ {
			out[i] = ff.ZERO
		}
		last := ff.MODULUS_MINUS1
		out = append(out, last)
		return out
	} else {
		return out
	}
}

// TODO test unhappy case
const maxRecoverAttempts = 10

func (fs *FFTSettings) ErasureCodeRecover(vals []*ff.Fr) ([]ff.Fr, error) {
	// Generate the polynomial that is zero at the roots of unity
	// corresponding to the indices where vals[i] is None
	positions := make([]uint64, 0, len(vals))
	for i := uint64(0); i < uint64(len(vals)); i++ {
		if vals[i] == nil {
			positions = append(positions, i)
		}
	}
	// TODO: handle len(positions)==0 case
	z := fs._zPoly(positions, fs.MaxWidth/uint64(len(vals)))
	//debug.DebugFrs("z", z)
	zVals, err := fs.FFT(z, false)
	if err != nil {
		return nil, err
	}
	//debug.DebugFrs("zvals", zVals)

	// Pointwise-multiply (vals filling in zero at missing spots) * z
	// By construction, this equals vals * z
	pTimesZVals := make([]ff.Fr, len(vals), len(vals))
	for i := uint(0); i < uint(len(vals)); i++ {
		if vals[i] == nil {
			// 0 * zVals[i] == 0
			pTimesZVals[i] = ff.ZERO
		} else {
			ff.MulModFr(&pTimesZVals[i], vals[i], &zVals[i])
		}
	}
	//debug.DebugFrs("p_times_z_vals", pTimesZVals)
	pTimesZ, err := fs.FFT(pTimesZVals, true)
	if err != nil {
		return nil, err
	}
	//debug.DebugFrs("p_times_z", pTimesZ)

	// Keep choosing k values until the algorithm does not fail
	// Check only with primitive roots of unity
	attempts := 0
	var kFr ff.Fr
	var tmp ff.Fr
	for k := uint64(2); attempts < maxRecoverAttempts; k++ {
		ff.AsFr(&kFr, k)
		// // TODO: implement this, translation of 'if pow(k, (modulus - 1) // 2, modulus) == 1:'
		//someOp(&tmp, &kFr)
		//if EqualOne(&tmp) {
		//	continue
		//}
		var invk ff.Fr
		ff.InvModFr(&invk, &kFr)
		// Convert p_times_z(x) and z(x) into new polynomials
		// q1(x) = p_times_z(k*x) and q2(x) = z(k*x)
		// These are likely to not be 0 at any of the evaluation points.
		pTimesZOfKX := pOfKX(pTimesZ, &kFr)
		//debug.DebugFrs("p_times_z_of_kx", pTimesZOfKX)
		pTimesZOfKXVals, err := fs.FFT(pTimesZOfKX, false)
		if err != nil {
			return nil, err
		}
		//debug.DebugFrs("p_times_z_of_kx_vals", pTimesZOfKXVals)
		zOfKX := pOfKX(z, &kFr)
		//debug.DebugFrs("z_of_kx", zOfKX)
		zOfKXVals, err := fs.FFT(zOfKX, false)
		if err != nil {
			return nil, err
		}
		//debug.DebugFrs("z_of_kx_vals", zOfKXVals)

		// Compute q1(x) / q2(x) = p(k*x)
		invZOfKXVals := multiInv(zOfKXVals)
		//debug.DebugFrs("inv_z_of_kv_vals", invZOfKXVals)
		pOfKxVals := make([]ff.Fr, len(pTimesZOfKXVals), len(pTimesZOfKXVals))
		for i := 0; i < len(pOfKxVals); i++ {
			ff.MulModFr(&pOfKxVals[i], &pTimesZOfKXVals[i], &invZOfKXVals[i])
		}
		//debug.DebugFrs("p_of_kx_vals", pOfKxVals)
		pOfKx, err := fs.FFT(pOfKxVals, true)
		if err != nil {
			return nil, err
		}
		//debug.DebugFrs("p_of_kx", pOfKx)

		// Given q3(x) = p(k*x), recover p(x)
		pOfX := make([]ff.Fr, len(pOfKx), len(pOfKx))
		if len(pOfKx) >= 1 {
			pOfX[0] = pOfKx[0]
		}
		if len(pOfKx) >= 2 {
			ff.MulModFr(&pOfX[1], &pOfKx[1], &invk)
			invKPowI := invk
			for i := 2; i < len(pOfKx); i++ {
				ff.CopyFr(&tmp, &invKPowI)
				ff.MulModFr(&invKPowI, &tmp, &invk)
				ff.MulModFr(&pOfX[i], &pOfKx[i], &invKPowI)
			}
		}
		output, err := fs.FFT(pOfX, false)
		if err != nil {
			return nil, err
		}

		// Check that the output matches the input
		success := true
		for i, inpd := range vals {
			if inpd == nil {
				continue
			}
			if !ff.EqualFr(inpd, &output[i]) {
				success = false
				break
			}
		}

		if !success {
			attempts += 1
			continue
		}
		// Output the evaluations if all good
		return output, nil
	}
	return nil, fmt.Errorf("max attempts reached: %d", attempts)
}
