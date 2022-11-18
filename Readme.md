# go-poly

Fork of the code from [go-kzg](github.com/protolambda/go-kzg). Currently works with [go-mcl](github.com/alinush/go-mcl).

## List of features
- FFT
- Polynomial operations
    - Mul
    - xGCD
    - Div
    - Subproduct tree

## To do
- [ ] Add gurvy
- [ ] Add back kilic
- [ ] Add back herumi/bls-eth-go-binary

## Run benchmarks

```bash
go test ./... -bench=. -run=^a  -timeout 240m
time go test ./... -bench=. -run=^a -benchtime=100x -timeout 240m
time go test ./... -bench=BenchmarkPoly -run=^a -benchtime=100x  -timeout 240m
```
