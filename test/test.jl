include("../src/includes.jl")


function main()
	# generate random data
	L = 5
	d = 2
	N = 100
	data = [rand(1:2, L) for i in 1:N]

	mps = cross_approximation(data, [d for i in 1:L])
end