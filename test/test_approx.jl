
module test_approx
using ApproXD
using Test


@testset "testing coefficient estimation on random data" begin

	d = Dict{Integer,Array{Float64,2}}()
	for i=1:3
		d[i] = rand(i+1,i+1)
	end
	 n1 = 2
	 n2 = 3
	 n3 = 4
	 # function values
	 y = zeros(n1*n2*n3)
	 for i=1:n1
	 	for j =1:n2
	 		for k=1:n3
	 			y[i + n1*(j-1 + n2*(k-1))] = d[1][i,1] + d[2][j,1] + d[3][k,1] + d[1][i,1]*d[2][j,1]* d[3][k,1]
	 		end
	 	end
	 end


	 ## build tensor product matrix
	 krons = kron(d[3],kron(d[2],d[1]))

	 ## notice: in kron the indexing corresponds to this:
	 ## uncomment to see
	 ##expand.grid( fastest=1:nrow(m1), middle=1:nrow(m2), slowest=1:nrow(m3) )

	 # check for equality
	 res1 = krons * y
	 myres = getTensorCoef(d,y)

	 @test isapprox(sum(abs.(res1 - myres)) , 0.0, atol=1e-10)

end




@testset "testing computation of coefficients" begin

	# 3D approximation example

	# bounds
	lb = [-1.1,2.5,0.1]
	ub = [3.2,5.4,9.1]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [12,9,8]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,2,1]

	# implies a number of knots for each spline
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:3)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:3)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:3)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:3
		d[i] = Array(getBasis(points[i],bsp[i]))
	end
	ds = Dict()
	for i=1:3
		ds[i] = getBasis(points[i],bsp[i])
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x,y,z)
		(x+y).^2 + z.^z
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j,k) for i in points[1] , j in points[2], k in points[3]]

	yvec = y[:]

	# compute coefs by constructing tensor product
	# where note the order is important!
	 ikrons = kron(id[3],kron(id[2],id[1]))
	 coef1 = ikrons * yvec

	 # get coefs using the function
	 mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 5e-6

	 # check coefs are the same
	@test isapprox(maximum(abs,coef1 - mycoef), 0.0, atol=tol)

	# approximate the function values on the original grid
	# by using the basis in d. this is just reverse of getting coefs!
	pred = getTensorCoef(d,mycoef)
	@test isapprox(maximum(abs,pred - yvec) , 0.0 , atol=tol)

	t1 = reshape(pred,npoints[1],npoints[2],npoints[3])
	@test isapprox(maximum(abs,t1 .- y) , 0.0 , atol=tol)

	# is that really doing what you want?
	# you want B * coefs
	B = kron(d[3],kron(d[2],d[1]))
	pred2 = B * mycoef
	@test isapprox(maximum(abs,pred2 - yvec) , 0.0 , atol=tol)

	# predict usign the predict function
	print(typeof(ds[3]))
	pred = ApproXD.evalTensor3(ds[3],ds[2],ds[1],mycoef)
	@test isapprox(maximum(abs,pred - yvec) , 0.0 , atol=tol)

end


@testset "testing 1D spline evaluating off grid" begin


	ndims = 1

	# bounds
	lb = [-2.1]
	ub = [2.2]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [11]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x)
		x^3
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i) for i in points[1]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 1e-2

	# predict values off grid


	rval1 = lb[1] + 0.03
	b1 = getBasis(rval1,bsp[1])
	@test isapprox((mycoef' * b1)[1] , f(rval1), atol=1e-7)

end

@testset "testing 1D spline extrapolation" begin


	ndims = 1

	# bounds
	lb = [-2.1]
	ub = [2.2]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [11]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x)
		x^3
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i) for i in points[1]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 1e-7

	# predict values out of grid
	rval1 = ub[1] + eps()
	b1 = getBasis(rval1,bsp[1])
	@test isapprox((mycoef' * b1)[1], f(rval1),atol=tol)

	rval1 = ub[1] + 1.0
	b1 = getBasis(rval1,bsp[1])
	@test isapprox((mycoef' * b1)[1], f(rval1),atol=tol)
	# @test (mycoef' * b1)[1] --> roughly(f(rval1),atol=tol)

	rval1 = ub[1] + 2.0
	b1 = getBasis(rval1,bsp[1])
	@test isapprox((mycoef' * b1)[1], f(rval1),atol=tol)
	# @test (mycoef' * b1)[1] --> roughly(f(rval1),atol=tol)

end


@testset "testing 2D tensorProduct evaluating off grid" begin

	# 2D approximation example

	ndims = 2

	# bounds
	lb = [-2.1,-2.5]
	ub = [2.2,2.6]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [11,12]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i], length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x,y)
		sin(sqrt(x^2+y^2))
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j) for i in points[1], j in points[2]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 1e-2

	# predict values off grid


	rval1 = lb[1] + 0.03
	rval2 = lb[2] + 0.21
	b2 = getBasis(rval2,bsp[2])
	b1 = getBasis(rval1,bsp[1])
	@test isapprox(ApproXD.evalTensor2(b2,b1,mycoef), f(rval1,rval2),atol=tol)

	rval1 = ub[1] - 0.3
	rval2 = ub[2] - 0.11
	b2 = getBasis(rval2,bsp[2])
	b1 = getBasis(rval1,bsp[1])
	@test isapprox(ApproXD.evalTensor2(b2,b1,mycoef) ,f(rval1,rval2),atol=tol)

end


@testset "testing 3D tensorProduct approximations" begin

	# 3D approximation example

	ndims = 3

	# bounds
	lb = [-2.1,-2.5,0.9]
	ub = [2.2,2.6,5.0]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [12,16,13]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end


	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x,y,z)
		sin(sqrt(x^2+y^2)) - z^2
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j,k) for i in points[1], j in points[2], k in points[3]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)
	 # testing tolerance
	 tol = 1e-2

	# predict values off grid

	rval1 = lb[1] + 0.09
	rval2 = lb[2] + 0.99
	rval3 = lb[3] + 2.01
	b3 = Array(getBasis(rval3,bsp[3]))[:]
	b2 = Array(getBasis(rval2,bsp[2]))[:]
	b1 = Array(getBasis(rval1,bsp[1]))[:]
	@test isapprox(ApproXD.evalTensor3(b3,b2,b1,mycoef),f(rval1,rval2,rval3),atol=tol)

	rval1 = ub[1] - 0.9
	rval2 = ub[2] - 0.09
	rval3 = ub[3] - 1.01
	b3 = Array(getBasis(rval3,bsp[3]))[:]
	b2 = Array(getBasis(rval2,bsp[2]))[:]
	b1 = Array(getBasis(rval1,bsp[1]))[:]
	@test isapprox(ApproXD.evalTensor3(b3,b2,b1,mycoef)[1],f(rval1,rval2,rval3),atol=tol)


end


@testset "testing 4D tensorProduct approximations" begin

	ndims = 4

	# bounds
	lb = [-1.1,-1.5,-0.9,-1.0]
	ub = [1.2,1.6,0.9,1]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [13,13,13,13]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3,3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x,y,z,w)
		sin(sqrt(x^2+y^2)) + (z-w)^3
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 5e-2

	# predict values off grid


	rval1 = lb[1] + 0.14
	rval2 = lb[2] + 1.01
	rval3 = lb[3] + 0.33
	rval4 = lb[4] + 0.45
	b4 = getBasis(rval4,bsp[4])
	b3 = getBasis(rval3,bsp[3])
	b2 = getBasis(rval2,bsp[2])
	b1 = getBasis(rval1,bsp[1])
	@test isapprox(ApproXD.evalTensor4(b4,b3,b2,b1,mycoef)[1] ,f(rval1,rval2,rval3,rval4),atol=tol)

	rval1 = ub[1] - 0.14
	rval2 = ub[2] - 1.01
	rval3 = ub[3] - 0.33
	rval4 = ub[4] - 0.45
	b4 = getBasis(rval4,bsp[4])
	b3 = getBasis(rval3,bsp[3])
	b2 = getBasis(rval2,bsp[2])
	b1 = getBasis(rval1,bsp[1])
	@test isapprox(ApproXD.evalTensor4(b4,b3,b2,b1,mycoef)[1] ,f(rval1,rval2,rval3,rval4),atol=tol)

end

@testset "testing getTensorCoef performance on 4D" begin

	ndims = 4

	# bounds
	lb = [-1.1,-1.5,-0.9,-1.0]
	ub = [1.2,1.6,0.9,1]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [3,3,3,17]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [1,1,1,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	# #  get a function
	# function f(x,y,z,w)
	# 	sin(sqrt(x^2+y^2)) + (z-w)^3
	# end

	# # compute function values so that
	# # k is fastest varying index
	# y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

	# yvec = y[:]

	yvec = rand(prod(npoints))
	yvec2 = rand(prod(npoints))
	yvec3 = rand(prod(npoints))
	yvec4 = rand(prod(npoints))

	# get coefs using the function
	t0 = time()
	for i in 1:(2^4 * 81 * 29)	# is,ihh,ih,itau,j,k,age
		mycoef = getTensorCoef(id,yvec);
		mycoef = getTensorCoef(id,yvec2);
		mycoef = getTensorCoef(id,yvec3);
		mycoef = getTensorCoef(id,yvec4);
	end
	@info("4D timing: $(time()-t0)")
end

if Sys.isapple()
	@testset "testing getTensorCoef! performance on 4D" begin

		ndims = 4

		# bounds
		lb = [-1.1,-1.5,-0.9,-1.0]
		ub = [1.2,1.6,0.9,1]

		# number of eval points and basis functions:
		# we require square basis matrices!
		npoints = [3,3,3,17]

		# number of basis funcs
		nbasis = npoints

		# splien degrees
		degs = [1,1,1,3]

		# implies a number of knots for each spline
		# remember the restriction that nknots == ncoefs
		nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

		# eval points
		points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

		# set up ApproXD
		bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

		# set of basis functions
		d = Dict{Integer,Array{Float64,2}}()
		for i=1:ndims
			d[i] = Array(getBasis(points[i],bsp[i]))
		end

		# set of INVERSE basis functions
		id = Dict{Integer,Array{Float64,2}}()
		for k in collect(keys(d))
			id[k] = inv(d[k])
		end

		# #  get a function
		# function f(x,y,z,w)
		# 	sin(sqrt(x^2+y^2)) + (z-w)^3
		# end

		# # compute function values so that
		# # k is fastest varying index
		# y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

		# yvec = y[:]

		yvec = rand(prod(npoints))
		yvec2 = rand(prod(npoints))
		yvec3 = rand(prod(npoints))
		yvec4 = rand(prod(npoints))

		# get coefs using the function
		t0 = time()
		for i in 1:(2^4 * 81 * 29)	# is,ihh,ih,itau,j,k,age
			getTensorCoef(id,yvec);
			getTensorCoef(id,yvec2);
			getTensorCoef(id,yvec3);
			getTensorCoef(id,yvec4);
		end
		@info("4D timing preallocated: $(time()-t0)")
	end



	@testset "testing getTensorCoef performance on 5D" begin

		ndims = 5

		# bounds
		lb = [-1.1,-1.5,-0.9,-1.0,1]
		ub = [1.2,1.6,0.9,1,9]

		# number of eval points and basis functions:
		# we require square basis matrices!
		npoints = [3,3,3,6,7]

		# number of basis funcs
		nbasis = npoints

		# splien degrees
		degs = [1,1,1,2,2]

		# implies a number of knots for each spline
		# remember the restriction that nknots == ncoefs
		nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

		# eval points
		points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

		# set up ApproXD
		bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

		# set of basis functions
		d = Dict{Integer,Array{Float64,2}}()
		for i=1:ndims
			d[i] = Array(getBasis(points[i],bsp[i]))
		end

		# set of INVERSE basis functions
		id = Dict{Integer,Array{Float64,2}}()
		for k in collect(keys(d))
			id[k] = inv(d[k])
		end

		# #  get a function
		# function f(x,y,z,w)
		# 	sin(sqrt(x^2+y^2)) + (z-w)^3
		# end

		# # compute function values so that
		# # k is fastest varying index
		# y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

		# yvec = y[:]

		yvec = rand(prod(npoints))
		yvec2 = rand(prod(npoints))
		yvec3 = rand(prod(npoints))
		yvec4 = rand(prod(npoints))

		# get coefs using the function
		t0 = time()
		for i in 1:(2^4 * 9 * 29)	# is,ihh,ih,itau,j,age
			mycoef = getTensorCoef(id,yvec);
			mycoef = getTensorCoef(id,yvec2);
			mycoef = getTensorCoef(id,yvec3);
			mycoef = getTensorCoef(id,yvec4);
		end
		@info("5D timing: $(time()-t0)")

	end

	@testset "testing getTensorCoef performance on 10D" begin


		# bounds
	#      (nJ, ns, nz, ny, np, na, nh, ntau,  nJ, nt-1 )
		lb = [1,1,-1,0,3.0, -1.1,0,1,1]
		ub = [9,2, 1,3,9,6.1,1  ,2,2,29]


		ndims = length(lb)

		# number of eval points and basis functions:
		# we require square basis matrices!
		npoints = [9,3,3,3,3,17,3,3,29]

		# number of basis funcs
		nbasis = npoints

		# splien degrees
		degs = [1,1,1,1,1,3,1,1,1]

		# implies a number of knots for each spline
		# remember the restriction that nknots == ncoefs
		nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

		# eval points
		points = Dict(i => collect(range(lb[i],stop = ub[i],length = npoints[i])) for i=1:ndims)

		# set up ApproXD
		bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

		# set of basis functions
		d = Dict{Integer,Array{Float64,2}}()
		for i=1:ndims
			d[i] = Array(getBasis(points[i],bsp[i]))
		end

		# set of INVERSE basis functions
		id = Dict{Integer,Array{Float64,2}}()
		for k in collect(keys(d))
			id[k] = inv(d[k])
		end

		# #  get a function
		# function f(x,y,z,w)
		# 	sin(sqrt(x^2+y^2)) + (z-w)^3
		# end

		# # compute function values so that
		# # k is fastest varying index
		# y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

		# yvec = y[:]

		yvec = rand(prod(npoints))

		# get coefs using the function
		t0 = time()
		for i=1:9
			mycoef = getTensorCoef(id,yvec);
		end
		@info("10D timing: $(time()-t0)")



	end
else
	@info("I am skipping performance tests on travis because too much memory required.")
end


end
