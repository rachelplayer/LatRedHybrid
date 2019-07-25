# Uncomment this in order to obtain Wunderer's original "over" estimates
# gp.read("bkzsim.gp")

# Wunderer's original cost model for the "over" estimates
# Computes the log of the estimated cost of BKZ as a function of the blocksize beta, number of rounds and dimension of the lattice
def bkzoperations_over(beta,dim,rounds):
	return RR(0.187281*beta*log(beta, 2)-1.0192*beta+ log((dim)*rounds,2) +16.1)

# Wunderer's original function to determine the "over" estimates
# Find the log of the cost of BKZ in dimension dim to achieve target root Hermite factor delta_target
def bkzcosts_mult_rounds(dim, delta_target):

    # Find the minimal block size achieving delta_target according to Chen's thesis
	beta =36
	while RR(((((pi*beta)**(1/beta))*beta)/(2*pi*e))**(1/(2*(beta-1))))>RR(delta_target):
		beta=beta+1

	# Find number of rounds using simulator
	rounds = ZZ(gp.simulate(dim,min(k,dim),x)[3]) 

	# Return the runtime, computed using bkzoperations_over()
	return bkzoperations_over(min(k,dim),dim,rounds)

# Core-sieve cost model for BKZ with blocksize beta
# Returns the log of the estimated cost of BKZ under the core-sieve cost model
def bkz_operations_coresieve(beta):
	return RR(0.292*beta)

# Determine the log of the (core-sieve) cost of BKZ to achieve target root Hermite factor delta_target
def bkz_costs_coresieve(dim, delta_target):

    # Find the minimal block size achieving delta_target according to Chen's thesis
	beta =36
	while RR(((((pi*beta)**(1/beta))*beta)/(2*pi*e))**(1/(2*(beta-1))))>RR(delta_target):
		beta=beta+1

    # Return the runtime, computed using bkz_operations_coresieve()
	return bkz_operations_coresieve(min(beta,dim))








	