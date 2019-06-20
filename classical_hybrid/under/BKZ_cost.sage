# Wunderer's original cost model for the "under" estimates
# Computes the log of the estimated cost of BKZ as a function of the blocksize beta, number of rounds and dimension of the lattice
def bkzoperations(beta,dim,rounds):
	return RR(0.187281*beta*log(beta, 2)-1.0192*beta+ log((dim+1-beta)*rounds,2) +16.1)

# Wunderer's original function to determine the "under" estimates
# Find the log of the cost of BKZ in dimension dim to achieve target root Hermite factor delta_target
def bkzcosts_one_round(dim, delta_target):

	# Find the minimal block size achieving delta_target according to Chen's thesis
	beta = 36
	while RR(((((pi*beta)**(1/beta))*beta)/(2*pi*e))**(1/(2*(beta-1))))>RR(delta_target):
		beta=beta+1

	# Return the runtime, computed using bkzoperations()	
	return bkzoperations(min(beta,dim),dim,1)

# Core-sieve cost model for BKZ with blocksize beta
# Returns the log of the estimated cost of BKZ under the core-sieve cost model
def bkz_operations_coresieve(beta):
	return RR(0.292*beta)

# Determine the log of the (core-sieve) cost of BKZ to achieve target root Hermite factor delta_target
def bkz_costs_coresieve(delta_target):

    # Find the minimal block size achieving delta_target according to Chen's thesis
	beta =36
	while RR(((((pi*beta)**(1/beta))*beta)/(2*pi*e))**(1/(2*(beta-1))))>RR(delta_target):
		beta=beta+1

    # Return the runtime, computed using bkz_operations_coresieve()
	return bkz_costs_coresieve(beta)