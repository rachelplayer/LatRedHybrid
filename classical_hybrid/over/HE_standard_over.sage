load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage

"""
This function is based on find_r_bin_LWE_over() in bin_LWE_over.sage
It takes an additional parameter: h, the Hamming weight of the secret
It can be used to estimate the cost of hybrid attack for typical HE parameters, 
according to Wunderer's analyis.
"""  
	
def find_r_HE_over(n, q, m, h, k_vals = [4], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}

	# Set the usual homomorphic encryption LWE standard deviation
	sigma = 3.2

	# Determine the scale factor nu for the Bai-Galbraith rebalancing
	nu = sqrt(n-r/h) * sigma

	# The determinant of the lattice
	det = nu**(n-r) * q**m

	for r in k_vals:

		c_minus1 = r * h / 4*n
		c_1 = r * h / 4*n

		# Define c_0 as the expected number of 0 entries in w'_g, a guess for half of w_g
		c_0 = r - c_minus1 - c_1

		# The dimension of the lattice reduction
		dim = m + n - r

		# The expected Euclidean norm Y of || w_l ||
		Y = sigma * sqrt(dim)

		# The probability p_c that w_g has exactly 2c_1 entries equal to 1 and 2c_minus1 entries equal to -1
		# This will be used to determine success probability as p_succ = p_c * p_NP
		p_c = 

		# We assume that |S|=1 in the LWE case
		size_S = 1
	
		# Denote delta by x 
		# Inner function with fixed r, c_1, c_minus1, det, dim, q, Y
		def f(x):
			f_x = RR(bkz_costs_coresieve(x) - log(nr_loops_HE(r, c1, c_minus1, dim, det, x, Y, q, size_S) * rt_NP_over(dim), 2))
			return f_x
	
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkz_costs_coresieve(delta_r) + 1 - log(p_c * prob_NP(dim, det, Y, delta_r, q), 2))
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map