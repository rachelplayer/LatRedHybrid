load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage

"""
This function is based on find_r_bin_LWE_over() in bin_LWE_over.sage
It takes an additional parameter: h, the Hamming weight of the secret
It can be used to estimate the cost of hybrid attack for typical HE parameters, 
according to Wunderer's analyis.
"""  
	
def find_r_HE_over(n, q, m, h, k_vals = [20], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}

	# Set the usual homomorphic encryption LWE standard deviation
	sigma = 3.2

	for r in k_vals:

		# Determine the scale factor nu for the Bai-Galbraith rebalancing
		nu = sqrt(n-r/h) * sigma

		# The determinant of the lattice
		det = nu**(n-r) * q**m

		# c_minus1 and c_1 are the expected number of -1 and 1 entries in w'_g, a guess for half of w_g
		# 2c_minus1 and 2c_1 are the expected number of -1 and 1 entries in w_g
		c_minus1 = round(r * h / 4*n)
		c_1 = round(r * h / 4*n)

		# Define c_0 as the expected number of 0 entries in w'_g, a guess for half of w_g
		c_0 = r - c_minus1 - c_1

		# We also need c_0_tilde, the expected number of 0 entries in w_g
		c_0_tilde = r - 2*c_1 - 2*c_minus1

		# The dimension of the lattice reduction
		dim = m + n - r

		# The expected Euclidean norm Y of || w_l ||
		Y = sigma * sqrt(dim)

		# To calculate p_c we'll need the expected number of nonzero entries in w_g
		h_wg = round(h * (r/n))

		# The probability p_c that w_g has exactly 2c_1 entries equal to 1 and 2c_minus1 entries equal to -1
		# This is equal to number of vectors with 2c_1 1s and 2c_minus1 -1s / number of ternary vectors with Hamming weight h
		# This will be used to determine success probability as p_succ = p_c * p_NP
		p_c = multinomial(2*c_minus1, c_0_tilde, 2*c_1) / (binomial(r, h_wg) * 2**h_wg)

		# We assume that |S|=1 in the LWE case
		size_S = 1
	
		# Denote delta by x 
		# Inner function with fixed r, c_1, c_minus1, det, dim, q, Y
		def f(x):
			# Default behaviour is to set probability p = 1, to calculate explicity set calc_p = True
			loops = nr_loops_HE(r, c_1, c_minus1, dim, det, x, Y, q, size_S)
			
			f_x = RR(bkz_costs_coresieve(dim, x) - log(loops * rt_NP_over(dim), 2))
			return f_x
	
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkz_costs_coresieve(dim, delta_r) + 1 - log(p_c * prob_NP(dim, det, Y, delta_r, q), 2))
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map