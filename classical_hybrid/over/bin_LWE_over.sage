load get_runtime.sage
load binary_search.sage
load BKZ_cost.sage

    
	
	
def find_r_bin_LWE_over(n, q, m, k_vals = [4], min_delta = 1.002, max_delta = 1.03):
	global time_map
	time_map = {}

	# The determinant of the lattice
	det = q**(m-n)

	for r in k_vals:

		# The c_i are given by 2c_0 = 2c_1 = r/2
		# In fact, we only need c_0 to determine p_c
		c_0 = r/4

		# The dimension of the lattice reduction
		dim = m-r

		# The expected Euclidean norm Y of || v_l || in {-0.5,0.5}^{m-r}
		Y = sqrt(dim/4)

		# The probability that v_g has exactly 2c_0 entries equal to 0 and 2c_1 entries equal to 1
		# This will be used to determine success probability as p_succ = p_c * p_NP
		p_c = (binomial(r,2*c_0))/(2**r)

		# We assume that |S|=1 in the LWE case
		size_S = 1
		
		# x = delta
		#inner function with fixed r, c_0, det, dim, q, Y
		def f(x):
			f_x = RR(bkzcosts_mult_rounds(dim,x) - log(nr_loops_bin(r,c_0,dim, det, x, Y, q, size_S)*rt_NP_over(dim), 2))
			return f_x
		
		delta_r = find_zero(f, "x", min_delta, max_delta, value=0, args = {}, eps = 1, print_steps = false)
		bit_hardness_r = RR(bkzcosts_mult_rounds(dim,delta_r) + 1 - log(p_c * prob_NP(dim, det, Y, delta_r, q), 2))
		
		time_map[r] = (delta_r, bit_hardness_r)
		print r, time_map[r]
	return time_map

	
