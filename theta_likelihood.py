from math import log

#recursive function that calculates P(S|theta) for a given sample size (from Hudson 1990)
def P(n, S, T):
	if n==2:
		return ((T/(1. + T))**S * 1/(1. + T))
	else:
		SUM = 0
		for i in range(S+1):
			if not P_calls.has_key((n-1, S-i, T)):
				P_calls[(n-1, S-i, T)] = P(n-1, S-i, T)
			if not Q_calls.has_key((n, i, T)):
				Q_calls[(n, i, T)] = Q(n, i, T)
			SUM+= P_calls[(n-1, S-i, T)] * Q_calls[(n, i, T)]
		return SUM

#function to calculate Q (from Hudson 1990)
def Q(n, S, T):
	return (T/(T+n-1.))**S * (n-1.)/(T+n-1)

#makes a range of theta values to test 
def theta_list_maker(start, end, increment):
	while start <= end:
		theta_list.append(start)
		start += increment



#BRA inputs key = gene name, value = list (Length, segregating sites, sample size)
inputs = {'EP0001': [418, 16, 278], 'EP0141':[530, 54, 278], 'EP0143':[717, 17, 278], 'EP0144':[245, 29, 278], 'EP0188':[291, 4, 278], 'EP0222':[387, 4, 278], 'EP0267':[367, 18, 278], 'EP0285':[677, 29, 278], 'EP0314':[538, 27, 278], 'EP0317':[540, 17, 278]}


#sort the keys into a list
input_keys = inputs.keys()
input_keys.sort()

#make theta list
theta_list = []
theta_list_maker(.00005, .05, .00002)

results = {'theta': theta_list}
outfile2 = open('theta_like.csv', 'w')
#Goes through each gene in the input and passes the values into function P for every value in theta_list
for key in input_keys:
	#outfile = open('bra_' + key + '.csv', 'w')
	#outfile.write('theta' + ',' + key + "\n") # write header lined
	results[key] = []
	sums = {}
	L = inputs[key][0]
	S = inputs[key][1]
	n = inputs[key][2]


	for T in theta_list:
		t_len = T * L	#converts theta/bp for each gene to theta total
		P_calls = {}	#necessary for memoP
		Q_calls = {}	#necessary for Q
		likelihood = P(n,S,t_len) #takes the ln of the likelihood calculated
		results[key].append(likelihood)
		#outfile.write(str(T) + ',' + str(likelihood) + "\n") #write to output

results_keys = results.keys()
results_keys.sort()

for key in results_keys:
	outfile2.write(key + ',')

outfile2.write("\n")



j = 0
while j < len(results['theta']):
	for key in results_keys:
		outfile2.write(str(results[key][j]) + ',')
	outfile2.write("\n")
	j += 1




