import galoisfield as gf
import random

class Decoder():
	
	def __init__(self, m ,t):
		self.m = m
		self.t = t
		self.field = dict()
		self.n = pow(2,m) - 1
		self.k = self.n - 2*t
		# self.gen_poly = list()
		self.GenerateGF()

	# create galois field (dictionary of alfas)
	def GenerateGF(self):
		g = gf.GaloisField(self.m, self.t)
		self.field = g.GenerateGF()
		self.char_dict = g.GetCharDict()

	# generate generating polynomial (x+a^1)(x+a^2)...(x+a^2t)
	def Generate_poly(self):
		n = []
		for i in range(self.t * 2):
			n.append([0, i+1])
		res = n[0]
		for i in range(self.t * 2 - 1):
			res = self.Multiply_poly(res, n[i+1])
		return res



	#															#
	#															#
	#		[	convert string of bits to polynomial	]		#
	#															#	
	#															#



	# create one list that contains message (divided by m bits) in bits
	# ['10000', '01001']
	def SplitMessageToMBits(self, message):
		msg = []
		for i in range(len(message)):
			# truncate message by every m elements
			if i % self.m == 0:
				msg.append(message[i:i+self.m])
		# "010" -> "00010" adjust to m bits at front
		if len(msg[-1]) != self.m:
			msg[-1] = "0" * (self.m - len(msg[-1])) + msg[-1]
		return msg

	# create list of lists that every list contains message (divided by m bits) in bits
	# len(message) <= m*k -> [['10000', '01001']]
	# len(message) >= m*k -> [['10000', '01001', ... ], ['10011', '11111']]
	# by default this method is used instead of 'SplitMessageToMBits()'
	def SplitMessage(self, message):
		message = message.replace(" ", "")
		if len(message) <= (self.n * self.m):
			splited_message = []
			splited_message.append(message[:(self.n * self.m)])
			splited_message[0] = self.SplitMessageToMBits(splited_message[0])
		else:
			splited_message = []
			for i in range(len(message)):
				if i % (self.n * self.m) == 0:
					splited_message.append(message[i:i+(self.n * self.m)])
			for i in range(len(splited_message)):
				splited_message[i] = self.SplitMessageToMBits(splited_message[i])
		return splited_message

	# code m bits values in list to alfa symbols using dictionary mapping
	# ['10000', '00101'] -> [4, 5]
	def MessageToPoly(self, msg):
		poly = []
		for message in msg:
			poly.append(self.Return_key(message))
		return poly

		# add zeros in front of message to make a packet of size m*n appropriate to send
	def AdjustMessage(self, message_in_bits):
		if len(message_in_bits) < (self.m * self.k):
			message_in_bits = "0" * ((self.m * self.n) - len(message_in_bits)) + message_in_bits
		return message_in_bits

	# convert polynomials to bits
	# [[11, 16, ...], [0, 4, ...]] -> "10110100101100101..."
	def PolyToBits(self, poly):
		string = str()
		str_tmp = str()
		list_tmp = list()
		for i in poly:
			str_tmp = str_tmp + self.field[i]
		str_tmp = self.AdjustMessage(str_tmp)
		list_tmp.append(str_tmp)
		str_tmp = ''
		for i in list_tmp:
			string += i
		return string



	#											#
	#											#
	#		[	arithmetic operations	]		#
	#											#
	#											#



	def XORonPolynomials(self, first_poly, second_poly):
		smaller = []
		product = []
		f_poly_len = len(first_poly) - 1
		s_poly_len = len(second_poly) - 1
		if f_poly_len <= s_poly_len:
			smaller = first_poly
			product = second_poly[:]
		else:
			smaller = second_poly
			product = first_poly[:]
		prod_len = len(product) - 1
		for i in range(len(smaller)):
			product[prod_len - i] = self.XOR(self.field[first_poly[f_poly_len - i]], self.field[second_poly[s_poly_len - i]])
		return product


	# addition in Galoa field (addition of polynomial's coefficients)
	# a xor a -> None
	# None xor a -> a
	def XOR(self, a, b):
		c = ''
		for x in range(len(a)):
			if a[x] != b[x]:
				c = c+'1'
			else:
				c = c+'0'
		return self.Return_key(c);

	def XOR_on_bits(self, a, b):
		c = ''
		for x in range(len(a)):
			if a[x] != b[x]:
				c = c+'1'
			else:
				c = c+'0'
		return c;

	# get key from value in galoa field dictionary
	def Return_key(self, n):
		for key, value in self.field.items():
			if n == value:
				return key;

	# reduce None type at front of polynomial
	# [None, None, 4, 2] -> [4, 2]
	def Adjust_poly(self, poly):
		while poly[0] == None:
			poly.pop(0)
			if poly == []:
				break
		return poly


	# polynomials multiplication
	# a = [4,None, 3] -> 4x^2 + 0x + 3
	def Multiply_poly(self,a,b):
		# self.Adjust_poly(a)
		# self.Adjust_poly(b)
		a = a[:]
		b = b[:]
		if(a == [] or b == []):
			return [None]
		product = [None] * (len(a) + len(b) - 1)
		for i in range(len(a)):
			for j in range(len(b)):
				if a[i] == None or b[j] == None:
					product[i+j] = self.XOR(self.field[product[i+j]], self.field[None])
				else:
					coef = (a[i] + b[j]) % (pow(2,self.m) - 1) # m = 5 -> % 31
					product[i+j] = self.XOR(self.field[product[i+j]], self.field[coef])
		return product


	# polynomial division
	# a = [4,None, 3] -> 4x^2 + 0x + 3
	def Divide_poly(self, dividend, divisor):
		if divisor == [None]:
			# division by 0 == exception?
			return [None]

		if len(dividend) < len(divisor):
			return [None]

		dividend = dividend[:]
		divisor = divisor[:]

		self.Adjust_poly(dividend)
		self.Adjust_poly(divisor)

		shift = len(dividend) - len(divisor)
		quotient = []
		while len(dividend) >= len(divisor):
			mult = (dividend[0] - divisor[0]) % (pow(2,self.m) - 1)
			quotient.append(mult)
			for index in range(len(divisor)):
				# 0 - alfa^n = alfa^n
				if divisor[index] == None or mult == None:
					dividend[index] = self.XOR(self.field[dividend[index]], self.field[None])
				else:
					coef = (mult + divisor[index]) % (pow(2,self.m) - 1)
					dividend[index] = self.XOR(self.field[dividend[index]], self.field[coef])
			self.Adjust_poly(dividend)
		# case when few last divisions are zeros
		if len(quotient) <= shift:
			quotient.append(None)
		return quotient

	# polynomial modulo
	# a = [4,None, 3] -> 4x^2 + 0x + 3
	def Modulo_poly(self, dividend, divisor):
		if divisor == [None]:
			# division by 0 == exception?
			return [None]

		if len(dividend) < len(divisor):
			return [None]

		dividend = dividend[:]
		divisor = divisor[:]

		self.Adjust_poly(dividend)
		self.Adjust_poly(divisor)

		shift = len(dividend) - len(divisor)
		quotient = []
		while len(dividend) >= len(divisor):
			mult = (dividend[0] - divisor[0]) % (pow(2,self.m) - 1)
			quotient.append(mult)
			for index in range(len(divisor)):
				if divisor[index] == None or mult == None:
					dividend[index] = self.XOR(self.field[dividend[index]], self.field[None])
				else:
					coef = (mult + divisor[index]) % (pow(2,self.m) - 1)
					dividend[index] = self.XOR(self.field[dividend[index]], self.field[coef])
			self.Adjust_poly(dividend)
			# print("dividen", dividend)
		return dividend



	# nie wiadomo czy będzie to potrzebne !!! --------
	def BinaryStringToInt(self, string):
		tmp_list = []
		for i in range(len(string)):
			if string[i] == '1':
				tmp_list.append(1)
			elif string[i] == '0':
				tmp_list.append(0)
		value = 0
		p = 0
		tl = len(tmp_list)
		for i in reversed(range(tl)):
			value += tmp_list[i] * (pow(2, p))
			p += 1
		return value




	#							#
	#							#
	#	[	decode message	]	#
	#							#
	#							#


	def DecodeMessage(self, message):
		self.GenerateGF()
		buffer = self.SplitMessage(message)
		# print(buffer)
		message_poly = self.GetMessagePoly(buffer)
		# print(message_poly)
		res = list()
		for poly in message_poly:
			res.append(self.SimpleDecodeAlgorithm(poly))
		bits = str()
		for i in res:
			bits += self.PolyToBits(i)
		return bits



	# czy getreminder jest potrzebne ???	-----------------	JEST !

	# calculate reminder of received message polynomial
	# return [[11, 16, 8, ...], [0, 4, None, ...]] or [[0, 4, ...]]
	# depending on message length
	def GetReminderOfMessagePoly(self, message_poly):
		gen_poly = self.Generate_poly()
		# print("gen_poly:", gen_poly)
		rem_poly = self.Modulo_poly(message_poly, gen_poly)
		# print("rem_poly: ", rem_poly)
		return rem_poly

	# convert alfas to list with polynomial coefficients in them
	# code m bits values in lists to alfa symbols
	# [['10000', '00101', ... ], ['10011', '11111']] -> [[4, 5, ...], [12, 15]]
	# used by default instead of 'MessageToPoly()'
	def GetMessagePoly(self, buffer):
		msg_poly = []
		for i in buffer:
			msg_poly.append(self.MessageToPoly(i))
		return msg_poly



	#							#
	#							#
	#	[	simple decoder	]	#
	#							#
	#							#



	def SimpleDecodeAlgorithm(self, polynomial):
		# print("polynomial origin:", polynomial)
		i = 0
		while True:
			syndrome = self.GetReminderOfMessagePoly(polynomial)
			w = self.WeightOfSyndrome(syndrome)
			# print("syn: ", syndrome)
			# print("w: ", w)
			if w <= self.t:
				polynomial = self.XORonPolynomials(polynomial, syndrome)
				polynomial = self.MoveITimesLeft(polynomial, i)
				# print("poly:", polynomial)
				# print("syndrome", syndrome)
				# print("skorygowany wektor")
				return polynomial
			else:
				# if i == self.n:
				if i == self.k:
					# print("bład niekorygowalny")
					return polynomial
				else:
					polynomial = self.MoveITimesRight(polynomial, 1)
					# print("poly:", polynomial)
					i += 1
					# print("i:", i)
		return polynomials

	# calculate weight of syndrome -> how many symbols was changed
	def WeightOfSyndrome(self, polynomial):
		weight = 0
		for i in polynomial:
			if i != None:
				weight += 1
		return weight

	# move element at last index to front
	def MoveITimesRight(self, poly, i):
		for n in range(i):
			poly.insert(0, poly.pop())
		return poly

	# move element at first index to end
	def MoveITimesLeft(self, poly, i):
		for n in range(i):
			poly.append(poly.pop(0))
		return poly




	#							#
	#							#
	#		[	tests	]		#
	#							#
	#							#


	def TestCorrectionEffectiveness(self, message, error_number, loop_number):
		self.GenerateGF()
		buffer = self.SplitMessage(message)
		message_poly = self.GetMessagePoly(buffer)
		# print(message_poly)
		# print("-------------------------------")
		# for i in message_poly:
		# 	for j in range(10):
		# 		res = self.MakeErrors(i, 1)
		# 		print(res)
		# print("-------------------------------")
		reminder_poly = self.GetReminderOfMessagePoly(message_poly)
		res = list()
		for poly in message_poly:
			negative_res = 0
			for i in range(loop_number):
				msg = self.MakeErrors(poly, error_number)
				one_res = self.TestSimpleDecodeAlgorithm(msg)
				if one_res == 1:
					negative_res += 1
			print("Not corrected:", negative_res)
			print("Corrected:", loop_number - negative_res)
			# res.append(self.SimpleDecodeAlgorithm(poly))







	# moduł random to jedyny moduł, który został zaimportowany z zewnętrznej biblioteki, ponieważ
	# własny generator liczb losowych nie był losowy
	def MakeErrors(self, message_poly, error_number):

		allowed_indexes = set()
		while len(allowed_indexes) != error_number:
			allowed_indexes.add(random.randrange(self.n - 1))
		allowed_indexes = list(allowed_indexes)
		random.shuffle(allowed_indexes)

		message_poly = message_poly[:]
		for i in range(error_number):
			rand_symbol = random.randrange(self.n - 1)
			while True:
				if message_poly[allowed_indexes[i]] == rand_symbol:
					rand_symbol = random.randrange(self.n - 1)
				else:
					message_poly[allowed_indexes[i]] = rand_symbol
					break
		return message_poly


	def MakeErrors_burst(self, message_poly, error_lenght):
		error_vector = [None] * self.n
		start = random.randrange(self.n - error_lenght - 1)
		end = start + error_lenght - 1
		rand_symbol = random.randrange(self.n - 1)
		if message_poly[start] == rand_symbol:
			while True:
				if message_poly[start] == rand_symbol:
					rand_symbol = random.randrange(self.n - 1)
				else:
					message_poly[start] = rand_symbol
					break
		else:
			message_poly[start] = rand_symbol
		rand_symbol = random.randrange(self.n - 1)
		if message_poly[end] == rand_symbol:
			while True:
				if message_poly[end] == rand_symbol:
					rand_symbol = random.randrange(self.n - 1)
				else:
					message_poly[end] = rand_symbol
					break
		else:
			message_poly[end] = rand_symbol
		for i in range(start+1, end):
			rand_symbol = random.randrange(self.n - 1)
			error_vector[i] = rand_symbol
		error_polynomial = self.XORonPolynomials(message_poly, error_vector)
		return error_polynomial


	def MakeErrors_bits(self, message_bits, error_number):

		allowed_indexes = set()
		while len(allowed_indexes) != error_number:
			allowed_indexes.add(random.randrange((self.n * self.m) - 1))
		allowed_indexes = list(allowed_indexes)
		random.shuffle(allowed_indexes)

		for i in range(error_number):
			# range [0, 1]
			rand_symbol = random.randrange(2)
			while True:
				if message_bits[allowed_indexes[i]] == str(rand_symbol):
					rand_symbol = random.randrange(2)
				else:
					print(rand_symbol)
					print(allowed_indexes[i])

					print(rand_symbol)
					message_bits = message_bits[:allowed_indexes[i] - 1] + str(rand_symbol) + message_bits[allowed_indexes[i]:]
					break
		return message_bits



		# wstawia na n-tej i n + długość błędu (z definicji mochnackiego str. 117) pozycji inny bit, a następnie losowo zmienia bity między krańcami błędu grupowego
	def MakeErrors_burst_bits(self, message_bits, error_lenght):
		error_vector = "0" * (self.n * self.m)
		start = random.randrange((self.n * self.m) - error_lenght - 1)
		end = start + error_lenght - 1
		if message_bits[start] == "0":
			start_bit = "1"
		else:
			start_bit = "0"
		if message_bits[end] == "0":
			end_bit = "1"
		else:
			end_bit = "0"
		error_vector = error_vector[:start-1] + start_bit + error_vector[start:]
		error_vector = error_vector[:end-1] + end_bit + error_vector[end:]
		for i in range(start+1, end-1):
			rand_bit = random.randrange(2)
			error_vector = error_vector[:i - 1] + str(rand_bit) + error_vector[i:]
		error_message = self.XOR_on_bits(message_bits, error_vector)
		return error_message






	def TestSimpleDecodeAlgorithm(self, polynomial):
		copy = polynomial[:]
		i = 0
		while True:
			syndrome = self.GetReminderOfMessagePoly(polynomial)
			w = self.WeightOfSyndrome(syndrome)
			if w <= self.t:
				print(f"copy: {copy}\npoly: {polynomial}\n")
				polynomial = self.XORonPolynomials(polynomial, syndrome)
				polynomial = self.MoveITimesLeft(polynomial, i)
				print(f"wynik: {polynomial}")
				print("-----------------------------------------------------------")
				return 0
			else:
				# if i == self.n:
				if i == self.k:
					# print("bład niekorygowalny")
					return 1
				else:
					polynomial = self.MoveITimesRight(polynomial, 1)
					i += 1
		# 0 to jest dobrze i nie ma dodawać do liczby błędnych poprawień w test correction effectiveness
		# 1 to źle i ma dodawać do liczby =======||==========





######################################################################################################################################################



	#												#
	#												#
	#		[	Berlekamp-Massey algorithm    ]		#
	#												#
	#												#


	# jeśli wynikiem jest [0] to znaczy, że wiadomość nie ma błędów
	#
	# ostatni element np. w lokalizatorze błędów [18, 20, 9, 10, 1, 0] musi być 0, czyli lambda_0 = 1
	#
	# syndromy (argument) muszą być ułożone ...s4, s3, s2, s1
	# czyli tak jak zwraca modulo, ale algorytm BM obraca je i działa na syndromach s1, s2, s3, s4...
	def Berlekamp_Massey(self, syndrome):
		# syndrome = [11, 29, 11, 5, 21, 20, 1, 27, 13, 11]
		# syndrome = [4, None, None, None]
		syndrome = self.AdjustReminder(syndrome)
		# syndrome = [6,9,7,3,6,4,0,3]
		# syndrome.reverse()
		# self.t = 4
		# self.m = 4
		# self.n = 15
		# self.k = 11
		# self.GenerateGF()

		print(self.field)
		# print(self.XOR(self.field[25], self.field[20]))

		Lr = 0
		r = 0

		error_locator = [0]
		B = [0]

		# while r < 2*self.t:
		while r <= 2*self.t:
			r += 1

			discrepancy = self.Discrepancy(error_locator, syndrome, r, Lr)
			# print("discrepancy:", discrepancy)

			if discrepancy == None:
				Lr = Lr

			else:

				# print(f"Lr = {Lr}")

				tmp_B = B[:]
				# delta * B(x)
				for i in range(len(tmp_B)):
					tmp_B[i] = self.Berlekamp_Massey_multiplication(tmp_B[i], discrepancy)

				# delta * xB(x)
				tmp_B.append(None)
				# print(f"Delta * xB(x): {tmp_B}")

				copy_error_loc = error_locator[:]

				# print(f"copy_error_loc = {error_locator}")
				# print(f"tmp_B = {tmp_B}")
				error_locator = self.XORonPolynomials(copy_error_loc, tmp_B)
				# error_locator = self.Multiply_poly(copy_error_loc, tmp_B)
				# print(f"error_locator^{r}: {error_locator}")
				
				# print(f"{2*Lr} <= {r-1}")
				if 2*Lr <= r - 1:
					Lr = r - Lr
					# print("Lr if", Lr)
					# B(x) = lambda / delta
					B = self.Divide_poly(copy_error_loc, [discrepancy])
					# print("B:", B)
					# print("copy error_locator", copy_error_loc)
				else:
					Lr = Lr
					B.append(None)
					# print("B:", B)
					# print("copy error_locator", copy_error_loc)

			# step 5
			self.Adjust_poly(error_locator)
			# if len(error_locator) > self.t:
			if Lr > self.t:

				print("\nFinish. Non correctable error occured")
				# print(f"error locator polynomial: {error_locator}")
				return error_locator

			if r == 2*self.t:
				print("\nSuccess !!!")
				# print(f"error locator polynomial: {error_locator}")
				return error_locator

			# print("-----------------------------------------------")




	def Discrepancy(self, error_locator_poly, syndrom_poly, r, Lr):
		# print("r:",r)
		error_locator = error_locator_poly
		syndrome = syndrom_poly
		# print(f"Error: {error_locator}\nSyndrome: {syndrome}")
		sum_poly = []

		if r == 1:
			first_discrepancy = self.Berlekamp_Massey_multiplication(error_locator[0], syndrome[0])
			return first_discrepancy
			# print(f"first_discrepancy: {first_discrepancy}")
			# sum_poly.append(first_discrepancy)
		else:
			l = len(error_locator) - 1
			for j in range(1, Lr+1):
				# print("j = ", j)
				# print("Lr+1 = ", Lr+1)

				# tmp = error_locator[iteration - j - 1] * syndrome[iteration - j] % self.n
				# print("error_locator:", error_locator)
				# print(f"Lr = {Lr} and j = {j}")
				# print(f"error-locator length = {error_locator}") 
				# print("error_locator[r - j]:", error_locator[l - j])


				# print("S[i]:",syndrome[r-j - 1])
				# print("syndrome:",syndrome)

				# pytanie czy dobrze jest indeksowanie
				#
				tmp = self.Berlekamp_Massey_multiplication(error_locator[l - j], syndrome[r - j - 1])
				# print(f"tmp: {tmp} in j = {j}")
				sum_poly.append(tmp)

		# print("sum_poly:", sum_poly)


		# sum everything
		s = 0
		res = [None]
		if len(sum_poly) > 1:	
			for i in range(len(sum_poly)):
				# print("---------------------------")
				# print(f"i = {i}")
				# print(f"sum[i] = {sum_poly[i]}")

				res[0] = self.XOR(self.field[res[0]], self.field[sum_poly[i]])
				# print(f"res = {res}")
				# print("---------------------------")

		else:
			res = sum_poly
		# print(f"res: {res[0]}")
		# print("syndrome[r-1]:", syndrom_poly[r-1])
		if r != 1:
			# problem kiedy syndrom_poly[iteration] = 0
			# wtedy 0 - 0 = 0
			# i kiedy 0 - None

			# co zrobić kiedy symdrom[iteration] = None ???
			if syndrom_poly[r-1] != None:
				# s = syndrom_poly[iteration] - s
				# print(syndrom_poly[r-1])
				# print(res[0])
				s = self.XOR(self.field[syndrom_poly[r-1]], self.field[res[0]])
				
			else:
				return 0
		return s


	def Berlekamp_Massey_multiplication(self, error_locator, syndrome):
		if error_locator == None or syndrome == None:
			return None
		if error_locator == 0:
			return syndrome
		if syndrome == 0:
			return error_locator
		res = (error_locator + syndrome) % self.n
		# res = (error_locator + syndrome) % 15

		return res

	def AdjustReminder(self, reminder_poly):
		length = len(reminder_poly)
		full_length = 2 * self.t
		rem = reminder_poly[:]
		if length != full_length:
			delta = full_length - length
			for i in range(delta):
				rem.insert(0, None)
		return rem


	################################
	#							   #
	#		{	 Forney		}	   #
	#							   #
	################################


	# syndromy też mają być odwrócone !!!
	# prawidłowe: ..., s^4, s^3, s^2, s^1
	#
	def Forney(self, error_locator, syndrome, chien_poly):
		syndrome.reverse()
		omega_poly = self.Omega(error_locator, syndrome)
		derivative_poly = self.Derivative(error_locator)
		error_values = []
		for i in range(len(chien_poly)):
			eval_omega = self.Omega_arg(omega_poly, chien_poly[i])
			eval_derivative = self.Derivative_arg(derivative_poly, chien_poly[i])
			div_res = self.Divide_GF(eval_omega, eval_derivative)
			inverse = self.Inverse_of_symbol(chien_poly[i])
			error = self.Multiply_GF(inverse, div_res)
			error_values.append(error)

		exponent_locations = self.Inverse_of_symbol_poly(chien_poly)
		if exponent_locations == []:
			return []
		max_exp = max(exponent_locations)
		# sort from max to min
		exponent_locations.sort(reverse=True)
		# print(f"exponent_locations = {exponent_locations}")

		error_polynomial = []
		# max_exp + 1 is equivalent to x^max_exp + ... + constant
		# for i in range(max_exp + 1):
			# error_polynomial.append(None)
		# print(f"error_polynomial = {error_polynomial}")

		for i in range(self.n):
			error_polynomial.append(None)

		positions = chien_poly[:]
		positions.sort()
		# print(f"positions = {positions}")

		for i in range(len(positions)):
			# print(f"i = {i}")
			# print(f"positions[{i}] - 1 = {positions[i] - 1}")
			# print(f"error_values[{i}] = {error_values[i]}")
			error_polynomial[positions[i] - 1] = error_values[i]

		# print(f"error = error_poly = {error_polynomial}")
		return error_polynomial



		# potrzebne?
	# def Adjust_forney_poly(self, error_values, exponent_locations):
	# 	pass


	# error evaluator
	def Omega(self, error_locator_poly, syndrom_poly):
		product = self.Multiply_poly(error_locator_poly, syndrom_poly)
		mod_arg = self.Omega_mod_argument_generator()
		omega = self.Modulo_poly(product, mod_arg)
		return omega


	# omega poly value at x
	def Omega_arg(self, omega_poly, arg):
		# print(f"omega_poly = {omega_poly}")
		exp = len(omega_poly)	# -1 ???
		result_xor = None
		for i in range(exp):
			# print(f"power = {(exp - i)}")
			# print(f"poly = {omega_poly[i]}")
			calc_arg = arg * (exp - i)
			# print(f"calc_arg = {calc_arg}")

			# print(f"polynomial[{i}] = {omega_poly[i]}")
			result_mul = self.Multiply_GF(omega_poly[i], calc_arg)
			# print(f"result_mul = {result_mul}")

			result_xor = self.XOR(self.field[result_xor], self.field[result_mul])
			# print(f"result_xor = {result_xor}")
			# print("__________________________")
		# print(result_xor)
		return result_xor


	# removes last index
	# it is a definition from fortney algorithm
	def Derivative(self, polynomial):
		del polynomial[-1]
		return polynomial


	def Derivative_arg(self, polynomial, arg):
		poly_len = len(polynomial)
		result_xor = None
		exp = poly_len - 1

		# pętla powinna być do t
		for i in range(poly_len):
			# print(f"power = {(arg - i)}")
			# print(f"poly = {polynomial[i]}")
			# calc_arg = self.Multiply_GF(arg, exp - i)
			calc_arg = arg * (exp - i)

			# print(f"calc_arg = {calc_arg}")

			result_mul = self.Multiply_GF(polynomial[i], calc_arg)
			# print(f"result_mul = {result_mul}")
			result_xor = self.XOR(self.field[result_xor], self.field[result_mul])
			# print(f"result_xor = {result_xor}")
		# print(result_xor)
		return result_xor


	# multiply powers of symbols
	# eg. alfa^2 * alfa^7 = alfa^9
	#
	def Multiply_GF(self, a, b):
		if a == None or b == None:
			return None
		power = (a + b) % self.n
		return power

	# divide two symbols
	# eg. alfa^14 / alfa^2 = alfa^12
	def Divide_GF(self, a, b):
		if a == None or b == None:
			return None
		sub = (a - b) % self.n
		return sub


	# inverse of given symbol
	# eg. for GF(2^5) where n = 31
	#	  alfa_4^-1 = alfa_27
	def Inverse_of_symbol(self, symbol):
		res = abs((symbol - self.n)) 
		return res

	# inverses of all given symbols in polynomials
	def Inverse_of_symbol_poly(self, symbols):
		symbols = symbols[:]
		for i in range(len(symbols)):
			symbols[i] = self.Inverse_of_symbol(symbols[i])
		return symbols

	# generate 2t+1 polynomial for omega function
	def Omega_mod_argument_generator(self):
		# Ω(x) mod x^2t+1, but exp 0 is 1, so [0, 8 * None] -> size = 2t+1
		modulo_poly = [0]
		for i in range(2*self.t):
			modulo_poly.append(None)
		return modulo_poly


	#############################################
	#
	#				chien search
	#
	#############################################

	def XOR_res_decimal(self, a ,b):
		c = ''
		for x in range(len(a)):
			if a[x] != b[x]:
				c = c+'1'
			else:
				c = c+'0'
		return self.Return_key(c);

		# multiply powers of symbols
	# eg. alfa^2 * alfa^7 = alfa^9
	#
	def Multiply_GF(self, a, b):
		if a == None or b == None:
			return None
		power = (a + b) % self.n
		return power


	# chien search finds roots of error locator polynomial by simply substituting symbols as arguments of a error locator polynomial
	def chien_search(self, error_locator):
		roots = []

		# przypadek dla x = 0
		xor = None
		for i in range(len(error_locator)):
			xor = self.XOR_res_decimal(self.field[xor], self.field[error_locator[i]])
		# print(f"xor = {xor}")
		if xor == None:
			roots.append(0)

		# przypadek dla alfa^x dla x >= 1
		x = 1
		err_len = len(error_locator) - 1
		for i in range(self.n - 1):
			result_xor = None
			for j in range(err_len):
				# result_xor = self.XOR(self.field[error_locator[j]], self.field[x * (err_len - j)])
				# print(f"power = {(x * (err_len - j))}")
				# print(f"error_locator = {error_locator[j]}")

				result_mul = self.Multiply_GF(error_locator[j], (x * (err_len - j)))
				# print(f"result_mul = {result_mul}")
				result_xor = self.XOR_res_decimal(self.field[result_xor], self.field[result_mul])
				# print(f"result_xor = {result_xor}")
			result_xor = self.XOR_res_decimal(self.field[result_xor], self.field[0])
			# print(f"result_xor ending = {result_xor}")
			# print("_____________________________")
			if result_xor == None:
				roots.append(x)
			x += 1
		return roots


	############################################
	#
	#				complete decoder
	#
	############################################


	# wynik forney'a zxorować z wiadomością

	# napisać pełne dekodowanie wiadomości bitowej, czyli przekopiować z funkcji poniżej 

	def Complete_decoder(self, message_poly):
		# ????
		syndromes = self.Find_syndromes(message_poly)
		#
		#
		print(f"syndromes = {syndromes}\n")
		check = [None] * (2*self.t)
		if syndromes == [] or syndromes == check:
			print(f"message received without errors!\n")
			return message_poly
		# print("######################################")

		# print(f"syndromes = {syndromes}")
		error_locator = self.Berlekamp_Massey(syndromes)
		print(f"error_locator Berlekamp_Massey = {error_locator}\n")
		# print("######################################")

		error_locator_roots = self.chien_search(error_locator)
		print(f"error locator roots = {error_locator_roots}\n")
		# print("######################################")
		error_correction = self.Forney(error_locator, syndromes, error_locator_roots)
		# print(f"syndromes = {syndromes}")

		print(f"error correction polynomial = {error_correction}\n")

		corrected_message = self.XORonPolynomials(message_poly, error_correction)
		# xorowanie

		print(f"corrected_message = {corrected_message}\n")
		return corrected_message


	def Find_syndromes(self, codeword):
		syndromes = []
		for arg in range(1, 2*self.t + 1):
			# print(f"codeword: {codeword}")
			# print(f"arg: {arg}")

			result_xor = None
			exp = len(codeword)
			for coef in range(len(codeword)):
				# print(f"coef: {coef}")
				if coef != len(codeword) - 1:
					# sprawdzić czy w omega_arg oraz derivative_arg nie potrzeba odjąć od potęgi -1 tak jak tutaj
					multiplied_exp = arg * (exp - coef - 1)
					# print(f"multiplied_exp = {multiplied_exp}")
					result_mul = self.Multiply_GF(multiplied_exp, codeword[coef])
					# print(f"result_mul = {result_mul}")
					# print(f"codeword[coef] = {codeword[coef]}")

					# print(f"result_xor before = {result_xor}")
					# print(self.XOR(self.field[result_xor], self.field[codeword[coef]]))
					result_xor = self.XOR(self.field[result_xor], self.field[result_mul])
					# print(f"result_xor after = {result_xor}")
					# print("__________________________")
				else:
					# print(f"codeword[coef] = {codeword[coef]}")

					# print(f"result_xor before = {result_xor}")
					result_xor = self.XOR(self.field[result_xor], self.field[codeword[coef]])
					# print(f"result_xor after = {result_xor}")
					# print("__________________________")
			# print(f"result_xor at {arg} = {result_xor}")

			syndromes.append(result_xor)
			# print("########################################################################")
		
		# print(syndromes)

		return syndromes



	def DecodeMessage_complete_decoder(self, message):
		# self.GenerateGF()
		buffer = self.SplitMessage(message)
		# print(buffer)
		message_poly = self.GetMessagePoly(buffer)
		# print(message_poly)
		res = list()
		for poly in message_poly:
			res.append(self.Complete_decoder(poly))
		bits = str()
		for i in res:
			bits += self.PolyToBits(i)
		return bits



	def Test_complete_decoder(self, message_without_errors, error_number, loop_number):
		buffer = self.SplitMessage(message_without_errors)
		message_poly = self.GetMessagePoly(buffer)
		res = list()
		for poly in message_poly:
			not_corrected = 0
			for i in range(loop_number):
				msg = self.MakeErrors(poly, error_number)
				one_res = self.Complete_decoder(msg)
				bits = self.PolyToBits(one_res)
				print(bits)
				if bits != message_without_errors:
					not_corrected += 1
			print("Not corrected:", not_corrected)
			print("Corrected:", loop_number - not_corrected)
		return not_corrected




	def TestBerlekamp_Massey(self, message_poly):
		syndromes = self.Find_syndromes(message_poly)
		#
		print(f"syndromes = {syndromes}\n")
		check = [None] * (2*self.t)
		if syndromes == [] or syndromes == check:
			print(f"message received without errors!\n")
			return message_poly

		error_locator = self.Berlekamp_Massey(syndromes)
		print(f"error_locator Berlekamp_Massey = {error_locator}\n")


	def MessageToPoly_one_message(self, message):
		# self.GenerateGF()
		buffer = self.SplitMessage(message)
		# print(buffer)
		message_poly = self.GetMessagePoly(buffer)
		# print(message_poly)
		return message_poly[0]

		


if __name__ == "__main__":
	# message [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 24, 1, 10, 19, 9, 16, None, 8, 23, 6, 28, 29]
	# message = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111100001010001001101101011011000000110101111010101011001101"


	# PRZYKŁAD Z PDF

	# wektor błędów
	# e(x) = αx^14 + α^2x^9 + x^4
	# [1, None, None, None, None, 2, None, None, None, None, 0, None, None, None, None]

	# polu = [6, 1, 0, 2, 1, 1, 2, 3, 7, 9, 5, 9, None, 0, 10]
	# poly = [1, None, None, None, None, 2, None, None, None, None, 0, None, None, None, None]

	# dla m = 4, t = 4
	# send_message = [6, 1, 0, 2, 1, 1, 2, 3, 7, 9, 5, 9, None, 0, 10]
	# send_message = "1100 0010 0001 0100 0010 0010 0100 1000 1011 1010 0110 1010 0000 0001 0111"


	# received_message = [11, 1, 0, 2, 1, 5, 2, 3, 7, 9, 10, 9, None, 0, 10]
	# received_message = "1110 0010 0001 0100 0010 0110 0100 1000 1011 1010 0111 1010 0000 0001 0111"



	# send_message = "110000100001010000100010010010001011101001101010000000010111"

	# received_message = "110000100001010000100010010010001011101001101010000000010111"



	# check = "000000000000000000001111000101111011101111100011000111011000"

	# received_message = "000000000000000000001111000101111011101111100011000111011000"




	# [None, None, None, None, None, 12, 0, 10, 7, 7, 11, 4, 0, 13, 3]
	# received_message = "0000 0000 0000 0000 0000 1111 0001 0111 1011 1011 1110 0011 0001 1101 1000"

	# [None, None, None, None, None, 12, 0, 8, 7, 7, 11, 4, 0, 13, 3]
	#										x
	# received_message = "0000 0000 0000 0000 0000 1111 0001 0111 1011 1011 1111 0011 0001 1101 1000"


	# send_message = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111100001010001001101101011011000000110101111010101011001001"

	# received_message = "00000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000111100001010001001101101011011000000110101111010101011001001"


	m = 5
	t = 5

	dec = Decoder(m, t)

	# poly = dec.MessageToPoly_one_message(received_message)
	# print(poly)

	# p = [12, None, 1, 14, 13, 13, 7, 2, 13, None, 12]

	# print(dec.TestBerlekamp_Massey(poly))

	# dec.Find_syndromes(received_message)


	
	# res = dec.DecodeMessage_complete_decoder(received_message)
	# print(res)
	# if send_message == res:
	# 	print("\noba są takie same")
	# else:
	# 	print("\nnie są")
	

	# dec.Test_complete_decoder(send_message, 1, 100)


	msg = "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
	print(dec.MakeErrors_burst_bits(msg, 3))


	poly = [None] * 31
	print(dec.MakeErrors_burst(poly, 3))

#####
##### dodać funkcję sprawdzającą działanie dekoderów przy pomocy burst errors (może to co jacek ma będzie dobre)
#####






	# a = [3, None, 5]
	# b = [2]
	# print(dec.Multiply_poly(a, b))

	# buffer = dec.SplitMessage(message)
	# print(buffer)
	# message_poly = dec.GetMessagePoly(buffer)
	# print(message_poly)
	# rem = []
	# for i in message_poly:
		# rem = dec.GetReminderOfMessagePoly(i)
		# rem = dec.AdjustReminder(rem)
	# print(rem)


	################################################################################
	# nie wiadomo czy wielomian syndromów zaczyna od x0 czy od x1
	# przykład. Ω(x) mod x9 = (α12x3+1)(α3x8 +x7 +α4x6 +α6x5 +α3x4 +α7x3 +α9x2 +α6x)
	# a6x jest ostatnie, a nie ma wolnego wyrazu, np. a2
	# w tym przypadku:
	# polynomial = [None,6,9,7,3,6,4,0,3]

	# polynomial = [6,9,7,3,6,4,0,3]
	# print(polynomial)
	# reverse, bo berlekamp_massey przyjmuje syndromy zaczynające się od najmniejszego indeksu
	# polynomial.reverse()

	# print(dec.Derivative(polynomial))
	# symbol = 6
	# print(f"inverse of {symbol} = {dec.Inverse_of_symbol(symbol)}")


	# result from chien search
	# Error positions = [1, 6, 11]

	#################################
	#
	#		forney tests
	#
	#################################

	# error_locator = [12, None, None, 0]
	# chien_search_poly = [1, 6, 11]
	# syndromes = polynomial

	# forney = dec.Forney(error_locator, syndromes, chien_search_poly)
	# print(f"forney = {forney}")



	# omega = dec.Omega(error_locator, syndromes)
	# print(f"omega = {omega}")
	# s = 3
	# print(f"Omega_arg for s = {s}: {dec.Omega_arg(omega, s)}")
	# print(f"multiply {first_poly}\nby {second_poly}")

	# product = dec.Multiply_poly(first_poly, second_poly)
	# print(f"omega mod 9 = {product}")

	# modulo_arg = dec.Omega_mod_argument_generator()
	# print(f"modulo = {modulo_arg}")
	# final_product = dec.Modulo_poly(product, modulo_arg)
	# print(f"final_product = {final_product}")

	# der = dec.Derivative(first_poly)
	# print(f"derivative of first_poly: {der}")
	# s = 2
	# print(f"derivative at s = {s}: {dec.Derivative_arg(der, s)}")

	# print(dec.XOR(dec.field[1], dec.field[6]))


	##################################

	# polynomial.reverse()

	# res = dec.Berlekamp_Massey(rem)
	# print(res)

	# rem.reverse()
	# print(rem)

















	# error_locator_poly = [None] * len(reminder)
	# error_locator_poly[0] = 1
	# dec.DiscrepancyOne(error_locator_poly, reminder, 3)



	# res = dec.TestCorrectionEffectiveness(message, 1, 1000)

	# poly = [1,2,3,4,5]
	# p = dec.MoveITimesLeft(poly, 1)
	# print(p)

	# for i in range(10):
	# 	res = dec.MakeErrors(err_message, 1)
	# 	print(res)

	# s = dec.BinaryStringToInt("11010")
	# print(s)

	# poly = [1] * 31
	# res = dec.MakeErrors(poly, 1)

	# print(res)

	####################################################################################
	#
	#	UWAGI!
	#
	#	problem z korekcją 25 pierwych bitów wiadomości
	#	z liczenia od i = 0 do k nie poprawia 5 pierwszych syndromów (25 bitów). 
	#	zamiana z liczenia od i = 0 do n poprawia 1 błąd w całej wiadomości
	#
	#	za duży rozrzut błędów od siebie sprawia, że dekoder nie poprawia błędów
	#
	####################################################################################




	# Zakres testów - dla  każdego testowanego dekodera należy:

    # generować losowe błędy w zakresie od jednego do t błędów (dla kodu BCH od 1 do t bitów błędnych; dla  RS od 1 do t błędnych symboli)
    # liczba generowanych błędów dla każdej krotności powinna zależeć od liczby możliwych błędów, tj. 
    #     dla błędów  pojedynczych i ewentualnie podwójnych można zrobić przegląd zupełny - proszę oszacować ile jest takich błędów  ile to będzie trwało i na tej podstawie zadecydować
    #     dla błędów pozostałej krotności należy wygenerować wiele losowych błędów  i sprawdzić czy są korygowane - liczba generowanych błędów zależy od tego czy wyniki jakie będziecie uzyskiwać zgadzają się z oczekiwaniami (analizą teoretyczną) i należy ich zrobić wystarczająco dużo aby taką analizę potwierdzić, albo obalić (wówczas pytanie czy na pewno macie dobrze zrealizowana implementacje albo czy dobrze przeprowadziliście analizę)
    # zrobić próby korekcji dla błędów o ekstremalnie dużej krotności tzn. większość symboli/bitów przekłamana. Testów  takich zrobic kilka tysięcy szukając przypadku w  którym dekoder powie, że naprawił błąd
    # testując  proszę zapisywać wyniki, które nie zgadzają się z teorią:
    #     jeśli dekoder nie koryguje błędu pojedynczego to należy go zapisać i przeanalizować dlaczego
    #     jeśli koryguje niektóre błędy określonej krotności a innych nie, to też należy je zapisać i znaleźć odpowiedź na pytanie dlaczego
    #     jeśli koryguje błędy ekstremalnie dużej krotności, to trzeba ustalić co się dzieje i dlaczego
    # wyniki testów należy przedstawić w postaci tabeli (krotność błędu, liczba testów, skuteczność korekcji (poprawne korekcje) [%], procent błędów skorygowanych niepoprawnie) oraz wykresów
    # uzyskane wyniki trzeba omówić - dlaczego tak wyszło, jak to uzasadniacie itp.