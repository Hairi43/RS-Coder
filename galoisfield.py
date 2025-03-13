class GaloisField():
	prim_poly = {
		1: [1,0],
		2: [2,1,0],
		3: [3,1,0],
		4: [4,1,0],
		5: [5,2,0],
		6: [6,1,0],
		7: [7,3,0],
		8: [8,4,3,2,0],
		9: [9,4,0],
		10:[10,3,0]
	}

	def __init__(self, m, t):
		self.m = m
		self.t = t
		self.n = (2 ** m) - 1
		self.k = self.n - (2*self.t)
		self.field = dict()
		

	# binary xor on two binary strings
	def XOR(self, a, b):
		c = ''
		for x in range(len(a)):
			if a[x] != b[x]:
				c = c+'1'
			else:
				c = c+'0'
		return c;

	def XOR_res_decimal(self, a ,b):
		c = ''
		for x in range(len(a)):
			if a[x] != b[x]:
				c = c+'1'
			else:
				c = c+'0'
		return self.Return_key(c);


	# get key from value's dictionary
	# eg. get 5 from 00101 in dict
	def Return_key(self, n):
		for key, value in self.field.items():
			if n == value:
				return key;
	#converts the list representation of polynomial to number
	def Poly_to_int(self, poly):
		value = 0
		for exp in poly:
			value |= 1 << exp # value OR 1 and shift by exp to the left
		return value

	# returns list that doesn't contain repeated elements in poly
	def OccurenceCounter(self, poly):
		tmp = poly[:]
		result = list()
		for i in range(len(poly)):
			counter = 0
			for j in range(len(tmp)):
				if poly[i] == tmp[j]:
					counter += 1
			if counter == 1:
				result.append(poly[i])
		return result

	# multiplies polynomial coefficients by alfa and returns
	# binary string that goes into galois field and
	# coef list that is the next argument
	def MultiplyByAlfa(self, m_poly):
		if m_poly == []:
			return [None]
		# "00000" -> None
		string = "0"*self.m
		for index in range(len(m_poly)):
			m_poly[index] += 1
		for i in range(len(m_poly)):
			if m_poly[i] == self.m:
				m_poly.pop(i)
				m_poly = m_poly + self.prim_poly[self.m][1:]
		unique_list = self.OccurenceCounter(m_poly)
		for index in range(len(unique_list)):
			string = self.XOR(string, self.field[unique_list[index]])
		return string, unique_list



	# generate galois field from class's M value
	def GenerateGF(self):
		self.field[None] = '0'*self.m
		# [0, m-1]
		for index in range(self.m):
			self.field[index] = '0'*self.m 	# 00000
			tmp = self.field[index]
			tmp = list(tmp)
			tmp[self.m-1-index] = '1'
			tmp = ''.join(tmp)				# 00001 -> 00010 ...
			self.field[index] = tmp
		alfa_m = self.prim_poly[self.m][1:]
		self.field[self.m] = '0'*self.m
		# [m]
		for value in alfa_m:
			self.field[self.m] = self.XOR(self.field[self.m], self.field[value])
		# [m+1, 2^m -2]
		m_poly = self.prim_poly[self.m][1:]
		for index in range(self.m + 1, pow(2,self.m) - 1):
			self.field[index], m_poly = self.MultiplyByAlfa(m_poly)
		return self.field

	# print galois field as a dictionary
	def ShowGF(self):
		print(self.field)

	def GetCharDict(self):
		char_dict = {
   			None: "00000",
    		"a": "00001",
    		"b": "00010",
    		"c": "00100",
    		"d": "01000",
    		"e": "10000",
    		"f": "00101",
    		"g": "01010",
    		"h": "10100",
    		"i": "01101",
    		"j": "11010",
    		"k": "10001",
    		"l": "00111",
    		"m": "01110",
    		"n": "11100",
    		"o": "11101",
    		"p": "11111",
    		"r": "11011",
    		"s": "10011",
    		"t": "00011",
    		"u": "00110",
    		"v": "01100",
    		"w": "11000",
    		"x": "10101",
    		"y": "01111",
    		"z": "11110",
    		".": "11001",
    		",": "10111",
    		"?": "01011",
    		"!": "10110",
    		" ": "01001",
    		"'": "10010"
		}
		return char_dict


	########################################
	#									   #
	#		{	chien search 	}		   #
	#									   #
	########################################

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
		print(f"xor = {xor}")
		if xor == None:
			roots.append(0)

		# przypadek dla alfa^x dla x >= 1
		x = 1
		err_len = len(error_locator) - 1
		for i in range(self.n - 1):
			result_xor = None
			for j in range(err_len):
				# result_xor = self.XOR(self.field[error_locator[j]], self.field[x * (err_len - j)])
				print(f"power = {(x * (err_len - j))}")
				print(f"error_locator = {error_locator[j]}")

				result_mul = self.Multiply_GF(error_locator[j], (x * (err_len - j)))
				print(f"result_mul = {result_mul}")
				result_xor = self.XOR_res_decimal(self.field[result_xor], self.field[result_mul])
				print(f"result_xor = {result_xor}")
			result_xor = self.XOR_res_decimal(self.field[result_xor], self.field[0])
			print(f"result_xor ending = {result_xor}")
			print("_____________________________")
			if result_xor == None:
				roots.append(x)
			x += 1
		return roots


if __name__ == "__main__":
	# gf = GaloisField(5, 5)
	gf = GaloisField(4, 4)

	# gf.GenerateGF()
	# gf.ShowGF()
	# error_locator = [1, 2, 3]  # Example coefficients for Λ(x) = 1 + 2x + 3x^2

	# przykład z dekodera
	# error_locator = [18, 20, 9, 10, 1, 0]

	# error_locator = [12, None, None, 0]


	# roots = gf.chien_search(error_locator) 
	# print(f"Error positions: {roots}") #Pozycje błędów