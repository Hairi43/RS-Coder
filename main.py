# Main file for RS error correction project

import galoisfield as gf
import decoder as decoder

class RS:

	def __init__(self, m, t):
		self.m = m
		self.t = t
		# galoa field (dictionary of alfas)
		self.field = dict()
		self.n = pow(2,m) - 1
		self.k = self.n - 2*t
		# self.char_dict = dict() # added alphabet 

	# create galois field (dictionary of alfas)
	def GenerateGF(self):
		g = gf.GaloisField(self.m, self.t)
		self.field = g.GenerateGF()
		self.char_dict = g.GetCharDict()

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
	# argument string
	# by default this method is used instead of 'SplitMessageToMBits()'
	def SplitMessage(self, message):
		message = message.replace(" ", "")
		if len(message) <= (self.k * self.m):
			splited_message = []
			splited_message.append(message[:(self.k * self.m)])
			splited_message[0] = self.SplitMessageToMBits(splited_message[0])
		else:
			splited_message = []
			for i in range(len(message)):
				if i % (self.k * self.m) == 0:
					splited_message.append(message[i:i+(self.k * self.m)])
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
			for j in i:
				str_tmp = str_tmp + self.field[j]
			str_tmp = self.AdjustMessage(str_tmp)
			list_tmp.append(str_tmp)
			str_tmp = ''
		for i in list_tmp:
			string += i
		return string

	# generate generating polynomial (x+a^1)(x+a^2)...(x+a^2t)
	def Generate_poly(self):
		n = []
		for i in range(self.t * 2):
			n.append([0, i+1])
		res = n[0]
		for i in range(self.t * 2 - 1):
			res = self.Multiply_poly(res, n[i+1])
		return res

	# polynomials multiplication
	# a = [4,None, 3] -> 4x^2 + 0x + 3
	def Multiply_poly(self,a,b):
		# self.Adjust_poly(a)
		# self.Adjust_poly(b)
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

	# polynomial division
	# a = [4,None, 3] -> 4x^2 + 0x + 3
	def Divide_poly(self, dividend, divisor):
		if divisor == [None]:
			# division by 0 == exception?
			return [None]

		if len(dividend) < len(divisor):
			return [None]

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
					dividend[index] = self.XOR(field[dividend[index]], field[None])
				else:
					coef = (mult + divisor[index]) % (pow(2,self.m) - 1)
					dividend[index] = self.XOR(field[dividend[index]], field[coef])
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

		print("dividend: ", dividend)
		print("divisor: ", divisor)

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
		return dividend

	# calculate reminder of message polynomial
	# return [[11, 16, 8, ...], [0, 4, None, ...]] or [[0, 4, ...]]
	# depending on message length
	def GetReminderOfMessagePoly(self, message_poly):
		gen_poly = self.Generate_poly()
		offset_poly = [0] + [None] * (2*self.t)
		tmp_poly = []
		for i in message_poly:
			tmp_poly.append(self.Multiply_poly(i, offset_poly))
		rem_poly = []
		for i in tmp_poly:
			rem_poly.append(self.Modulo_poly(i, gen_poly))
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


	def Fill_with_zeros(self, message_in_bits):
		full_length = self.m * self.n
		message_lenght = len(message_in_bits)
		if message_lenght != full_length:
			delta = full_length - message_lenght
			full_message = "0" * delta + message_in_bits
			return full_message
		return message_in_bits

	# code message (make full package with correction)
	# return string of bits of size m*n ready to send to decoder
	def CodeMessage(self, message):
		self.GenerateGF()
		buffer = self.SplitMessage(message)
		print(buffer)
		message_poly = self.GetMessagePoly(buffer)
		reminder_poly = self.GetReminderOfMessagePoly(message_poly)
		coded_msg = []
		for i in range(len(message_poly)):
			coded_msg.append(message_poly[i] + reminder_poly[i])
		#
		print("coded msg:", coded_msg)
		#
		message_in_bits = self.PolyToBits(coded_msg)
		message_in_bits = self.Fill_with_zeros(message_in_bits)
		return message_in_bits


if __name__ == "__main__":

	# message = "1001001101010110110100001 1001001101010110110100001 1001001101010110110100001 1001001101010110110100001 10010 10100 10010 10010"

	#			    4     29		3
	# message = "10000 01001 01000"

	#				   24    1
	message = "11110 00010"
	# message = "1111 0001"


	#				 29
	# message = "01001"

	# message = "1"
	m = 4
	t = 4

	# full codeword send to decoder
	# m = 4
	# t = 4
	# c(x) = α6x14 + αx13 + x12 + α2x11 + αx10 + αx9 +α2x8 + α3x7 + α7x6 + α9x5 + α5x4 + α9x3 + x + α10
	# [6, 1, 0, 2, 1, 1, 9, 2, 3, 7, 9, 5, 9, 0, 10]
	# message = "1100 0010 0001 0100 0010 0010 1010 0100 1000 1011 1010 0110 1010 0001 0111"


	rs = RS(m, t)

	#
	#	[  encode message  ]
	#
	
	message_in_bits = rs.CodeMessage(message)

	print(message_in_bits, "\n")


	# res = rs.Multiply_poly([8], [9, None])
	# print(res)
	

	#	[	decode message  ]
	#



	# dec = decoder.Decoder(m, t)

	# print("decoded_msg_reminder:", decoded_msg)


	#
	#	można zrobić dowolny dekoder
	#

