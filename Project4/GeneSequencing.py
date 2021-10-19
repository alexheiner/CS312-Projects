#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


# !!for back pointers!!
## diagonal(match) = 0, diagonal(sub) = 1, insert(left) = 2, delete(top) = 3

class GeneSequencing:

	def __init__( self ):
		pass

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.align_length = align_length

		if banded:
			score, al1, al2 = self.rest_dist(seq1[:align_length], seq2[:align_length])
		else:
			score, al1, al2 = self.norest_dist(seq1[:align_length], seq2[:align_length])

		alignment1 = '{}'.format(al1[:100])
		alignment2 = '{}'.format(al2[:100])

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	# time and space O(km)
	def rest_dist(self, seq1, seq2):
		if len(seq1) > len(seq2) and len(seq1) - len(seq2) > 7:
			return math.inf, "No Alignment Possible", "No Alignment Possible"
		elif len(seq2) > len(seq1) and len(seq2) - len(seq1) > 7:
			return math.inf, "No Alignment Possible", "No Alignment Possible"
		self.matrix_band = [[math.inf for x in range(7)] for y in range(len(seq1) + 1)]
		self.pointers_band = [[math.inf for x in range(7)] for y in range(len(seq1) + 1)]
		self.pointers_band[0][0] = -math.inf
		offset = 1
		d = 3
		for i in range(len(seq1) + 1):
			if i < d+1:
				for j in range(d + i + 1):
					# compute path from 0 to d + 2 normally
					if i == 0:
						if j == 0:
							self.matrix_band[i][j] = 0
						else:
							# populate base cases
							self.matrix_band[i][j] = 5 * j
					else:
						if j == 0:
							# populate base cases
							self.matrix_band[i][j] = 5 * i
						else:
							# get characters from seq1 and 2
							ch1 = seq1[i - 1]
							ch2 = seq2[j - 1]
							# if characters don't match get the minimum path
							if ch1 != ch2:
								if j == d + i + 1:
									# cant check up
									self.matrix_band[i][j] = self.get_min_res(i, j, True, False, False)
								else:
									# can check up,left, and diagonally
									self.matrix_band[i][j] = self.get_min_res(i, j, True, True, False)
							else:
								# if they do match we need to subtract three from the diagonal
								#self.matrix_band[i][j] = (self.matrix_band[i - 1][j - 1]) - 3
								#self.pointers_band[i][j] = 0
								self.matrix_band[i][j] = self.get_min_res(i, j, True, True, True)
			else:
				self.offsetEngaged = True
				spotsFilled = 0
				# keep track of the offset so we know where we are in sequence 2
				for j in range(offset, len(seq2) + 1):
					# if we reach our limit, break
					if spotsFilled == 7:
						break
					else:
						ch1 = seq1[i - 1]
						ch2 = seq2[j - 1]
						if ch1 != ch2:
							if j == offset:
								# we can only have diagonal or up
								self.matrix_band[i][j-offset] = self.get_min_shifted(i, j-offset, False, True, False)
							elif spotsFilled == 6:
								# we are at the end, can only have diagonal or left
								self.matrix_band[i][j - offset] = self.get_min_shifted(i, j-offset, True, False, False)
							else:
								self.matrix_band[i][j - offset] = self.get_min_shifted(i, j-offset, True, True, False)
						else:
							# match
							#self.matrix_band[i][j-offset] = self.matrix_band[i - 1][j-offset] - 3
							#self.pointers_band[i][j-offset] = 0
							if spotsFilled == 6:
								self.matrix_band[i][j - offset] = self.get_min_shifted(i, j-offset, True, False, True)
							else:
								self.matrix_band[i][j - offset] = self.get_min_shifted(i, j-offset, True, True, True)


					spotsFilled += 1
				offset += 1
		i = len(seq1)
		j = 6
		num = self.matrix_band[i][j]
		while num == math.inf:
			j -= 1
			num = self.matrix_band[i][j]
		mod_seq1, mod_seq2 = self.trace_band_path(seq1, j, seq2)
		return num, mod_seq1, mod_seq2

	# time O(n) space O(n)
	def trace_band_path(self, seq1, j, seq2):
		d = 3
		i = len(seq1)
		mod_seq1 = seq1[:self.align_length]
		mod_seq2 = seq2[:self.align_length]
		num = self.pointers_band[i][j]
		while num != -math.inf:

			# match from diagonal or substitute from diagonal
			if num == 0:
				# if our row is after we have to keep track of an offset diagonal is different
				if i > d:
					i -= 1
					num = self.pointers_band[i][j]
				else:
					i -= 1
					j -= 1
					num = self.pointers_band[i][j]
			# insert from left
			elif num == 1:
				j -= 1
				num = self.pointers_band[i][j]
				# insert dash into sequence 1 string
				mod_seq1 = mod_seq1[:i] + "-" + mod_seq1[i:]
			# delete from top
			elif num == 2:
				# if our row is after we have to keep track of an offset up is different
				if i > d:
					i -= 1
					j += 1
					num = self.pointers_band[i][j]
				else:
					i -= 1
					num = self.pointers_band[i][j]
				if self.offsetEngaged and i > 4:
					mod_seq2 = mod_seq2[:i - d + j] + "-" + mod_seq2[i - d + j:]
				else:
					mod_seq2 = mod_seq2[:i - 1] + "-" + mod_seq2[i - 1:]

			else:
				num = -math.inf
		return mod_seq1, mod_seq2


	# time O(1) space O(1)
	def get_min_res(self, i, j, canLeft, canTop, isMatch):
		# get min for banded algorithm with restrictions on top and left- want to stay in bounds
		# diagonal and up are still the same, just needed a function to implement restrictions
		if canLeft:
			left = (self.matrix_band[i][j - 1]) + 5
		else:
			left = math.inf
		if canTop:
			top = (self.matrix_band[i - 1][j]) + 5
		else:
			top = math.inf
		if isMatch:
			diag = self.matrix_band[i - 1][j - 1] - 3
		else:
			diag = (self.matrix_band[i - 1][j - 1]) + 1

		min_num = min(diag, top, left)

		# order to return for tie break
		if left == min_num:
			self.pointers_band[i][j] = 1
			return left
		elif top == min_num:
			self.pointers_band[i][j] = 2
			return top
		else:
			self.pointers_band[i][j] = 0
			return diag

	# time O(1) space O(1)
	def get_min_shifted(self, i, j, canLeft, canTop, isMatch):
		# get min for banded algorithm with restrictions on top and left- want to stay in bounds
		# diagonal is now above, left is still same, top is diagonal to the right
		if canLeft:
			left = (self.matrix_band[i][j - 1]) + 5
		else:
			left = math.inf
		if canTop:
			top = (self.matrix_band[i - 1][j + 1]) + 5
		else:
			top = math.inf
		if isMatch:
			diag = self.matrix_band[i - 1][j] - 3
		else:
			diag = (self.matrix_band[i - 1][j]) + 1

		min_num = min(diag, top, left)

		# order to return for tie break
		if left == min_num:
			self.pointers_band[i][j] = 1
			return left
		elif top == min_num:
			self.pointers_band[i][j] = 2
			return top
		else:
			self.pointers_band[i][j] = 0
			return diag


	# Time and space O(mn)
	def norest_dist(self, seq1, seq2):
		# init first row and column
		self.matrix = [[math.inf for x in range(len(seq2) + 1)] for y in range(len(seq1) + 1)]
		self.back_pointers = [[math.inf for x in range(len(seq2) + 1)] for y in range(len(seq1) + 1)]

		self.back_pointers[0][0] = -math.inf
		self.matrix[0][0] = 0

		# set up base cases for sequence 2
		for i in range(len(seq2) + 1):
			if i == 0:
				continue
			else:
				self.matrix[0][i] = INDEL * i

		# set up base cases for sequence 1
		for i in range(len(seq1) + 1):
			if i == 0:
				continue
			else:
				self.matrix[i][0] = INDEL * i

		for i in range(1, len(seq1) + 1):
			ch1 = seq1[i - 1]
			for j in range(1, len(seq2) + 1):
				ch2 = seq2[j - 1]
				# if characters don't match get the minimum path
				if ch1 != ch2:
					self.matrix[i][j] = self.get_min(i, j, False)
				else:
					# if they do match we need to subtract three from the diagonal
					self.matrix[i][j] = self.get_min(i, j, True)

		mod_seq1, mod_seq2 = self.trace_path(seq1, seq2)
		return self.matrix[len(seq1)][len(seq2)], mod_seq1, mod_seq2


	# time O(n) space O(1)
	def trace_path(self, seq1, seq2):
		i = len(seq1)
		j = len(seq2)
		mod_seq1 = seq1
		mod_seq2 = seq2
		num = self.back_pointers[i][j]
		while num != -math.inf:
			# match from diagonal
			if num == 0:
				i -= 1
				j -= 1
				num = self.back_pointers[i][j]
			# insert from left
			elif num == 1:
				j -= 1
				num = self.back_pointers[i][j]
				# insert dash into sequence 1 string
				mod_seq1 = mod_seq1[:i] + "-" + mod_seq1[i:]
			# delete from top
			elif num == 2:
				i -= 1
				num = self.back_pointers[i][j]
				mod_seq2 = mod_seq2[:j] + "-" + mod_seq2[j:]
			else:
				num = -math.inf
		return mod_seq1, mod_seq2

	# time and space O(1)
	def get_min(self, i, j, isMatch):
		top = (self.matrix[i-1][j]) + 5
		left = (self.matrix[i][j-1]) + 5
		if isMatch:
			diag = self.matrix[i - 1][j - 1] - 3
		else:
			diag = (self.matrix[i - 1][j - 1]) + 1
		min_num = min(diag, top, left)

		# order to return for tie break
		if left == min_num:
			self.back_pointers[i][j] = 1
			return left
		elif top == min_num:
			self.back_pointers[i][j] = 2
			return top
		else:
			self.back_pointers[i][j] = 0
			return diag

