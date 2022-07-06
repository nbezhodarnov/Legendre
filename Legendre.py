import numpy as np
import matplotlib.pyplot as plt
import math
	
class Legendre():
	x = np.array([], dtype = float)
	y = np.array([], dtype = float)
	n = 0
	c = np.array([], dtype = float)
	c_size = 0
	
	def __init__(self, x_array_input, y_array_input, size):
		self.x = x_array_input
		self.y = y_array_input
		self.n = self.x.size
		self.c_size = size
		self.c = np.zeros(size, dtype = float)
		Gram_matrix = np.zeros((size, size), dtype = float)
		free_member_vector = np.zeros(size, dtype = float)
		for i in range(self.c_size):
			for k in range(self.n):
				Gram_matrix[i][i] += self.__Legendre_formula(i, self.x[k]) * self.__Legendre_formula(i, self.x[k])
			for j in range(i + 1, self.c_size):
				for k in range(self.n):
					Gram_matrix[i][j] += self.__Legendre_formula(i, self.x[k]) * self.__Legendre_formula(j, self.x[k])
				Gram_matrix[j][i] = Gram_matrix[i][j]
			for l in range(self.n):
				free_member_vector[i] += self.__Legendre_formula(i, self.x[l]) * self.y[l]
		print(Gram_matrix)
		for i in range(self.c_size - 1):
			for j in range(i + 1, self.c_size):
				if (Gram_matrix[j][i] != 0):
					for k in range(i + 1, self.c_size):
						Gram_matrix[j][k] -= Gram_matrix[i][k] * Gram_matrix[j][i] / Gram_matrix[i][i]
					free_member_vector[j] -= free_member_vector[i] * Gram_matrix[j][i] / Gram_matrix[i][i]
			for k in range(i + 1, self.c_size):
				Gram_matrix[k][i] = 0
		temp = self.c_size - 1
		self.c[temp] = free_member_vector[temp] / Gram_matrix[temp][temp]
		for i in range(1, self.c_size + 1):
			index = self.c_size - i
			temp = free_member_vector[index]
			for j in range(index + 1, self.c_size):
				temp -= Gram_matrix[index][j] * self.c[j]
			self.c[index] = temp / Gram_matrix[index][index]
			
		
	def __Legendre_formula(self, k, x_input):
		if (k == 0):
			return 1
		if (k == 1):
			return 2 * (x_input - self.x[0]) / (self.x[self.n - 1] - self.x[0]) - 1
		return ((2 * k - 1) * (2 * (x_input - self.x[0]) / (self.x[self.n - 1] - self.x[0]) - 1) * self.__Legendre_formula(k - 1, x_input) + (k - 1) * self.__Legendre_formula(k - 2, x_input)) / k
		
	def Legendre_calculate(self, x_input):
		if (self.c_size == 0):
			return 1
		if (self.c_size == 1):
			return self.c[0]
		if (self.c_size == 2):
			return self.c[0] + self.c[1] * (2 * (x_input - self.x[0]) / (self.x[self.n - 1] - self.x[0]) - 1)
		x_modified = 2 * (x_input - self.x[0]) / (self.x[self.n - 1] - self.x[0]) - 1
		pre_last = 1
		last = x_modified
		result = self.c[0] + self.c[1] * x_modified
		temp = 0
		for i in range(2, self.c_size):
			temp = ((2 * i - 1) * x_modified * last + (i - 1) * pre_last) / i
			pre_last = last
			last = temp
			result += self.c[i] * temp
		return result
	
def main():
	x = np.array([0, 1.75, 3.5, 5.25, 7])
	y = np.array([0, -1.307, -2.211, -0.927, -0.871])
	polynomial = Legendre(x, y, 4)
	x_plot_table = np.linspace(0, 7, 50, dtype = float)
	y_plot_table = np.linspace(0, 0, 50, dtype = float)
	y_plot_table_of_original = np.linspace(0, 0, 50, dtype = float)
	x = np.append(x, 2.555)
	y = np.append(y, polynomial.Legendre_calculate(2.555))
	for i in range(50):
		y_plot_table[i] = polynomial.Legendre_calculate(x_plot_table[i])
		y_plot_table_of_original[i] = math.cos(x_plot_table[i]) - 2 ** (0.1 * x_plot_table[i])
	fig, ax = plt.subplots()
	plt.plot(x_plot_table, y_plot_table, 'b-', label = 'Legendre')
	plt.plot(x_plot_table, y_plot_table_of_original, 'm--', label = 'Original')
	plt.plot(x, y, 'r*')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.annotate('f(x) ~ ' + str(y[y.size - 1]), xy=(x[x.size - 1], y[y.size - 1]), xytext=(x[x.size - 1], y[y.size - 1] - 0.32),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )
	plt.legend(loc='upper right')
	plt.show()
	
if __name__ == '__main__':
    main()
