import sympy


class ODE:

    moultonCoeff = [[1, 0, 0, 0, 0, 0, 0, 0],
                    [1 / 2, 1 / 2, 0, 0, 0, 0, 0, 0],
                    [5 / 12, 2 / 3, -1 / 12, 0, 0, 0, 0, 0],
                    [3 / 8, 19 / 24, -5 / 24, 1 / 24, 0, 0, 0, 0],
                    [251 / 720, 323 / 360, -11 / 30, 53 / 360, -19 / 720, 0, 0, 0],
                    [95 / 288, 1427 / 1440, -133 / 240, 241 / 720, -173 / 1440, 3 / 160, 0, 0],
                    [19087 / 60480, 2713 / 2520, -15487 / 20160, 586 / 945, -6737 / 20160, 263 / 2520, -863 / 60480, 0],
                    [5257 / 17280, 139849 / 120960, -4511 / 4480, 123133 / 120960, -88547 / 120960, 1537 / 4480, -11351 / 120960, 275 / 24192]]

    bashforthCoeff = [[1, 0, 0, 0, 0, 0, 0, 0],
                      [3 / 2, -1 / 2, 0, 0, 0, 0, 0, 0],
                      [23 / 12, -4 / 3, 5 / 12, 0, 0, 0, 0, 0],
                      [55 / 24, -59 / 24, 37 / 24, -3 / 8, 0, 0, 0, 0],
                      [1901 / 720, -1387 / 360, 109 / 30, -637 / 360, 251 / 720, 0, 0, 0],
                      [4277 / 1440, -2641 / 480, 4991 / 720, -3649 / 720, 959 / 480, -95 / 288, 0, 0],
                      [198721 / 60480, -18637 / 2520, 235183 / 20160, -10754 / 945, 135713 / 20160, -5603 / 2520, 19087 / 60480, 0],
                      [16083 / 4480, -1152169 / 120960, 242653 / 13440, -296053 / 13440, 2102243 / 120960, -115747 / 13440, 32863 / 13440, -5257 / 17280]]

    def __init__(self, y, t, h, func, order=1):
        self.y = y
        self.t = t
        self.h = h
        self.func = func
        self.order = order

    def euler(self):
    	return self.y + self.h * self.func(self.t, self.y)

    def backward_euler(self):
	    y1 = euler(self.y, self.t, self.h, self.steps, self.func)
	    return self.y + self.h * self.func(self.t + self.h, y1)

    def modified_euler(self):
	    y1 = euler(self.y, self.t, self.h, self.steps, self.func)
	    return self.y + (self.h / 2) * (self.func(self.t, self.y) + self.func(self.t + self.h, y1))

    def runge_kutta(self):
	    k1 = self.func(self.t, self.y)
	    k2 = self.func(self.t + self.h / 2, self.y + k1 * self.h / 2)
	    k3 = self.func(self.t + self.h / 2, self.y + k2 * self.h / 2)
	    k4 = k2 = self.func(self.t + self.h, self.y + k3 * self.h)
	    return self.y + (k1 + 2 * k2 + 2 * k3 + k4) * self.h / 6

    def adam_bashforth(self):
        return self.bashforthCoeff[self.order][i] * self.func(self.t, self.y) * self.h

    def adam_multon(self):
        return self.moultonCoeff[self.order][i] * self.func(self.t, self.y) * self.h


class Solver:

	x_axis = []
	y_axis = []

	def __init__(self, method, steps, ode):
		self.method = method
		self.steps = steps
		self.ode = ode

	def solve(self)
		if self.method == "euler":
			pass
		elif self.method == "euler_inverso":
			pass
		elif self.method == "euler_aprimorado":
			pass
		elif self.method == "runge_kutta":
			pass
		elif self.method == "adam_bashforth":
			pass
		elif self.method == "adam_bashforth_by_euler":
			pass
		elif self.method == "adam_bashforth_by_euler_inverso":
			pass
		elif self.method == "adam_bashforth_by_euler_aprimorado":
			pass
		elif self.method == "adam_bashforth_by_runge_kutta":
			pass
		elif self.method == "adam_multon":
			pass
		elif self.method == "adam_multon_by_euler":
			pass
		elif self.method == "adam_multon_by_euler_inverso":
			pass
		elif self.method == "adam_multon_by_euler_aprimorado":
			pass
		elif self.method == "adam_multon_by_runge_kutta":
			pass
		elif self.method == "formula_inversa":
			pass
		elif self.method == "formula_inversa_by_euler":
			pass
		elif self.method == "formula_inversa_by_euler_inverso":
			pass
		elif self.method == "formula_inversa_by_euler_aprimorado":
			pass
		elif self.method == "formula_inversa_by_runge_kutta":
			pass

		def euler(self):
			points = []
			points.append(ode.t, ode.y)
			for i in range(steps):
				ode.y = ode.euler
				points.append(ode.t, )
				ode.t += h
				ode.y =


def main():

	# Defining equation symbols
	t, y = sympy.symbols("t y")

    # Reading input file
    inFile = open("entrada.txt", "r")
    outFile = open("saida.txt", "w")
    lines = read.readlines()

    for line in lines:
        words = line.split()

        method = words[0]
        y0 = float(words[1])
        t0 = float(words[2])
        h = float(words[3])
        steps = int(words[4])
        # Converts the function string into a lambda function
        f = sympy.lambdify([t,y], sympy.sympify(words[5]), "math")
        if len(words) > 6:
        	order = int(words[6])
        	problem = ODE(y0, t0, h, f, order)
        else:
        	problem = ODE(y0, t0, h, f)

        calc = Solver(method, steps, problem)
