class ODE:
    def __init__(self, y, t, h, func, order=1):
        self.y = y
        self.t = t
        self.h = h
        self.steps = steps
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

    # def adam_bashforth(self):
#	return y + h*func

    # def adam_multon(self):
#	return y + h*func


def main():
    # Lendo o arquivo de entrada
    inFile = open("entrada.txt", "r")
    outFile = open("saida.txt", "w")
    lines = read.readlines()
