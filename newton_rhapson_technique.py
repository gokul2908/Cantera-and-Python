# newton rhapson technique


def f(x):
	return pow(x,3) + 4*pow(x,2) + 2*x -5

def fprime(x):
	return 3*pow(x,2) + 8*x + 2


x_guess = 10
tolerance = 1e-22
alpha = 1.29
iterations = 0

while (f(x_guess)>tolerance):
	x_guess = x_guess - alpha*(f(x_guess)/fprime(x_guess))
	iterations = iterations + 1

print(iterations)
print(x_guess)

