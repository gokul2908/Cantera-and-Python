'''
Adiabatic flame temperature for methane
CH4 + 2(O2 + 3.76N2) = 0.8CO2 + 1.6H20 + 7.52N2 + O2
for lean mixture
CH4 + 2/ER(O2+3.76N2) -> CO2 + 2H2O + (7.52/ER)N2 + (2/Φ- 2)O2
for rich mixture
CH4 + 2/ER(O2+3.76N2) -> ( 4/ER–3)CO2 + 2H2O + (7.52/ER)N2 + (4- 4/ER)CO
'''
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


def h(T, coeff): # fuction to determine the enthalpy
	R = 8.314 # jol/mol-k
	a1 = coeff[0]
	a2 = coeff[1]
	a3 = coeff[2]
	a4 = coeff[3]
	a5 = coeff[4]
	a6 = coeff[5]
	a7 = coeff[6]

	return (a1 + a2*T/2 + a3*pow(T,2)/3 + a4*pow(T,3)/4 + a5*pow(T,4)/5 + a6/T)*R*T

def f(T,ER): # fuction to determine total enthalpy of product & reactant

	#enthalpy of reactants
	t_sat   = 298.15
	h_CH4_r = h(t_sat,CH4_coefficient_l)
	h_O2_r  = h(t_sat,O2_coefficient_l)
	h_N2_r  = h(t_sat,N2_coefficient_l)

	H_reactants = h_CH4_r + (2/ER)*(h_O2_r + 3.76*h_N2_r)	

	#enthalpy of produce
	h_CO2_p = h(T,CO2_coefficient_h)
	h_N2_p  = h(T,N2_coefficient_h)
	h_H2O_p = h(T,H2O_coefficient_h)
	#h_CH4_p = h(T,CH4_coefficient_h)
	h_O2_p  = h(T,O2_coefficient_h)
	h_CO_p  = h(T,CO_coefficient_h)

	if ER==1: #Stoichiometric 
		H_product   = h_CO2_p + 7.52*h_N2_p + 2*h_H2O_p + 0.4*h_O2_p
	else:
	 	if ER>1: #rich mixture
	 		H_product   = ((4/ER)-3)*h_CO2_p + (7.52/ER)*h_N2_p + 2*h_H2O_p + (4-(4/ER))*h_CO_p
	 	else: 	 #lean mixture
	 		H_product   = h_CO2_p + (7.52/ER)*h_N2_p + 2*h_H2O_p + ((2/ER)-2)*h_O2_p

	return H_product-H_reactants

def fprime(T,ER): #Numerical differntiation
	return (f(T+1e-6,ER)-f(T,ER))/1e-6

#Methane Coefficient from Nasa polynomials data
CH4_coefficient_l = [ 5.14987613E+00, -1.36709788E-02,  4.91800599E-05, -4.84743026E-08,  1.66693956E-11, -1.02466476E+04, -4.64130376E+00]
#Oxygen 
O2_coefficient_l  = [ 3.78245636E+00, -2.99673416E-03,  9.84730201E-06, -9.68129509E-09,  3.24372837E-12, -1.06394356E+03,  3.65767573E+00]
O2_coefficient_h  = [ 3.28253784E+00,  1.48308754E-03, -7.57966669E-07,  2.09470555E-10, -2.16717794E-14, -1.08845772E+03,  5.45323129E+00]
#Nitrogen
N2_coefficient_h  = [ 0.02926640E+02,  0.14879768E-02, -0.05684760E-05,  0.10097038E-09, -0.06753351E-13, -0.09227977E+04,  0.05980528E+02]
N2_coefficient_l  = [ 0.03298677E+02,  0.14082404E-02, -0.03963222E-04,  0.05641515E-07, -0.02444854E-10, -0.10208999E+04,  0.03950372E+02]
#Carbon-di-oxide
CO2_coefficient_h = [ 3.85746029E+00,  4.41437026E-03, -2.21481404E-06,  5.23490188E-10, -4.72084164E-14, -4.87591660E+04,  2.27163806E+00]
#Carbon-mono-oxide
CO_coefficient_h  = [ 2.71518561E+00,  2.06252743E-03, -9.98825771E-07,  2.30053008E-10, -2.03647716E-14, -1.41518724E+04,  7.81868772E+00]
#Water vapor
H2O_coefficient_h = [ 3.03399249E+00,  2.17691804E-03, -1.64072518E-07, -9.70419870E-11,  1.68200992E-14, -3.00042971E+04,  4.96677010E+00]

ER_x      = np.linspace(0.1,2,100) #Equivalence ratio     
T_guess = 2000
tol     = 1e-5
alpha   = 1
temp = []  # empty array for newton rhapson technique
temp_cantera = [] # empty array for cantera
gas = ct.Solution('gri30.xml') # centera

for ER in ER_x:
	while (abs(f(T_guess,ER))>tol ):
		T_guess = T_guess - alpha*(f(T_guess,ER)/fprime(T_guess,ER)) # newton rhapson technique
	
	gas.TPX = 298.15,101325,{"CH4":1,"O2":2/ER,"N2":(2*3.76/ER)}
	gas.equilibrate("HP","auto")	
	print("Adiabatic Flame temperature From newton rhapson technique is ",+ T_guess)
	print("Adiabatic Flame temperature From Cantera is ",+ gas.T)
	temp.append(T_guess)
	temp_cantera.append(gas.T)
		

plt.plot(ER_x,temp,color="red")
plt.plot(ER_x,temp_cantera,color="blue")
plt.xlabel("Equivalence ratio")
plt.ylabel("Adiabatic flame temperature")
plt.grid("on")
plt.title("Adiabatic flame temperature")
plt.legend(["Newton rhapson technique","cantera"])
plt.show()

print(temp_cantera.index(max(temp_cantera)))
aa = temp_cantera.index(max(temp_cantera))
ER = ER_x[aa]
print(ER)
