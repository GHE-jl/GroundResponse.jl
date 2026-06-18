# Exemple of use of the gST_ANN function
# TODO : - Gérer le plot quand g_EWT et t_EWT sont des [NaN] (retourne une erreur)

# 1 - Input parameters
ks=2.13
Cs=2.0e6
kg=1.65
Cg=2.0e6
kp=0.4
Cp=1.9e6
Cf=4.2e6
ri=0.017
ro=0.022
rb=0.08
H=150
Vdot=23.7/1000/60
D=0.029
dt=60 
tf = 1000 

# 2 - Transfer function
# Creating the Entering Water Temperature's timestep and transfer function
t_EWT, g_EWT = gST_ANN(ks, Cs, kg, Cg, kp, Cp, Cf, ri, ro, rb, H, Vdot, D, dt, tf)

# 3 - Creation of the figure
plot(t_EWT/(H^2/(9*ks/Cs)), g_EWT, legend=false)

# quand les paramètres sont out of range, ça retourne une erreur
plot!(xscale=:log10, minorgrid=true) # x axis in log
xlabel!("t/t_s (-)")
ylabel!("g (-)")