import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# shooting method d'una edo amb valors inicials [x, x, - - ] i finals [x, x, -, -]

def odef(t, y):
    return np.array((y[1], y[2], y[3], np.exp(t)*y[1] - y[3]))

def RK4(odef, x, y, h):
    k1 = odef(x, y)
    k2 = odef(x + h/2, y + (h/2)*k1)
    k3 = odef(x + h/2, y + (h/2)*k2)
    k4 = odef(x + h, y + h*k3)
    return y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

def RK4_solve(odef, initialTime, finalTime, y0, nOfSteps):
    h = (finalTime - initialTime) / nOfSteps
    t = np.linspace(initialTime, finalTime, nOfSteps+1)
    
    y = np.zeros((nOfSteps+1, len(y0)))
    y[0,:] = y0

    for i in range(nOfSteps):
        y[i+1,:] = RK4(odef, t[i], y[i,:], h)

    return t, y

initialTime = 0.
finalTime = 1.
nOfSteps = 100 

# Funció que volem fer zero. 
# 'guess' conté les incògnites inicials: [y''(0), y'''(0)]
def shooting_objective(guess):
    # Construïm el vector inicial complet:
    # Sabem: y(0)=0, y'(0)=1
    # No sabem: y''(0)=guess[0], y'''(0)=guess[1]
    y0_attempt = np.array([0., 1., guess[0], guess[1]])
    _, y_sol = RK4_solve(odef, initialTime, finalTime, y0_attempt, nOfSteps)
    
    # Agafem els valors al final de l'interval (última fila)
    y_final_calc = y_sol[-1, 0]  # y(1) calculat
    yp_final_calc = y_sol[-1, 1] # y'(1) calculat
    
    # Condicions de contorn desitjades a x=1: y(1)=0 i y'(1)=-1
    # Retornem la diferència (Residu)
    return [y_final_calc - 0., yp_final_calc - (-1.)]

# A) Trobar les condicions inicials correctes
# Fem una suposició inicial per [y''(0), y'''(0)]
initial_guess = [0., 0.] 
roots = fsolve(shooting_objective, initial_guess) # funció a minimitzar, ini vals

print(f"Valors optimitzats trobats: y''(0)={roots[0]}, y'''(0)={roots[1]}")

# B) Executar la simulació final amb els valors correctes
y0_opt = np.array([0., 1., roots[0], roots[1]])
t, y = RK4_solve(odef, initialTime, finalTime, y0_opt, nOfSteps)

print(y[-1, :])

# 4. RESULTATS I GRÀFICA
# Trobar y(0.5). Com que nOfSteps=100, l'índex 50 correspon exactament a t=0.5
idx_05 = int(nOfSteps / 2) 
val_05 = y[idx_05, 0]

print(f"La solució y(0.5) amb 3 xifres significatives és: {val_05:.3g} (Valor complet: {val_05:.5f})")

plt.figure(figsize=(8, 5))
plt.title(f"Solució del problema de contorn (Shooting Method)\ny(0.5) = {val_05:.4f}")
plt.plot(t, y[:, 0], label="y(x)", linewidth=2)
plt.plot(t, y[:, 1], label="y'(x)", linestyle="--")
plt.axvline(0.5, color='gray', linestyle=':', label='x=0.5')
plt.scatter([0.5], [val_05], color='red', zorder=5) # Punt clau
plt.xlabel("x")
plt.ylabel("Solució")
plt.legend()
plt.grid(True)
plt.show()
