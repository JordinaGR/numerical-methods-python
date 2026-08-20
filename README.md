# Càlcul Numèric

Aquest repositori conté diferents implementacions d'algoritmes de **càlcul numèric** en Python i un notebook d'**àlgebra lineal numèrica**. Els continguts estan organitzats en els següents blocs:

### 1. Equacions no lineals

1. **Mètode de Newton-Raphson** — resolució numèrica d'equacions no lineals.

   * Fitxer: `newton_raphson_eq_no_lineals.py`

2. **Mètode de Newton per a sistemes d'equacions no lineals** — resolució de sistemes mitjançant el vector de funcions i la seva matriu jacobiana.

   * Fitxer: `metode_newton_eq_no_lineals.py`

3. **Mètode de Newton-Raphson per a sistemes no lineals i Newton-Raphson modificat**.

   * Fitxer: `newton_raphson_eq_no_lineals.py`

### 2. Mínims quadrats i interpolació

4. **Mètode de mínims quadrats** — aproximació d'un conjunt de dades mitjançant polinomis de diferents graus.

5. **Interpolació polinòmica** — interpolació a partir d'un conjunt de punts.

6. **Comparació entre mínims quadrats i interpolació**, incloent-hi el càlcul i estudi dels errors.

   * Fitxer: `minims_quadrats_i_interpolacio.py`

### 3. Splines

7. **Splines cúbics** — construcció de splines cúbics a partir d'un conjunt de punts.

8. **Càlcul de derivades als punts** — aproximació de les derivades mitjançant diferències finites.

9. **Spline cúbic de tipus Hermite (`C1`)** — construcció a partir dels valors de la funció i les seves derivades.

10. **Spline cúbic (`C2`)** — construcció mitjançant `CubicSpline`.

* Fitxer: `splines.py`

### 4. Quadratura numèrica

11. **Quadratura de Gauss** — càlcul dels nodes i pesos de Gauss per a diferents nombres de punts d'integració.

* Fitxer: `quadraturaGauss.py`

12. **Quadratures numèriques** — implementació i estudi de diferents tècniques de quadratura.

* Fitxer: `quadratures.py`

13. **Regla del trapezi composta**.

14. **Regla de Simpson composta**.

15. **Quadratura de Gauss composta**.

* Fitxer: `quadratures_trapeziCompost_SimpsonComp_GaussComp.py`

### 5. Equacions diferencials ordinàries (EDOs)

16. **Mètode d'Euler explícit** — resolució numèrica de problemes de valor inicial.

17. **Mètode d'Euler enrere** — implementació implícita del mètode d'Euler.

18. **Mètode de Runge-Kutta de quart ordre (RK4)**.

* Fitxer: `EDOs_RungeKutta4_Euler.py`

19. **Mètode d'Euler per a sistemes d'EDOs**.

* Fitxer: `EulerEdos.py`

20. **Mètode de Heun** — mètode de segon ordre basat en una predicció i una correcció.

* Fitxer: `Edo_Heune.py`

21. **Mètode de tir (Shooting Method)** — resolució d'un problema de contorn mitjançant la conversió a un problema de valor inicial i l'ajust de les condicions inicials.

* Fitxer: `shootingMethod.py`

### 6. Equacions diferencials parcials (EDPs)

22. **Equació de difusió 1D** — discretització espacial mitjançant diferències finites i evolució temporal.

23. **Mètode d'Euler explícit per a EDPs**.

24. **Condicions de Dirichlet i Neumann** — tractament de diferents condicions de contorn.

25. **Mètode explícit per a condicions de Dirichlet**.

26. **Mètode implícit per a condicions de Dirichlet**.

* Fitxer: `EDPs_Neumann_Dirichlet.py`

27. **Discretització de l'equació de difusió amb condicions de Dirichlet i Neumann**.

28. **Mètode d'Euler cap endavant per a l'evolució temporal**.

* Fitxer: `EDPs_euler.py`

29. **Mètode de diferències finites per a l'equació de difusió 1D**.

* Fitxer: `FD1D_diffusionEq.py`

---

# Àlgebra lineal numèrica

## 1. Resolució de sistemes triangulars

1. **Substitució enrere (`triU`)** — per a matrius triangulars superiors.
2. **Substitució endavant (`triL`)** — per a matrius triangulars inferiors.
3. **Substitució endavant amb diagonal unitària (`triL_uns_diagonal`)**.

### 2. Eliminació de Gauss i factoritzacions

4. **Eliminació de Gauss sense pivotatge (`elimGaussresolucio`)** — transforma el sistema en triangular i després aplica `triU`.
5. **Factorització LU (`factLU`)** — obté les matrius `L` i `U`.
6. **Factorització de Doolittle (`doolittle`)** — factorització `A = LU`, amb diagonal de `L` igual a 1.
7. **Factorització PA = LU amb pivotatge parcial (`factPALU`)**.
8. **Pivotatge parcial esglaonat (`pivotatge_parcial_esglaonat`)** — inclou també el recompte de permutacions per calcular el determinant.
9. **Càlcul de la inversa mitjançant LU (`inversaLU`)** — resolent diversos sistemes triangulars.

### 3. Factoritzacions per a matrius especials

10. **Factorització de Cholesky (`cholesky`)** — per a matrius simètriques definides positives.
11. **Factorització `A = LDLᵀ`**.

### 4. Factorització QR

12. **Gram-Schmidt modificat (`gram_schmidt`)** — calcula `A = QR`, amb `Q` ortonormal i `R` triangular superior.
13. **Transformacions de Householder (`matHouseholder`, `factHouseholder`)** — altra manera d'obtenir la factorització QR.
14. **Householder per resoldre sistemes (`Householder`)** — aplica les transformacions a `A` i `b`, obté un sistema triangular i el resol amb `triU`.

### 5. Normes

15. **Norma 1 vectorial (`norma1`)**.
16. **Norma 2 vectorial (`norma2`)**.
17. **Norma infinit vectorial (`normainf`)**.
18. **Norma p (`normap`)**.
19. **Norma 1 matricial (`norma1Matricial`)**.
20. **Norma infinit matricial (`normaInfMatricial`)**.
21. **Norma de Frobenius (`normaFrobMatricial`)**.

### 6. Mètodes iteratius per a sistemes lineals

22. **Mètode de Jacobi (`jacobi`)** — mètode iteratiu per resoldre sistemes lineals.
23. **Mètode de Gauss-Seidel (`gaussseidel`)** — utilitza immediatament els valors actualitzats.
24. **Mètode de relaxació / SOR (`omegarelaxacio`)** — introdueix el paràmetre `ω`.
25. **Mètode del gradient (`Gradient`)** — per a matrius simètriques definides positives.
26. **Mètode del gradient conjugat** — per a matrius simètriques definides positives.
