# Algoritmes per al càlcul numèric

### 1. Resolució de sistemes triangulars

1. **Substitució enrere (`triU`)** — per a matrius triangulars superiors.
2. **Substitució endavant (`triL`)** — per a matrius triangulars inferiors.
3. **Substitució endavant amb diagonal unitària (`triL_uns_diagonal`)**.

### 2. Eliminació de Gauss i factoritzacions

4. **Eliminació de Gauss sense pivotatge** (`elimGaussresolucio`) — transforma el sistema en triangular i després aplica `triU`.
5. **Factorització LU** (`factLU`) — obté les matrius `L` i `U`.
6. **Factorització de Doolittle** (`doolittle`) — factorització `A = LU`, amb diagonal de `L` igual a 1.
7. **Factorització PA = LU amb pivotatge parcial** (`factPALU`).
8. **Pivotatge parcial esglaonat** (`pivotatge_parcial_esglaonat`) — inclou també el recompte de permutacions per calcular el determinant.
9. **Càlcul de la inversa mitjançant LU** (`inversaLU`) — resolent diversos sistemes triangulars.

### 3. Factoritzacions per a matrius especials

10. **Factorització de Cholesky** (`cholesky`) — per a matrius simètriques definides positives.
11. **Factorització `A = LDLᵀ`**

### 4. Factorització QR

12. **Gram-Schmidt modificat** (`gram_schmidt`) — calcula `A = QR`, amb `Q` ortonormal i `R` triangular superior.
13. **Transformacions de Householder** (`matHouseholder`, `factHouseholder`) — altra manera d'obtenir la factorització QR.
14. **Householder per resoldre sistemes** (`Householder`) — aplica les transformacions a `A` i `b`, obté un sistema triangular i el resol amb `triU`.

### 5. Normes

15. **Norma 1 vectorial** (`norma1`)
16. **Norma 2 vectorial** (`norma2`)
17. **Norma infinit vectorial** (`normainf`)
18. **Norma p** (`normap`)
19. **Norma 1 matricial** (`norma1Matricial`)
20. **Norma infinit matricial** (`normaInfMatricial`)
21. **Norma de Frobenius** (`normaFrobMatricial`)

### 6. Mètodes iteratius per a sistemes lineals

22. **Mètode de Jacobi** (`jacobi`) — mètode iteratiu per resoldre sistemes lineals.
23. **Mètode de Gauss-Seidel** (`gaussseidel`) — utilitza immediatament els valors actualitzats.
24. **Mètode de relaxació / SOR** (`omegarelaxacio`) — introdueix el paràmetre `ω`.
25. **Mètode del gradient** (`Gradient`) — per a matrius simètriques definides positives.
26. **Mètode del gradient conjugat** — per a matrius simètriques definides positives.
