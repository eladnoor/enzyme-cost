enzyme-cost
===========

Kinetic modeling of pathway fluxes for estimating the specific cost in enzymes
for sustaining a given flux.

Based on the idea of taking a kinetic model, and assuming fluxes are fixed but
the enzyme and metabolite levels are variables. This way, we minimize the total
enzyme level to get the minimal enzyme cost needed to achieve a certain 
pathway flux.

Dependencies:
- Ubuntu packages: (sudo apt-get install X)
    - liblapack-dev (3.6.0)
- PyPI packages:   (sudo pip install X)
    - cvxopt (1.1.9)
    - cvxpy  (0.4.8)
    - ecos   (2.0.4)
    
- SBtab: (git clone https://github.com/derHahn/SBtab.git)
