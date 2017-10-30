enzyme-cost
===========

Kinetic modeling of pathway fluxes for estimating the specific cost in enzymes
for sustaining a given flux.

Based on the idea of taking a kinetic model, and assuming fluxes are fixed but
the enzyme and metabolite levels are variables. This way, we minimize the total
enzyme level to get the minimal enzyme cost needed to achieve a certain 
pathway flux.

Dependencies:
- Ubuntu packages: `sudo apt install X`
    - glpk-utils (4.63-1)

- PyPI packages:   `sudo pip install X`
    - numpy (1.13.3)
    - scipy (1.0.0)
    - matplotlib (2.1.0)
    - pulp (1.6.8)
    - pandas (0.21.0)
    - tablib (0.12.1)
    
- Component Contribution (optional):
```
git clone https://github.com/eladnoor/component-contribution.git
```


Example
-------
Try running the example script:
```
python -m example.test
```
