"""
    G
    
Gravitational constant in the code units:
- Mass in solar masses (MSun),
- Length in milliparsecs (mpc),
- Time in millions of years (Myr).

Set as `4498502.15157528`.
This value comes from `astropy`:
```python-repl
>>> astropy.__version__
'4.0.1.post1'
```
 
```python-repl
>>> from astropy import units as u
>>> from astropy.constants import G
>>> uG = (u.mpc ** 3)/(u.M_sun * u.Myr ** 2)
>>> G.to(uG).value
4498502.151575286
```
"""
const G = 4498502.151575286

"""
    cvel
    
Speed of light in the code units:
- Mass in solar masses (MSun),
- Length in milliparsecs (mpc),
- Time in millions of years (Myr).

Set as `306601393.7879527`.
This value comes from `astropy`:
```python-repl
>>> astropy.__version__
'4.0.1.post1'
```
 
```python-repl
>>> from astropy import units as u
>>> from astropy.constants import c
>>> uc = (u.mpc/u.Myr)
>>> c.to(uc).value
306601393.7879527
```
"""
const cvel = 306601393.7879527 
