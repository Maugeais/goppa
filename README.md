<a name="readme-top"></a>

## About The Project

Implementation in Python of a complete liabrary handling goppa codes, from creation of the code to decoding. It is developped for pedagogical purposes, emphasis is therefore made on readebility instead of speed.

The main library

```
    from goppa import goppa
    import numpy as np
    import time

    print("Génération du code de Goppa...")

    G = goppa.goppa(6, 5, 2**6)
```
