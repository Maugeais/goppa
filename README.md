## About The Project

Implementation in Python of a complete liabrary handling binary goppa codes, from creation of the code to decoding. 

This project is developped for pedagogical purposes, emphasis is therefore made on readebility instead of speed.

It is thought to be as self contained as possible, and contains libraries that handle finite fields, and polynomials over finite fields.

### Built With

* Python
* numpy

### Dependencies therefore

The library goppa needs is built on modules detailed in the following graph. Each one can be used on its own, and usage examples are given in the folder examples.

                                         goppa.py
                                          │   │
                                          │   │
                          ┌───────────────┘   └────────────────┐
                          │                                    │
                          │                                    │
                          │                                    │
              polynomials_finite_field.py                    error_control
                          │
                          │
                          │
                     finite_field.py
                          │
                          │
                          │
                polynomials_prime_field.py


## Getting Started

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/Maugeais/goppa.git
   ```
   or download the zip archive from https://github.com/Maugeais/goppa
   
2. Copy the subfolder goppa into your working directory

<!-- USAGE EXAMPLES -->
## Usage

```
    from goppa import goppa
    
    # Generation of the gopa code over the field GF(2**7), with L containing 2**6 elemnts 
    # chosen randomly in GF(2**7) and an irreducible polynomial of degree 8

    G = goppa.goppa(7, 8, 2**6)
    
    # Creating a code and an error
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
    

    error = np.zeros(G.H.shape[1], dtype = int)
    I = np.random.choice(range(G.H.shape[1]), G.g.deg, replace = False)
    error[I] = 1
        
    # Correcting the code with Patterson algorithm and decoding it      
    mp = G.decode((c+error))
```


<!-- ROADMAP -->
## Roadmap

- [ ] Adding full description of the modules
- [ ] Optimising the code
- [ ] Writting the documentation explainning teh theory

See the [open issues](https://github.com/Maugeais/goppa/issues) for a full list of proposed features (and known issues).


<!-- CONTACT -->
## Contact

Sylvain Maugeais  - sylvain.maugeais@univ-lemans.fr

Project Link: [https://github.com/Maugeais/goppa](https://github.com/Maugeais/goppa)

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* Most of the error correction code theory comes from 
  
  Van Lint. 1998. Introduction to Coding Theory (3rd. ed.). Springer-Verlag, Berlin, Heidelberg.
  
* The most complete version of Patterson's algorithm was found in 
  
  Bernstein, D.J. (2011). List Decoding for Binary Goppa Codes. In: Chee, Y.M., et al. Coding and Cryptology. IWCC 2011. Lecture Notes in Computer Science, vol 6639. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-642-20901-7_4


