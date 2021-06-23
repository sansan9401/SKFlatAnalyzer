First version of electron corrections. 

```cpp
Aepcor ec; // e.g. as a class member
ec.init("path/to/e_18UL.txt"); //e.g. in a class constructor

//for each electron, k-factor for momentum
double kData = ec.kScaleDT(pt, eta, phi, r9, run); // for data electron

double urnd = gRandom->Rndm(); // uniform between 0 and 1
double kMC0 = ec.kSpreadMC(pt, eta, phi, r9, urnd, genPt); // for MC electron with matched dressed electron (recommended)
double kMC1 = ec.kSmearMC(pt, eta, phi, r9, urnd); // for MC without matched dressed electron
```
