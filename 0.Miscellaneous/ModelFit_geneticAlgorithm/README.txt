DiscreteRelaxationSpectra.exe fits the eta_i and lambda_i parameters to
experimental loss modulus data. The following is the implemented model:

        n
      =====
      \         eta  . omega
       \           i
G'' =   \   --------------------
        /            2         2
       /    1 + omega  . lambda
      /                        i
      =====
      i = 1

Work instructions:
0. Create (and save) a CSV file with two colums.
   lossModulus_masterCurve.csv is provided as an example.
   1st column shall contain frequency (omega) data.
   2nd column shall contain loss modulus (G'') data.
   Columns can or cannot contain headers/column names.
   More columns can be added, but they will be ignored.
1. Run (double click) DiscreteRelaxationSpectra.exe.
   (it takes some seconds to load)
2. Select the newly created CSV.
3. Indicate the number of Maxwell elements to use in the model.
4. Repeat step 3 until satisfied.
5. To close the program, input a negative value or cancel the process.

Outputs:
* For each indicated number of maxwell elements, a plot will be created in the
  same directory as DiscreteRelaxationSpectra.exe.
* The output has the following nomenclature:
  >>> <number of Maxwell elements n> Maxwell Elements
  eta_i:
   [ eta(i) eta(i+1) eta(i+2) ... eta(n-1) eta(n) ]
  lambda_i:
   [ lambda(i) lambda(i+1) lambda(i+2) ... lambda(n-1) lambda(n) ]
  Root Mean Squared Error: <RMSE>
  R-squared:               <R**2>