# mrpbounds

Computing sharp bounds on multilevel regression with poststratification (MRP) estimates when the joint distribution of poststratification variables is only partially known. Given marginal or joint probability tables, the package uses linear programming to find the minimum and maximum possible values of MRP quantities, cell probabilities, and mean squared error.

### Install R Package

You can install the latest version from Github using:

```R
remotes::install_github("sumtxt/mrpbounds")
```

You may need to install the `remotes` package first.

### Usage

```R
library(mrpbounds)

# Fit a model on survey data
svy <- data.frame(
  A = rbinom(1000, 1, 0.3),
  B = rbinom(1000, 1, 0.1),
  C = rbinom(1000, 1, 0.5)
)
model <- glm(A ~ B + C, data = svy, family = binomial())

# Define marginal/joint constraints from the target population
joint_ab <- data.frame(A = c(0, 0, 1, 1), B = c(0, 1, 0, 1),
                       prob = c(0.1, 0.2, 0.3, 0.4))
margin_c <- data.frame(C = c(0, 1), prob = c(0.6, 0.4))

# Compute bounds on the MRP estimate
bounds <- mrp_bounds(model, list(joint_ab, margin_c))
```
