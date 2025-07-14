
# Modeling Mortgage



https://github.com/user-attachments/assets/37375059-65b5-45ca-a6ef-f0c79a052db5

This repository contains Python code for modeling mortgage products using stochastic interest rate models and various mortgage payment schemes. In this project, we simulate interest rate paths using the Hull–White model and use these paths to compute swap rates and mortgage profiles under different assumptions. Two primary mortgage models are implemented: the Bullet Mortgage and the Annuity Mortgage, each incorporating a stochastic prepayment mechanism.

---

## Table of Contents

- [Introduction](#introduction)
- [Overview](#overview)
- [Mathematical Models and Formulas](#mathematical-models-and-formulas)
  - [Hull–White Interest Rate Model](#hull–white-interest-rate-model)
  - [Zero-Coupon Bond Pricing](#zero-coupon-bond-pricing)
  - [Swap Rate Calculation](#swap-rate-calculation)
  - [Mortgage Models](#mortgage-models)
    - [Bullet Mortgage Model](#bullet-mortgage-model)
    - [Annuity Mortgage Model](#annuity-mortgage-model)
  - [Stochastic Prepayment Modeling](#stochastic-prepayment-modeling)
- [Implementation Details](#implementation-details)
- [Installation and Requirements](#installation-and-requirements)
- [Usage](#usage)
- [Results](#results)
- [License](#license)

---

## Introduction

The objective of this project is to model mortgage products by integrating stochastic interest rate dynamics with mortgage payment structures. The interest rates are simulated using the Hull–White model, and these rates are then used to derive swap rates. Mortgage profiles (i.e., outstanding notional balances) are computed for two types of mortgage products:
- **Bullet Mortgage:** where only interest and stochastic prepayments are made until a lump-sum principal repayment at maturity.
- **Annuity Mortgage:** where fixed installment payments (comprising both interest and scheduled principal repayment) are made, and additional prepayments may occur.

A stochastic prepayment rate is determined by an incentive function, which captures the borrower's incentive to prepay based on market conditions.

---

## Overview

The repository is structured into several key components:

- **Interest Rate Simulation:**  
  The `GeneratePathsHWEuler` function simulates multiple interest rate paths using the Euler discretization of the Hull–White model.

- **Hull–White Model Functions:**  
  Functions such as `HW_theta`, `HW_A`, `HW_B`, and `HW_ZCB` compute various components of the Hull–White model, including the short rate, zero-coupon bond prices, and adjustment factors.

- **Swap Rate Calculation:**  
  The `SwapRateHW` function calculates the swap rate based on simulated interest rate paths.

- **Mortgage Modeling:**  
  Two functions, `Bullet` and `Annuity`, simulate the evolution of the outstanding mortgage balance over time under different payment schemes, incorporating scheduled payments and stochastic prepayments.

- **Prepayment Incentive:**  
  A prepayment incentive function is defined to model the borrower's behavior, influencing the prepayment rate.

---

## Mathematical Models and Formulas

### Hull–White Interest Rate Model

The Hull–White model is used to simulate the evolution of the short-term interest rate, $r_t$, according to the stochastic differential equation:
$$
dr_t = \lambda \left( \theta(t) - r_t \right) dt + \eta \, dW_t,
$$
where:
- $\lambda$ is the mean-reversion speed,
- $\theta(t)$ is the time-dependent long-term mean level,
- $\eta$ is the volatility,
- $W_t$ is a standard Brownian motion.

The function $\theta(t)$ is defined as:
$$
\theta(t) = \frac{1}{\lambda} \frac{d f(0,t)}{dt} + f(0,t) + \frac{\eta^2}{2\lambda^2} \left( 1 - e^{-2\lambda t} \right),
$$
with $f(0,t)$ representing the instantaneous forward rate.

### Zero-Coupon Bond Pricing

Within the Hull–White framework, the price of a zero-coupon bond at time $t$ maturing at time $T$ is given by:
$$
P(t,T) = \exp\left( A(t,T) - B(t,T) r_t \right),
$$
where the functions $A(t,T)$ and $B(t,T)$ are calculated as:
- $$
  B(t,T) = \frac{1}{\lambda} \left( e^{-\lambda (T-t)} - 1 \right)
  $$
- $$
  A(t,T) = \text{(a function of } \theta(t), \eta, \lambda \text{, and the time difference } T-t \text{)}
  $$
The exact expression for $A(t,T)$ is computed numerically in the code using integration.

### Swap Rate Calculation

The swap rate, $S(t)$, is calculated from the prices of zero-coupon bonds as follows:
$$
S(t) = \frac{P(t, T_i) - P(t, T_m)}{\sum_{j=1}^{n} \tau \, P(t, T_j)},
$$
where:
- $P(t, T)$ is the zero-coupon bond price,
- $T_i$ and $T_m$ are the start and end dates of the swap,
- $n$ is the number of payment dates,
- $\tau$ is the time interval between payments.

### Mortgage Models

#### Bullet Mortgage Model

In the Bullet Mortgage, the borrower makes only interest payments along with any prepayments until the maturity date when the entire remaining principal is due. The key formulas are:

- **Interest Payment at time $t$:**
$$
I_t = r \cdot N_{t-1},
$$
where $r$ is the interest rate and $N_{t-1}$ is the outstanding notional at the previous time step.

- **Prepayment at time $t$:**
$$
P_t = \text{CPR}_t \times N_{t-1},
$$
where $\text{CPR}_t$ is the prepayment rate at time $t$.

- **Outstanding Notional Update:**
$$
N_t = N_{t-1} - P_t.
$$

At maturity, the full remaining principal is repaid.

#### Annuity Mortgage Model

The Annuity Mortgage involves fixed installment payments that cover both interest and principal. The formulas used are:

- **Fixed Installment Payment:**
$$
C_t = \frac{r \cdot N_{t-1}}{1 - (1 + r)^{-m}},
$$
where $m$ is the number of remaining periods.

- **Interest Component:**
$$
I_t = r \cdot N_{t-1}.
$$

- **Scheduled Principal Repayment:**
$$
Q_t = C_t - I_t.
$$

- **Prepayment:**
$$
P_t = \text{CPR}_t \times \left( N_{t-1} - Q_t \right).
$$

- **Notional Update:**
$$
N_t = N_{t-1} - Q_t - P_t.
$$

### Stochastic Prepayment Modeling

Prepayments are modeled stochastically using an incentive function that depends on the difference between a target rate $K$ and the current swap rate $S(t)$. The incentive function is defined as:
$$
\text{Incentive}(x) = 0.04 + \frac{0.1}{1 + e^{-200x}},
$$
where:
$$
x = K - S(t).
$$
This function is designed to capture the idea that a lower current rate relative to $K$ provides a higher incentive to prepay.

---

## Implementation Details

- **Interest Rate Simulation:**  
  The function `GeneratePathsHWEuler` simulates multiple interest rate paths using the Euler method applied to the Hull–White SDE. It computes the forward rate, the mean reversion level $\theta(t)$, and updates the short rate path over discrete time steps.

- **Hull–White Model Components:**  
  Functions like `HW_theta`, `HW_A`, `HW_B`, and `HW_ZCB` are used to calculate the parameters and bond prices according to the Hull–White framework.

- **Swap Rate and Mortgage Modeling:**  
  The `SwapRateHW` function computes the swap rate from the simulated interest rate paths, which in turn feeds into the mortgage models (`Bullet` and `Annuity`). These mortgage functions compute the cash flows and outstanding balances over time while incorporating both scheduled payments and stochastic prepayments.

- **Prepayment Incentive:**  
  The prepayment rate, often denoted as $\text{CPR}_t$, is determined by the incentive function. This rate adjusts the mortgage balance according to the borrower’s prepayment behavior under changing market conditions.

---

## Installation and Requirements

### Requirements

- Python 3.x
- [NumPy](https://numpy.org/)
- [SciPy](https://www.scipy.org/)
- [Matplotlib](https://matplotlib.org/)
- (Optional) [Streamlit](https://streamlit.io/) for interactive visualizations

### Installation

Clone the repository and install the required packages:
```bash
git clone https://github.com/shubh123a3/Modeling-Mortgage.git
cd Modeling-Mortgage
pip install -r requirements.txt
```

---

## Usage

To run the simulations and generate the plots, simply execute the main application:
```bash
python app.py
```
If using Streamlit for interactive exploration, run:
```bash
streamlit run app.py
```

The application will:
- Simulate interest rate paths using the Hull–White model.
- Calculate swap rates from the simulated paths.
- Generate and plot the mortgage profiles under both the Bullet and Annuity models.
- Display the impact of stochastic prepayments driven by the incentive function.

---

## Results

The outputs of the simulation include:
- **Incentive Function Plots:** Illustrating how the prepayment incentive changes with the difference between the target rate $K$ and the current swap rate $S(t)$.
- **Swap Rate Distribution:** Histograms showing the distribution of swap rates at the end of the simulation horizon.
- **Mortgage Notional Profiles:** Time series plots comparing the outstanding notional under the Bullet and Annuity mortgage models, both with and without prepayments.

---

## License

This project is licensed under the [MIT License](LICENSE).

Feel free to contribute, open issues, or suggest improvements!

----
