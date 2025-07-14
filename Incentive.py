import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def Annuity(rate, notonal, periods, CPR):
    # it returns a matrix M such that
    # M=  [0,       1,            2,           3,                4              5]
    # M = [t  notional(t)  prepayment(t)  notional_quote(t)  interest_(t)  installment(t)]
    # WARNING! here "rate" and "periods" are quite general, the choice of getting year/month/day.. steps, depends on the rate
    # that the function receives. So, it is necessary to pass the correct rate to the function
    M = np.zeros((periods + 1, 6))
    M[:, 0] = np.arange(periods + 1)  # time
    M[0, 1] = notonal
    for t in range(1, periods):
        remaing_period = periods - (t - 1)
        # installmet C(t_i)
        M[:, 5] = rate * M[t - 1, 1] / (1 - 1 / (1 + rate) ** remaing_period)
        # intrest payemnt I(t_i)=K*N(t_{i})
        M[:, 4] = rate * M[t - 1, 1]
        # notoanl payment Q(t_i)=C(t_i)-I(t-i)
        M[t, 3] = M[t, 5] - M[t, 4]

        # prepament= P(t_i)=lambda*(N(t_i)-Q_(ti))
        M[t, 2] = CPR * (M[t - 1, 1] - M[t, 3])
        # notional, N(t_{i+1}) = N(t_{i}) - lambda * (Q(t_{i} + P(t_i)))
        M[t, 1] = M[t - 1, 1] - M[t, 3] - M[t, 2]

    return M


def main():
    IncentiveFunction = lambda x: 0.04 + 0.1 / (1 + np.exp(115 * (0.02 - x)))

    oldRate = 0.05
    newRate = np.linspace(-0.05, 0.15, 25)

    epsilon = oldRate - newRate
    incentive = IncentiveFunction(epsilon)

    plt.figure(1)
    plt.plot(newRate, incentive)
    plt.xlabel('S(t)')
    plt.ylabel('Incentive')
    plt.grid()

    plt.figure(2)
    plt.plot(epsilon, incentive)
    plt.xlabel('epsilon= K - S(t)')
    plt.ylabel('Incentive')
    plt.grid()


main()



