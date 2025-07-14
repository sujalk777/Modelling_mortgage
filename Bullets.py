import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def Bullet(rate, notonal, periods, CPR):
    # it returns a matrix M such that
    # M=  [0,       1,            2,           3,                4              5]
    # M = [t  notional(t)  prepayment(t)  notional_quote(t)  interest_(t)  installment(t)]
    # WARNING! here "rate" and "periods" are quite general, the choice of getting year/month/day.. steps, depends on the rate
    # that the function receives. So, it is necessary to pass the correct rate to the function
    M = np.zeros((periods + 1, 6))
    M[:, 0] = np.arange(periods + 1)  # time
    M[0, 1] = notonal
    for t in range(1, periods):
        M[t, 4] = rate * M[t - 1, 1]
        M[t, 3] = 0
        Scheduled_outstanding = M[t - 1, 1] - M[t, 3]
        M[t, 2] = Scheduled_outstanding * CPR
        M[t, 1] = Scheduled_outstanding - M[t, 2]  # notional(t) = notional(t-1) - (repayment + prepayment)
        M[t, 5] = M[t, 4] + M[t, 2] + M[t, 3]
    M[periods, 4] = rate * M[periods - 1, 1]
    M[periods, 3] = M[periods - 1, 1]
    M[periods, 5] = M[t, 4] + M[t, 2] + M[t, 3]

    return M


def main():
    notonal = 1000000
    r = 0.05
    T_end = 30
    Lambda = 0.01
    M = Bullet(r, notonal, T_end, Lambda)
    for i in range(0, T_end + 1):
        print(
            "Ti={0}, Notional={1:.0f}, Prepayment={2:.0f}, Notional Repayment={3:.0f}, Interest Rate={4:.0f}, Installment={5:.0f} ".format(
                M[i, 0], M[i, 1], M[i, 2], M[i, 3], M[i, 4], M[i, 5]))
    plt.figure(1)
    plt.title("notonal of Bullets")
    plt.plot(M[:, 0], M[:, 1], '-r')
    plt.grid()
    plt.xlabel('time')
    plt.ylabel('notional')


main()



