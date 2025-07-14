import streamlit as st
import matplotlib.pyplot as plt

import scipy.integrate as integrate
import numpy as np


def GeneratePathsHWEuler(NoOfPaths, NoOfSteps, T, P0T, lambd, eta):
    # time-step needed for differentiation
    dt = 0.0001
    f0T = lambda t: - (np.log(P0T(t + dt)) - np.log(P0T(t - dt))) / (2 * dt)

    # Initial interest rate is a forward rate at time t->0
    r0 = f0T(0.00001)
    theta = lambda t: 1.0 / lambd * (f0T(t + dt) - f0T(t - dt)) / (2.0 * dt) + f0T(t) + eta * eta / (
                2.0 * lambd * lambd) * (1.0 - np.exp(-2.0 * lambd * t))

    # theta = lambda t: 0.1 +t -t
    # print("changed theta")

    Z = np.random.normal(0.0, 1.0, [NoOfPaths, NoOfSteps])
    W = np.zeros([NoOfPaths, NoOfSteps + 1])
    R = np.zeros([NoOfPaths, NoOfSteps + 1])
    R[:, 0] = r0
    time = np.zeros([NoOfSteps + 1])

    dt = T / float(NoOfSteps)
    for i in range(0, NoOfSteps):
        # making sure that samples from normal have mean 0 and variance 1
        if NoOfPaths > 1:
            Z[:, i] = (Z[:, i] - np.mean(Z[:, i])) / np.std(Z[:, i])
        W[:, i + 1] = W[:, i] + np.power(dt, 0.5) * Z[:, i]
        R[:, i + 1] = R[:, i] + lambd * (theta(time[i]) - R[:, i]) * dt + eta * (W[:, i + 1] - W[:, i])
        time[i + 1] = time[i] + dt

    # Outputs
    paths = {"time": time, "R": R}
    return paths


def HW_theta(lambd, eta, P0T):
    dt = 0.0001
    f0T = lambda t: - (np.log(P0T(t + dt)) - np.log(P0T(t - dt))) / (2 * dt)
    theta = lambda t: 1.0 / lambd * (f0T(t + dt) - f0T(t - dt)) / (2.0 * dt) + f0T(t) + eta * eta / (
                2.0 * lambd * lambd) * (1.0 - np.exp(-2.0 * lambd * t))
    return theta


def HW_A(lambd, eta, P0T, T1, T2):
    tau = T2 - T1
    zGrid = np.linspace(0.0, tau, 250)
    B_r = lambda tau: 1.0 / lambd * (np.exp(-lambd * tau) - 1.0)
    theta = HW_theta(lambd, eta, P0T)
    temp1 = lambd * np.trapz(theta(T2 - zGrid) * B_r(zGrid), zGrid)

    temp2 = eta * eta / (4.0 * np.power(lambd, 3.0)) * (
                np.exp(-2.0 * lambd * tau) * (4 * np.exp(lambd * tau) - 1.0) - 3.0) + eta * eta * tau / (
                        2.0 * lambd * lambd)

    return temp1 + temp2


def HW_B(lambd, eta, T1, T2):
    return 1.0 / lambd * (np.exp(-lambd * (T2 - T1)) - 1.0)


def HW_ZCB(lambd, eta, P0T, T1, T2, rT1):
    n = np.size(rT1)

    if T1 < T2:
        B_r = HW_B(lambd, eta, T1, T2)
        A_r = HW_A(lambd, eta, P0T, T1, T2)
        return np.exp(A_r + B_r * rT1)
    else:
        return np.ones([n])


def SwapRateHW(t, Ti, Tm, n, r_t, P0T, lambd, eta):
    # CP- payer of receiver
    # n- notional
    # K- strike
    # t- today's date
    # Ti- beginning of the swap
    # Tm- end of Swap
    # n- number of dates payments between Ti and Tm
    # r_t -interest rate at time t

    if n == 1:
        ti_grid = np.array([Ti, Tm])
    else:
        ti_grid = np.linspace(Ti, Tm, n)
    tau = ti_grid[1] - ti_grid[0]

    # overwrite Ti if t>Ti
    prevTi = ti_grid[np.where(ti_grid < t)]
    if np.size(prevTi) > 0:  # prevTi != []:
        Ti = prevTi[-1]

    # Now we need to handle the case when some payments are already done
    ti_grid = ti_grid[np.where(ti_grid > t)]

    temp = np.zeros(np.size(r_t));

    P_t_TiLambda = lambda Ti: HW_ZCB(lambd, eta, P0T, t, Ti, r_t)

    for (idx, ti) in enumerate(ti_grid):
        if ti > Ti:
            temp = temp + tau * P_t_TiLambda(ti)

    P_t_Ti = P_t_TiLambda(Ti)
    P_t_Tm = P_t_TiLambda(Tm)

    swapRate = (P_t_Ti - P_t_Tm) / temp

    return swapRate


def Bullet(rate, notional, periods, CPR):
    # it returns a matrix M such that
    # M = [t  notional(t)  prepayment(t)  notional_quote(t)  interest_(t)  installment(t)]
    # WARNING! here "rate" and "periods" are quite general, the choice of getting year/month/day.. steps, depends on the rate
    # that the function receives. So, it is necessary to pass the correct rate to the function
    M = np.zeros((periods + 1, 6))
    M[:, 0] = np.arange(periods + 1)  # we define the times
    M[0, 1] = notional
    for t in range(1, periods):
        M[t, 4] = rate * M[t - 1, 1]  # interest quote
        M[t, 3] = 0  # notional quote, 0 for bullet mortgage
        scheduled_oustanding = M[t - 1, 1] - M[t, 3]
        M[t, 2] = scheduled_oustanding * CPR[t]  # prepayment
        M[t, 1] = scheduled_oustanding - M[t, 2]  # notional(t) = notional(t-1) - (notional quote + prepayment)
        M[t, 5] = M[t, 4] + M[t, 2] + M[t, 3]

    M[periods, 4] = rate * M[periods - 1, 1]  # interest quote
    M[periods, 3] = M[periods - 1, 1]  # notional quote
    M[periods, 5] = M[periods, 4] + M[periods, 2] + M[periods, 3]
    return M


def Annuity(rate, notional, periods, CPR):
    # it returns a matrix M such that
    # M = [t  notional(t)  prepayment(t)  notional_quote(t)  interest_quote(t)  installment(t)]
    # WARNING! here "rate" and "periods" are quite general, the choice of getting year/month/day.. steps, depends on the rate
    # that the function receives. So, it is necessary to pass the correct rate to the function
    M = np.zeros((periods + 1, 6))
    M[:, 0] = np.arange(periods + 1)  # we define the times
    M[0, 1] = notional
    for t in range(1, periods + 1):
        # we are computing the installment at time t knowing the oustanding at time (t-1)
        remaining_periods = periods - (t - 1)

        # Installment, C(t_i)
        M[t, 5] = rate * M[t - 1, 1] / (1 - 1 / (1 + rate) ** remaining_periods)

        # Interest rate payment, I(t_i) = r * N(t_{i})
        M[t, 4] = rate * M[t - 1, 1]

        # Notional payment, Q(t_i) = C(t_i) - I(t_i)
        M[t, 3] = M[t, 5] - M[t, 4]

        # Prepayment, P(t_i)= Lambda * (N(t_i) -Q(t_i))
        M[t, 2] = CPR[t] * (M[t - 1, 1] - M[t, 3])

        # notional, N(t_{i+1}) = N(t_{i}) - lambda * (Q(t_{i} + P(t_i)))
        M[t, 1] = M[t - 1, 1] - M[t, 3] - M[t, 2]
    return M

def Bullet_price(rate, notonal, periods, CPR):
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
def Annuity_price(rate, notonal, periods, CPR):
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
    st.title("Modeling Mortgage ")
    st.write("This app is used to model the mortgage using different methods")
    st.write("The methods are:")


    st.sidebar.write("Select the method you want to use from the sidebar")
    st.sidebar.write("The methods are:")
    morgage_type = st.sidebar.selectbox("Select the method", ["Annuity", "Bullet", "Stochastic Amortizing Swap"])

    co1, co2 = st.columns(2)
    with co1:

        notional = st.sidebar.number_input("Enter the notional", min_value=1000, max_value=1000000000, value=1000000)
        r = st.sidebar.number_input("Enter the interest rate", 0.01, 1.0, 0.05)
    with co2:

        Lambda = st.sidebar.number_input("Enter the lambda (PrePayment rate) ", 0.0, 1.0, 0.01)
        T_end = st.sidebar.slider("Enter the end time", 1, 100, 30, 1)

    if morgage_type == "Annuity":
        M=Annuity_price(r, notional, T_end, Lambda)
        st.write(
            "An annuity mortgage is a type of home loan in which the borrower pays a fixed rate of interest on the principal amount for a specific period. The fixed rate of interest is determined at the time of the loan origination and remains the same throughout the loan term. The borrower pays a fixed amount of money every month to the lender, and this amount includes both the interest and the principal amount.")
        st.subheader("Annuity")
        plt.figure(1)
        plt.title("notonal of Annuity")
        plt.plot(M[:, 0], M[:, 1], '-r')
        plt.grid()
        plt.xlabel('time')
        plt.ylabel('notional')
        st.pyplot(plt)
        for i in range(0, T_end + 1):
            st.write(
            "Ti={0}, Notional={1:.0f}, Prepayment={2:.0f}, Notional Repayment={3:.0f}, Interest Rate={4:.0f}, Installment={5:.0f} ".format(
                M[i, 0], M[i, 1], M[i, 2], M[i, 3], M[i, 4], M[i, 5]))
    elif morgage_type == "Bullet":
        st.write(
            "A bullet  represents a structure where the borrower pays periodic interest only, and the entire principal is repaid as a lump sum (the 'bullet') at maturity. This setup is useful in analyzing cash flows and risk profiles in fixed-income and mortgage models.")

        M=Bullet_price(r, notional, T_end, Lambda)
        st.subheader("Bullet")

        plt.figure(1)
        plt.title("notonal of Bullets")
        plt.plot(M[:, 0], M[:, 1], '-r')
        plt.grid()
        plt.xlabel('time')
        plt.ylabel('notional')
        st.pyplot(plt)
        for i in range(0, T_end + 1):
            st.write(
                "Ti={0}, Notional={1:.0f}, Prepayment={2:.0f}, Notional Repayment={3:.0f}, Interest Rate={4:.0f}, Installment={5:.0f} ".format(
                    M[i, 0], M[i, 1], M[i, 2], M[i, 3], M[i, 4], M[i, 5]))

    elif morgage_type == "Stochastic Amortizing Swap":
        st.write(
            "An amortization swap is a derivative contract where one leg involves fixed interest payments while the other involves floating payments, with the notional amount reducing over time according to an amortization schedule. It helps in hedging interest rate risks on diminishing balances.")
        st.subheader("Stochastic Amortizing Swap")
        # incentive function
        Irrational = lambda x: 0.04 + 0.1 / (1 + np.exp(200 * (-x)))
        Rational = lambda x: 0.04 * (x > 0.0)
        incentive_functions = {"Irrational": Irrational, "Rational": Rational}
        Morgagse_types_selction={"Annuity":Annuity,"Bullet":Bullet}
        selete_morgage_type=st.sidebar.radio("Select the Morgage Type",list(Morgagse_types_selction.keys()))
        selected_key = st.sidebar.radio("Select the Incentive Function", list(incentive_functions.keys()))
        IncentiveFunction = incentive_functions[selected_key]
        K=r
        newRate = np.linspace(-0.1, 0.1, 150)
        NoOfPaths=2000
        NoOfSteps=T_end
        lambd = 0.05
        eta = 0.01
        Tend=T_end

        col1, col2 = st.columns(2)
        with col1:
            epsilon = K - newRate
            incentive = Irrational(epsilon)
            st.write('Irrational Incentive Function')
            plt.figure(1)
            plt.plot(newRate, incentive)
            plt.xlabel('S(t)')
            plt.ylabel('Incentive')
            plt.grid()
            st.pyplot(plt)
        with col2:
            epsilon = K - newRate
            incentive = Rational(epsilon)
            st.write('Rational Incentive Function')
            plt.figure(2)
            plt.plot(epsilon, incentive)
            plt.xlabel('epsilon= K - S(t)')
            plt.ylabel('Incentive')
            plt.grid()
            st.pyplot(plt)
        P0T = lambda T: np.exp(-0.05 * T)
        paths = GeneratePathsHWEuler(NoOfPaths, NoOfSteps, Tend, P0T, lambd, eta)
        R = paths["R"]

        tiGrid = paths["time"]
        # computing swap rate at point of assmuming
        S = np.zeros([NoOfPaths, NoOfSteps + 1])
        for (i, ti) in enumerate(tiGrid):
            S[:, i] = SwapRateHW(ti, ti, Tend + ti, 30, R[:, i], P0T, lambd, eta)





        #with col3:
        col3, col4 = st.columns(2)
        with col3:
            incentive = IncentiveFunction(epsilon)
            st.write("Incentive for prepayment given stochastis S(t)")
            plt.figure(3)
            plt.plot(epsilon, incentive, '.r')
            plt.xlabel('epsilon= K - S(t)')
            plt.ylabel('Incentive')
            plt.grid()
            plt.title('Incentive for prepayment given stochastis S(t)')
            st.pyplot(plt)
        with col4:
            st.write("Swap distribution at Tend")
            plt.figure(4)
            plt.hist(S[:, -1], bins=50)
            plt.grid()
            plt.title('Swap distribution at Tend')
            st.pyplot(plt)

        MortgageProfile = Morgagse_types_selction[selete_morgage_type]
        st.subheader("Mortgage Profile")
        N = np.zeros([NoOfPaths, NoOfSteps + 1])
        for i in range(0, NoOfPaths):
            epsilon = K - S[i, :]
            Lambda = IncentiveFunction(epsilon)
            NotionalProfile = MortgageProfile(K, notional, NoOfSteps, Lambda)
            N[i, :] = NotionalProfile[:, 1]

        plt.figure(6)
        plt.grid()
        plt.xlabel('time')
        plt.ylabel('notional')

        n = st.slider("Enter the number of paths", min_value=1, max_value=1000, value=100, step=1)
        for k in range(0, n):
            plt.plot(tiGrid, N[k, :], '-b')

        AnnuityProfile_NoPrepayment = MortgageProfile(K, notional, NoOfSteps, np.zeros(NoOfSteps + 1))
        plt.plot(tiGrid, AnnuityProfile_NoPrepayment[:, 1], '--r')
        plt.title('Notional profile')
        st.pyplot(plt)








if __name__ == "__main__":
    main()
