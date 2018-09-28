RR = RealField(128)

# |\vartheta(x) - x| <= eta*x/log^k0(x) for log(x) >= x_log_lower_bound
# Dusart
eta = 0.01
x_log_lower_bound = log(7713133853)
k0 = 2

LOG_LOG_BRIGGS_HEIGHT = 13.099

def prime_upper_bound(n):
    if n >= 688383:
        return n*(log(n) +log(log(n)) - 1 + (log(log(n)) - 2)/log(n))
    else:
        return +Infinity

def prime_lower_bound(n):
    if n >= 3:
        return n*(log(n) +log(log(n)) - 1 + (log(log(n)) - 2.1)/log(n))
    else:
        return 2

RIEMANN_HEIGHT = 30610046000

def Buethe_check(x):
    return bool(0.492 * sqrt(x/log(x)) <= RIEMANN_HEIGHT)

def primorial_log_lower_bound(n):
    nth_prime = prime_lower_bound(n)
    if nth_prime >= 10^15:
        Dusart = n*(log(n) + log(log(n)) -1 + (log(log(n))-2.04)/log(n))
    elif nth_prime >= 10^11:
        Dusart = n*(log(n) + log(log(n)) -1 + (log(log(n))-2.050735)/log(n))
    else:
        Dusart = n*(log(n) + log(log(n)) -1 + (log(log(n))-2.1454)/log(n))
    if Buethe_check(nth_prime):
        Buethe = (nth_prime
                  - 0.125*sqrt(nth_prime)*log(nth_prime)*(log(nth_prime)-2)/pi)
        return max(Buethe, Dusart)
    else:
        return Dusart

def primorial_log_upper_bound(n):
    nth_prime = prime_upper_bound(n)
    Dusart = n*(log(n) + log(log(n)) - 1 + (log(log(n))-2)/log(n))
    if Buethe_check(nth_prime):
        Buethe = (nth_prime
                  + 0.125*sqrt(nth_prime)*log(nth_prime)*(log(nth_prime)-2)/pi)
        return min(Buethe, Dusart)
    else:
        return Dusart

def primorial_log_lower_bound(n):
    nth_prime = prime_lower_bound(n)
    Dusart = n*(log(n) + log(log(n)) - 1 + (log(log(n))-2.1454)/log(n))
    if Buethe_check(nth_prime):
        Buethe = (nth_prime
                  - 0.125*sqrt(nth_prime)*log(nth_prime)*(log(nth_prime)-2)/pi)
        return max(Buethe, Dusart)
    else:
        return Dusart

def prod_upper_bound(eta, x_log_lower_bound, k):
#   Calculate constant in Lemma 5
    return exp(eta * (1/k + (k+2)/((k+1)*x_log_lower_bound)))-1

def Robin_upper_bound(n,t):
#   compute g_t(n)
    nth_prime = prime_upper_bound(n)
    prod_bound = prod_upper_bound(eta, x_log_lower_bound, k0)
    numerator = (exp(2/prime_lower_bound(n))*log(nth_prime)
                 *(1 + prod_bound/(log(prime_lower_bound(n))^k0)))
    denominator = zeta(t)*log(primorial_log_lower_bound(n))
    return numerator/denominator

def min_t(n):
#   compute t so that R_t(N_n) < exp(gamma) using upper bound of R_t
    nth_prime = prime_upper_bound(n)
    prod_bound = prod_upper_bound(eta, x_log_lower_bound, k0)
    numerator = (exp(2/prime_lower_bound(n))*log(nth_prime)
                 *(1 + prod_bound/log(prime_lower_bound(n))^k0))
    denominator = log(primorial_log_lower_bound(n))
    t = 2
    while numerator < zeta(t)*denominator:
        t += 1
    t -= 1
    print("At worst, Robin's inequality holds for {}-free integers larger than "
          "N_{}.".format(t, n))

def bin_search(t,N,s):
    Lwr = N*2^(s-1)
    Upr = N*2^s
    Mdl = (Lwr + Upr)/2
    while RR(Upr - Lwr) > 1:
        if RR(Robin_upper_bound(Mdl, t)) < 1:
            Upr = Mdl
        else:
            Lwr = Mdl
        Mdl = floor((Lwr + Upr)/2)
    if RR(Robin_upper_bound(Upr, t)) < 1:
        return Upr
    elif RR(Robin_upper_bound(Lwr, t)) < 1:
        return Lwr
    else:
        return "Rounding error?"

def double_check(k, t):
    if RR(Robin_upper_bound(k, t)) >= 1:
        return "Error: k is too small."
    i = 0
    while RR(Robin_upper_bound(k, t)) < 1:
        k = k - 1
        i = i + 1
        if i == 1000:
            print "We're not getting any younger."
            return k
    return k + 1

def eta_check(n):
    n_log_lower_bound = log(prime_lower_bound(n))
    if n_log_lower_bound >= x_log_lower_bound:
        print "We got \eta_{} to kick in at the right time.".format(k0)
    else:
        print "We need \eta_{} to kick in at b = {}.".format(k0,
                                                    round(n_log_lower_bound, 3))

def t_free(t):
    if RIEMANN_HEIGHT == +Infinity:
        print "Assuming the Riemann Hypothesis,"
#   Lower bound for eta to kick in
    N = 355231136
    s = 1
    while RR(Robin_upper_bound(N*2^s, t)) >= 1:
        s = s+1
    n = bin_search(t,N,s)
    n = double_check(n, t)
    eta_check(n)
    log_log_height = round(log(primorial_log_upper_bound(n)/log(10),10), 3)
    if log_log_height <= LOG_LOG_BRIGGS_HEIGHT:
        string = "and N_n is about 10^10^{}".format(log_log_height)
    else:
        string1 = "but N_n is too big. Can we get the 'Briggs height' to "
        string2 = "10^10^{}?".format(log_log_height)
        string = string1 + string2
    print("The {}th primorial is sufficient for {}-free RI, with "
          "R(N_n)/exp(gamma) < {},".format(n, t, RR(Robin_upper_bound(n, t))),
                                           string)
