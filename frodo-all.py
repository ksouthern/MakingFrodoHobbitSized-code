from proba_util import *
from scipy.special import binom
import math
from math import floor, ceil, log
import copy
from tqdm import tqdm
import pickle
import sys

gauss_640 = {-11: 4.57763671875e-05, -10: 0.0001983642578125, -9: 0.0007171630859375, -8: 0.002197265625,
             -7: 0.005859375, -6: 0.013702392578125, -5: 0.02813720703125, -4: 0.0506744384765625,
             -3: 0.0800933837890625, -2: 0.111083984375, -1: 0.1351470947265625, 0: 0.144287109375,
             1: 0.1351470947265625, 2: 0.111083984375, 3: 0.0800933837890625, 4: 0.0506744384765625,
             5: 0.02813720703125, 6: 0.013702392578125, 7: 0.005859375, 8: 0.002197265625, 9: 0.0007171630859375,
             10: 0.0001983642578125, 11: 4.57763671875e-05}
gauss_976 = {-10: 1.52587890625e-05, -9: 9.1552734375e-05, -8: 0.0004425048828125, -7: 0.001800537109375,
             -6: 0.00604248046875, -5: 0.0167999267578125, -4: 0.0388336181640625, -3: 0.074493408203125,
             -2: 0.118621826171875, -1: 0.1568145751953125, 0: 0.172088623046875, 1: 0.1568145751953125,
             2: 0.118621826171875, 3: 0.074493408203125, 4: 0.0388336181640625, 5: 0.0167999267578125,
             6: 0.00604248046875, 7: 0.001800537109375, 8: 0.0004425048828125, 9: 9.1552734375e-05,
             10: 1.52587890625e-05}
gauss_1344 = {-6: 2 * (2 ** -16), -5: 40 * (2 ** -16), -4: 364 * (2 ** -16), -3: 2023 * (2 ** -16),
              -2: 6876 * (2 ** -16), -1: 14320 * (2 ** -16), 0: 18286 * (2 ** -16), 1: 14320 * (2 ** -16),
              2: 6876 * (2 ** -16), 3: 2023 * (2 ** -16), 4: 364 * (2 ** -16), 5: 40 * (2 ** -16), 6: 2 * (2 ** -16)}


def failure_rate_a(ps):
    # print(ps)
    chi_s = ps['chi']
    chi_e1 = ps['chi']
    chi_e2 = ps['chi']
    Ua = law_convolution(chi_e2, build_mod_switching_error_law(ps['q'], ps['p']))
    Ub = chi_e1
    R = law_product(chi_s, Ua)
    R = iter_law_convolution(R, ps['n'])
    S = law_product(chi_s, Ub)
    S = iter_law_convolution(S, ps['n'])
    S = law_convolution(R, S)
    return S


def failure_rate(ps, S):
    chi_e1 = ps['chi']
    chi_e2 = ps['chi']
    if ps['q'] == ps['t']:
        V = chi_e2
    else:
        V = law_convolution(chi_e2, build_mod_switching_error_law(ps['q'], ps['t']))
    F = law_convolution(S, V)
    return F


def prob_e_from_n(d_bit, e, m):
    return binom(m, e) * (d_bit ** e) * ((1 - d_bit) ** (m - e))


def d_ct_from_d_bit(d_bit, e, m):
    d_ct = 0
    for i in range(e + 1, m + 1):
        d_ct += prob_e_from_n(d_bit, i, m)
    return d_ct


def d_ct_from_f(f, q, e, m, b, thresh):
    d_bit = 0
    for i in f:
        if min(i, q - i) > int(thresh) or min(i, q - i) < int(-thresh):
            d_bit += f[i]
    d_ct = d_ct_from_d_bit(d_bit, e, math.ceil(m / b))
    return d_ct


def get_size(ps):
    c1 = math.log(ps['p'], 2) * ps['n'] * ps['m_bar']
    c2 = math.log(ps['t'], 2) * ps['n_bar'] * ps['m_bar']
    c3 = math.log(ps['t'], 2) * ps['w']
    return [math.ceil((c1 + c2) / 8), math.ceil((c1 + c3) / 8)]


def minimize_sum(n):
    """Finds a pair of integers so that their product is at least n and the sum is minimal.

    Args:
       n: An integer.

    Returns:
      A pair of integers.
    """
    opt_k1, opt_k2 = 1, n
    k1 = int(floor(sqrt(n))) + 1
    while k1 > 0 and k1 + n * 1. / k1 < opt_k1 + opt_k2:
        k2 = int(ceil(n * 1. / k1))
        if k1 + k2 < opt_k1 + opt_k2:
            opt_k1, opt_k2 = k1, k2
        k1 -= 1
    return opt_k1, opt_k2


pairs = [
    [0, 0],
    [1, 9],
    [2, 18],
    [3, 27],
    [4, 511 - 475],
    [5, 511 - 466],
    [6, 511 - 457],
    [7, 511 - 448],
    [8, 511 - 439],
    [9, 511 - 430],
    [10, 511 - 421],
    [11, 511 - 412],
    [12, 511 - 403],
    [13, 511 - 394],
    [14, 511 - 385],
    [15, 511 - 376],
    [16, 511 - 367],
    [18, 511 - 358],
    [19, 511 - 349],
    [20, 511 - 340],
    [21, 511 - 331],
    [22, 511 - 322],
    [23, 511 - 313],
    [25, 511 - 304],
    [26, 511 - 295],
    [27, 511 - 286],
    [28, 511 - 277],
    [29, 511 - 268],
    [30, 511 - 259]
]


def frodo_alts(target_size, target_d_ct, ps, mb=False, name="none"):
    resut_pairs = []
    lens = []
    fails = []
    sets = []
    best_ps = ps
    best_size = target_size
    mb_range = [1, 2, 3, 4, 5, 6, 7, 8]
    if not mb:
        mb_range = [ps['b']]
    for c in (range(10, int(log(ps['q'], 2)) + 1)):
        print("q = 2**" + str(c))
        q = 2 ** c
        for b in range(c - 6, c + 1):
            p = 2 ** b
            print("p = 2**" + str(b))
            test_ps = {'n': ps['n'], 'q': q, 'p': p, 'chi': ps['chi']}
            f_a = failure_rate_a(test_ps)
            for t in [2 ** a for a in range(1, b + 1)]:
                test_ps['t'] = t
                print("t = 2**" + str(int(log(t, 2))))
                f = failure_rate(test_ps, f_a)
                for mb_v in mb_range:
                    test_ps['b'] = mb_v
                    for bch in pairs:
                        test_ps['e'] = bch[0]
                        test_ps['m'] = ps['m'] + bch[1]
                        test_ps['w'] = math.ceil(test_ps['m'] / mb_v)
                        (n_bar, m_bar) = minimize_sum(math.ceil(test_ps['m'] / mb_v))
                        test_ps['n_bar'] = n_bar
                        test_ps['m_bar'] = m_bar
                        test_d_ct = d_ct_from_f(f, test_ps['q'], test_ps['e'], test_ps['m'], test_ps['b'],
                                                int(test_ps['q'] / (2 * (2 ** test_ps['b']))))
                        test_size = get_size(test_ps)
                        lens.append(test_size[0])
                        lens.append(test_size[1])
                        fails.append(test_d_ct)
                        fails.append(test_d_ct)
                        resut_pairs.append([test_size[0], test_d_ct])
                        resut_pairs.append([test_size[1], test_d_ct])
                        sets.append(copy.copy(test_ps))
                        # sets.append(copy.copy(test_ps))
                        # print(sets)
                        # print(lens)
                        # print(fails)
                        # print(resut_pairs)
    with open('results' + name + '.pkl', 'wb') as file:
        pickle.dump([lens, fails, resut_pairs, sets], file)


def main(run_time):
    if run_time == 1:  # 640
        frodo_640_ps = {'n': 640, 'q': 2 ** 15, 't': 2 ** 15, 'chi': gauss_640, 'n_bar': 8, 'm_bar': 8, 'p': 2 ** 15,
                        'e': 0, 'm': 128, 'b': 2, 'w': 64}
        # frodo_640_d_ct = d_ct_from_f(failure_rate(frodo_640_ps),frodo_640_ps['q'],frodo_640_ps['e'],frodo_640_ps['m'],frodo_640_ps['b'],int(2**15/2**3))
        frodo_640_size = get_size(frodo_640_ps)
        print("frodo_640")
        print("CT: " + str(frodo_640_size))
        # print("d_ct: " + str(math.log(frodo_640_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_640_size, 2 ** -5, frodo_640_ps, True, "640M=128")

    elif run_time == 2:  # 976
        frodo_976_ps = {'n': 976, 'q': 2 ** 16, 't': 2 ** 16, 'chi': gauss_976, 'n_bar': 8, 'm_bar': 8, 'p': 2 ** 16,
                        'e': 0, 'm': 192, 'b': 3, 'w': 64}
        # frodo_976_d_ct = d_ct_from_f(failure_rate(frodo_976_ps),frodo_976_ps['q'],frodo_976_ps['e'],frodo_976_ps['m'],frodo_976_ps['b'],int(2**16/2**4))
        frodo_976_size = get_size(frodo_976_ps)
        print("frodo_976")
        print("CT: " + str(frodo_976_size))
        # print("d_ct: " + str(math.log(frodo_976_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_976_size, 2 ** -5, frodo_976_ps, True, "976M=192")

    elif run_time == 3:  # 1344
        frodo_1344_ps = {'n': 1344, 'q': 2 ** 16, 't': 2 ** 16, 'chi': gauss_1344, 'n_bar': 8, 'm_bar': 8, 'p': 2 ** 16,
                         'e': 0, 'm': 256, 'b': 4, 'w': 64}
        # frodo_1344_d_ct = d_ct_from_f(failure_rate(frodo_1344_ps),frodo_1344_ps['q'],frodo_1344_ps['e'],frodo_1344_ps['m'],frodo_1344_ps['b'],int(2**16/2**5))
        frodo_1344_size = get_size(frodo_1344_ps)
        print("frodo_1344")
        print("CT: " + str(frodo_1344_size))
        # print("d_ct: " + str(math.log(frodo_1344_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_1344_size, 2 ** -5, frodo_1344_ps, True, "1344M=256")

    elif run_time == 4:  # 640, 256bit pt
        frodo_640_ps = {'n': 640, 'q': 2 ** 15, 'chi': gauss_640, 'n_bar': 11, 'm_bar': 12, 'p': 2 ** 15, 't': 2 ** 15,
                        'e': 0, 'm': 256, 'b': 2, 'w': 128}
        # frodo_640_d_ct = d_ct_from_f(failure_rate(frodo_640_ps),frodo_640_ps['q'],frodo_640_ps['e'],frodo_640_ps['m'],frodo_640_ps['b'],int(2**15/2**3))
        frodo_640_size = get_size(frodo_640_ps)
        print("frodo_640")
        print("CT: " + str(frodo_640_size))
        # print("d_ct: " + str(math.log(frodo_640_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_640_size, 2 ** -5, frodo_640_ps, True, "640M=256")

    elif run_time == 5:  # 976, 128bit pt
        frodo_976_ps = {'n': 976, 'q': 2 ** 16, 'chi': gauss_976, 'n_bar': 7, 'm_bar': 7, 'p': 2 ** 16, 't': 2 ** 15,
                        'e': 0, 'm': 128, 'b': 3, 'w': 43}
        # frodo_976_d_ct = d_ct_from_f(failure_rate(frodo_976_ps),frodo_976_ps['q'],frodo_976_ps['e'],frodo_976_ps['m'],frodo_976_ps['b'],int(2**16/2**4))
        frodo_976_size = get_size(frodo_976_ps)
        print("frodo_976")
        print("CT: " + str(frodo_976_size))
        # print("d_ct: " + str(math.log(frodo_976_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_976_size, 2 ** -5, frodo_976_ps, True, "976M=128")

    elif run_time == 6:  # 976, 256bit pt
        frodo_976_ps = {'n': 976, 'q': 2 ** 16, 'chi': gauss_976, 'n_bar': 9, 'm_bar': 10, 'p': 2 ** 16, 't': 2 ** 15,
                        'e': 0, 'm': 256, 'b': 3, 'w': 86}
        # frodo_976_d_ct = d_ct_from_f(failure_rate(frodo_976_ps),frodo_976_ps['q'],frodo_976_ps['e'],frodo_976_ps['m'],frodo_976_ps['b'],int(2**16/2**4))
        frodo_976_size = get_size(frodo_976_ps)
        print("frodo_976")
        print("CT: " + str(frodo_976_size))
        # print("d_ct: " + str(math.log(frodo_976_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_976_size, 2 ** -5, frodo_976_ps, True, "976M=256")

    elif run_time == 7:  # 1344, 128 bit pt
        frodo_1344_ps = {'n': 1344, 'q': 2 ** 16, 'chi': gauss_1344, 'n_bar': 6, 'm_bar': 6, 'p': 2 ** 16, 't': 2 ** 15,
                         'e': 0, 'm': 128, 'b': 4, 'w': 32}
        # frodo_1344_d_ct = d_ct_from_f(failure_rate(frodo_1344_ps),frodo_1344_ps['q'],frodo_1344_ps['e'],frodo_1344_ps['m'],frodo_1344_ps['b'],int(2**16/2**5))
        frodo_1344_size = get_size(frodo_1344_ps)
        print("frodo_1344")
        print("CT: " + str(frodo_1344_size))
        # print("d_ct: " + str(math.log(frodo_1344_d_ct,2)))
        print("###########################################")
        frodo_alts(frodo_1344_size, 2 ** -5, frodo_1344_ps, True, "1344M=128")
    quit()


if __name__ == '__main__':
    i = (int(sys.argv[1]))
    main(i)
