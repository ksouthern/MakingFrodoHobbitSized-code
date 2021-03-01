from proba_util import *
from scipy.special import binom
import math

def failure_rate(param_set):
    (n,k,q,chi1,chi2,du,dv,e,m,mb_v) = param_set
    chi_s = build_centered_binomial_law(chi1)
    chi_e1 = build_centered_binomial_law(chi1)
    chi_e2 = build_centered_binomial_law(chi2)
    Ua = build_mod_switching_error_law(q,2**du)
    Ub = build_mod_switching_error_law(q,2**dv)
    Rs = law_convolution(chi_e2,Ua)
    Re = law_convolution(chi_e2,Ub)
    B1 = law_product(chi_s,chi_e1)
    B2 = law_product(chi_s,Rs)
    V1 = iter_law_convolution(B1,n*k)
    V2 = iter_law_convolution(B2,n*k)
    V = law_convolution(V1,V2)
    F = law_convolution(V,Re)
    return F

def prob_e_from_n(d_bit,e,m):
    return binom(m,e)*(d_bit**e)*((1-d_bit)**(m-e))

def d_ct_from_d_bit(d_bit,e,m):
    d_ct = 0
    for i in range(e+1,m+1):
        d_ct += prob_e_from_n(d_bit,i,m)
    return d_ct

def d_ct_from_f(f,e,m,thresh):
    d_bit = 0
    for i in f:
        if abs(i) > int(thresh):
            d_bit += f[i]
    d_ct = d_ct_from_d_bit(d_bit,e,m)
    return d_ct

def get_size(ps):
    (n,k,q,chi1,chi2,du,dv,e,m,mb_v) = ps
    return n*k*du + m*dv

pairs=[
[0,0],
[1,9],
[2,18],
[3,27],
[4,511-475],
[5,511-466],
[6,511-457],
[7,511-448],
[8,511-439],
[9,511-430],
[10,511-421],
[11,511-412],
[12,511-403],
[13,511-394],
[14,511-385],
[15,511-376],
[16,511-367],
[18,511-358],
[19,511-349],
[20,511-340],
[21,511-331],
[22,511-322],
[23,511-313],
[25,511-304],
[26,511-295],
[27,511-286],
[28,511-277],
[29,511-268],
[30,511-259]
] # pair[a,b] corrects a errors with redundancy b


def kyber_alts(target_size,target_d_ct,starting_ps,mb = False):
    (n,k,q,chi1,chi2,du,dv,e,m,mb_v)=(starting_ps)
    best_ps = ()
    best_size = target_size
    if mb:
        mb_range = [2,3,4]
    else:
        mb_range = [1]
    for mb_v in mb_range:
        for (du,dv) in [(12,4),(12,3),(12,2),(11,4),(11,3),(11,2),(10,4),(10,3),(10,2),(9,5),(9,4),(9,3),(9,2),(8,5),(8,4),(8,3),(8,2),(7,5),(7,4),(7,3),(7,2),(6,4),(6,3),(6,2)]:
            test_ps = (n,k,q,chi1,chi2,du,dv,e,math.ceil(m/mb_v),mb_v)
            f = failure_rate(test_ps)
            for bch in pairs:
                e = bch[0]
                m = n+bch[1]
                test_ps = (n,k,q,chi1,chi2,du,dv,e,math.ceil(m/mb_v),mb_v)
                test_d_ct = d_ct_from_f(f,e,math.ceil(m/mb_v),int(3329/(2*(2**mb_v))))
                if test_d_ct < target_d_ct:
                    test_size = get_size(test_ps)
                    if test_size < best_size:
                        print(str(test_size))
                        best_size = test_size
                        best_ps = test_ps
    (n,k,q,chi1,chi2,du,dv,e,m,mb_v) = best_ps
    print("###########################################")
    print("n: " + str(n))
    print("k: " + str(k))
    print("q: " + str(q))
    print("chi_1 :" + str(chi1))
    print("chi_2 :" + str(chi2))
    print("du: " + str(du))
    print("dv: " + str(dv))
    print("e: " + str(e))
    print("m: " + str(m))
    print("mb: " + str(mb_v))
    print("ct: " + str(best_size))
    print("d_ct: " + str(math.log(d_ct_from_f(failure_rate(best_ps),e,m,int(3329/(2*(2**mb_v)))),2))) 
    print("###########################################")

#                n  k  q  chi1 chi2 du dv e  m mb_v
kyber_512_ps = (256,2,3329,3,   2,   10,4,0,256,1)
kyber_512_d_ct = d_ct_from_f(failure_rate(kyber_512_ps),0,256,int(3329/4))
kyber_512_size = get_size(kyber_512_ps)
print("KYBER-512")
print("CT: " + str(kyber_512_size))
print("d_ct: " + str(math.log(kyber_512_d_ct,2)))
print("###########################################")
kyber_alts(kyber_512_size,kyber_512_d_ct,kyber_512_ps)
kyber_alts(kyber_512_size,kyber_512_d_ct,kyber_512_ps,True)

kyber_768_ps = (256,3,3329,2,2,10,4,0,256,1)
kyber_768_d_ct = d_ct_from_f(failure_rate(kyber_768_ps),0,256,int(3329/4))
kyber_768_size = get_size(kyber_768_ps)
print("KYBER-768")
print("CT: " + str(kyber_768_size))
print("d_ct: " + str(math.log(kyber_768_d_ct,2)))
print("###########################################")
kyber_alts(kyber_768_size,kyber_768_d_ct,kyber_768_ps)
kyber_alts(kyber_768_size,kyber_768_d_ct,kyber_768_ps,True)

kyber_1024_ps = (256,4,3329,2,2,11,5,0,256,1)
kyber_1024_d_ct = d_ct_from_f(failure_rate(kyber_1024_ps),0,256,int(3329/4))
kyber_1024_size = get_size(kyber_1024_ps)
print("KYBER-1024")
print("CT: " + str(kyber_1024_size))
print("d_ct: " + str(math.log(kyber_1024_d_ct,2)))
print("###########################################")
kyber_alts(kyber_1024_size,kyber_1024_d_ct,kyber_1024_ps)
kyber_alts(kyber_1024_size,kyber_1024_d_ct,kyber_1024_ps,True)


