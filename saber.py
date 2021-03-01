from proba_util import *
from scipy.special import binom
import math

def failure_rate(param_set):
    (n,k,q,chi,p,t,e,m,mb_v) = param_set
    chi_s = build_centered_binomial_law(chi)
    Ua = build_mod_switching_error_law(q,2**p)
    Ub = build_mod_switching_error_law(q,2**t)
    B = law_product(chi_s,Ua)
    V = iter_law_convolution(B,2*n*k)
    F = law_convolution(V,Ub)
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
    (n,k,q,chi,p,t,e,m,mb_v) = ps
    return n*k*p + m*t

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


def saber_alts(target_size,target_d_ct,starting_ps,mb = False):
    (n,k,q,chi,p,t,e,m,mb_v)=(starting_ps)
    best_ps = ()
    best_size = target_size
    if mb:
        mb_range = [2,3,4]
    else:
        mb_range = [1]
    for mb_v in mb_range:
        seta = [12,11,10,9,8,7,6]
        setb = [7,6,5,4,3,2]
        sets = [(a,b) for a in seta for b in setb if a > b and ((mb_v>1) or (a<p))]
        for (p,t) in sets:
            test_ps = (n,k,q,chi,p,t,e,math.ceil(m/mb_v),mb_v)
            test_size = get_size(test_ps)
            #print("SET: " + str(p)+","+str(t))
            #print("SIZE: " + str(math.ceil(test_size/8)))
            if test_size < best_size:
                f = failure_rate(test_ps)
                for bch in pairs:
                    e = bch[0]
                    m = n+bch[1]
                    test_ps = (n,k,q,chi,p,t,e,math.ceil(m/mb_v),mb_v)
                    test_d_ct = d_ct_from_f(f,e,math.ceil(m/mb_v),int(q/(2*(2**mb_v))))
                    if test_d_ct < target_d_ct:
                        test_size = get_size(test_ps)
                        if test_size < best_size:
                            print(str(test_size))
                            best_size = test_size
                            best_ps = test_ps
                        else:
                            break
    (n,k,q,chi,du,dv,e,m,mb_v) = best_ps
    print("###########################################")
    print("n: " + str(n))
    print("k: " + str(k))
    print("q: " + str(q))
    print("chi :" + str(chi))
    print("p: " + str(p))
    print("t: " + str(t))
    print("e: " + str(e))
    print("m: " + str(m))
    print("mb: " + str(mb_v))
    print("ct: " + str(best_size))
    print("d_ct: " + str(math.log(d_ct_from_f(failure_rate(best_ps),e,m,int(2**13/(2*(2**mb_v)))),2))) 
    print("###########################################")

#                  n  k   q  chi p r e m mb_v
light_saber_ps = (256,2,2**13,5,10,3,0,256,1)
light_saber_d_ct = d_ct_from_f(failure_rate(light_saber_ps),0,256,int(2**13/4))
light_saber_size = get_size(light_saber_ps)
print("light_saber")
print("CT: " + str(light_saber_size))
print("d_ct: " + str(math.log(light_saber_d_ct,2)))
print("###########################################")
#saber_alts(light_saber_size,light_saber_d_ct,light_saber_ps)
#saber_alts(light_saber_size,light_saber_d_ct,light_saber_ps,True)

saber_ps = (256,3,2**13,4,10,4,0,256,1)
saber_d_ct = d_ct_from_f(failure_rate(saber_ps),0,256,int(2**13/4))
saber_size = get_size(saber_ps)
print("saber")
print("CT: " + str(saber_size))
print("d_ct: " + str(math.log(saber_d_ct,2)))
print("###########################################")
#saber_alts(saber_size,saber_d_ct,saber_ps)
#saber_alts(saber_size,saber_d_ct,saber_ps,True)

fire_saber_ps = (256,4,2**13,3,10,6,0,256,1)
fire_saber_d_ct = d_ct_from_f(failure_rate(fire_saber_ps),0,256,int(2**13/4))
fire_saber_size = get_size(fire_saber_ps)
print("fire_saber")
print("CT: " + str(fire_saber_size))
print("d_ct: " + str(math.log(fire_saber_d_ct,2)))
print("###########################################")
saber_alts(fire_saber_size,fire_saber_d_ct,fire_saber_ps)
saber_alts(fire_saber_size,fire_saber_d_ct,fire_saber_ps,True)


