from proba_util import *
from scipy.special import binom
import math
from math import floor, ceil, log
import pickle
from tqdm import tqdm

def failures(ps):
    #functions
    def failure_rate_a(ps):
        #print(ps)
        if not 'k' in ps:
            ps['k'] = 1
        if not 'LWR' in ps:
            ps['LWR'] = False
        chi_s = ps['chis']
        chi_e1 = ps['chie1']
        chi_e2 = ps['chie2']
        Ua = law_convolution(chi_e2,build_mod_switching_error_law(ps['q'],ps['p']))
        if ps['LWR']:
            Ub = law_convolution(chi_e1,build_mod_switching_error_law(ps['q'],ps['p']))
        else:
            Ub = chi_e1
        R = law_product(chi_s,Ua)
        R = iter_law_convolution(R,ps['n']*ps['k'])
        S = law_product(chi_s,Ub)
        S = iter_law_convolution(S,ps['n']*ps['k'])
        S = law_convolution(R,S)
        return S

    def failure_rate(ps,S):
        chi_e1 = ps['chie1']
        chi_e2 = ps['chie2']
        if ps['q'] == ps['t']:
            V = chi_e2
        else:
            V = law_convolution(chi_e2,build_mod_switching_error_law(ps['q'],ps['t']))
        F = law_convolution(S,V)
        return F

    def prob_e_from_n(d_bit,e,m):
        return binom(m,e)*(d_bit**e)*((1-d_bit)**(m-e))

    def d_ct_from_d_bit(d_bit,e,m):
        d_ct = 0
        for i in range(e+1,m+1):
            d_ct += prob_e_from_n(d_bit,i,m)
        return d_ct

    def d_ct_from_f(f,q,e,m,b,thresh):
        d_bit = 0
        for i in f:
            if min(i,q-i) > int(thresh) or min(i,q-i) < int(-thresh):
                d_bit += f[i]
        d_ct = d_ct_from_d_bit(d_bit,e,math.ceil(m/b))
        return d_ct
    #finding failure rate
    print("Starting "+ps['name'])
    f_a = failure_rate_a(ps)
    f = failure_rate(ps,f_a)
    d_ct = d_ct_from_f(f,ps['q'],ps['e'],ps['m'],ps['b'],int(ps['q']/(2*(2**ps['b']))))
    print(ps['name'] +": " + str(log(d_ct,2)))                        

def test_ps():
    #Frodo ps
    #-148.81718112355094
    #-199.55708684917505
    #-252.49239945959872
    gauss_640 = {-11: 4.57763671875e-05, -10: 0.0001983642578125, -9: 0.0007171630859375, -8: 0.002197265625, -7: 0.005859375, -6: 0.013702392578125, -5: 0.02813720703125, -4: 0.0506744384765625, -3: 0.0800933837890625, -2: 0.111083984375, -1: 0.1351470947265625, 0: 0.144287109375, 1: 0.1351470947265625, 2: 0.111083984375, 3: 0.0800933837890625, 4: 0.0506744384765625, 5: 0.02813720703125, 6: 0.013702392578125, 7: 0.005859375, 8: 0.002197265625, 9: 0.0007171630859375, 10: 0.0001983642578125, 11: 4.57763671875e-05}
    gauss_976 = {-10: 1.52587890625e-05, -9: 9.1552734375e-05, -8: 0.0004425048828125, -7: 0.001800537109375, -6: 0.00604248046875, -5: 0.0167999267578125, -4: 0.0388336181640625, -3: 0.074493408203125, -2: 0.118621826171875, -1: 0.1568145751953125, 0: 0.172088623046875, 1: 0.1568145751953125, 2: 0.118621826171875, 3: 0.074493408203125, 4: 0.0388336181640625, 5: 0.0167999267578125, 6: 0.00604248046875, 7: 0.001800537109375, 8: 0.0004425048828125, 9: 9.1552734375e-05, 10: 1.52587890625e-05}
    gauss_1344 = {-6:2*(2**-16), -5:40*(2**-16), -4:364*(2**-16), -3:2023*(2**-16), -2:6876*(2**-16), -1:14320*(2**-16), 0:18286*(2**-16),1:14320*(2**-16),2:6876*(2**-16),3:2023*(2**-16),4:364*(2**-16),5:40*(2**-16),6:2*(2**-16)}

    frodo_640_ps =  {'name':"Frodo640", 'expected':-148.81718112355094,'n':640, 'q':2**15,'t':2**15,'chis':gauss_640,'chie1':gauss_640,'chie2':gauss_640,'n_bar':8,'m_bar':8,'p':2**15,'e':0,'m':128,'b':2, 'w':64}
    frodo_976_ps =  {'name':"Frodo976", 'expected':-199.55708684917505,'n':976, 'q':2**16,'t':2**16,'chis':gauss_976,'chie1':gauss_976,'chie2':gauss_976,'n_bar':8,'m_bar':8,'p':2**16,'e':0,'m':192,'b':3, 'w':64}
    frodo_1344_ps = {'name':"Frodo1344", 'expected':-252.49239945959872,'n':1344,'q':2**16,'t':2**16,'chis':gauss_1344,'chie1':gauss_1344,'chie2':gauss_1344,'n_bar':8,'m_bar':8,'p':2**16,'e':0,'m':256,'b':4, 'w':64}

    #Kyber ps
    kyber512_ps =  {'name':"Kyber512", 'expected':-139, 'n':256, 'k':2, 'q':3329, 'p':2**10, 't':2**4, 'chis':build_centered_binomial_law(3),'chie1':build_centered_binomial_law(3), 'chie2':build_centered_binomial_law(2),'e':0,'b':1,'m':256}
    kyber768_ps =  {'name':"Kyber768", 'expected':-164, 'n':256, 'k':3, 'q':3329, 'p':2**10, 't':2**4, 'chis':build_centered_binomial_law(2),'chie1':build_centered_binomial_law(2), 'chie2':build_centered_binomial_law(2),'e':0,'b':1,'m':256}
    kyber1024_ps = {'name':"Kyber1024",'expected':-174, 'n':256, 'k':4, 'q':3329, 'p':2**11, 't':2**5, 'chis':build_centered_binomial_law(2),'chie1':build_centered_binomial_law(2), 'chie2':build_centered_binomial_law(2),'e':0,'b':1,'m':256}

    #Saber ps
    LightSaber_ps = {'name':"LightSaber", 'expected':-120, 'n':256, 'k':2, 'q':2**13, 'p':2**10, 't':2**3, 'chis':build_centered_binomial_law(5),'chie1':{0:1}, 'chie2':{0:1},'e':0,'b':1,'m':256,'LWR':True}
    Saber_ps =      {'name':"Saber",      'expected':-136, 'n':256, 'k':3, 'q':2**13, 'p':2**10, 't':2**4, 'chis':build_centered_binomial_law(4),'chie1':{0:1}, 'chie2':{0:1},'e':0,'b':1,'m':256,'LWR':True}
    FireSaber_ps =  {'name':"FireSaber",  'expected':-165, 'n':256, 'k':4, 'q':2**13, 'p':2**10, 't':2**6, 'chis':build_centered_binomial_law(3),'chie1':{0:1}, 'chie2':{0:1},'e':0,'b':1,'m':256,'LWR':True}


    param_sets = [frodo_640_ps,frodo_976_ps,frodo_1344_ps,kyber512_ps,kyber768_ps,kyber1024_ps,LightSaber_ps,Saber_ps,FireSaber_ps]
    for i in param_sets:
        failures(i)
        print("Expected: " + str(i['expected']))

def check_correctness(lens,fails,sets):
    listlen = len(lens)
    if len(lens)!= len(fails) or len(fails)!=len(sets):
        print("len error")
    #lens
    #fails
    #sets
    for i in tqdm(range(listlen)):
        for j in range(i,listlen):
            if sets[i]['e'] == sets[j]['e']  and sets[i]['b'] == sets[j]['b'] and sets[i]['m'] == sets[j]['m'] and sets[i]['q'] >= sets[j]['q'] and sets[i]['p'] >= sets[j]['p'] and sets[i]['t'] >= sets[j]['t']:
                if fails[i]< fails[j] or lens[i]<lens[j]:
                    print("FALSE")
    
        
def check_file(name):
    lists = pickle.load(open(name+".pkl","rb"))
    lens = lists[0]
    fails = lists[1]
    result_pairs = lists[2]
    lens1 = [lens[i] for i in range(0,len(lens)) if i%2 == 0]
    fails1 = []
    sets = lists[3]
    for i in range(0,len(fails)):
        if i%2 == 0:
            if fails[i] != fails[i+1]:
                print("ERROR")
            if fails[i] == 0:
                fails1.append(1500)
            else:
                fails1.append(-log(fails[i],2))
    check_correctness(lens1,fails1,sets)

def check_files():
    files = ["results640M=128","results976M=192","results1344M=256","results640M=256","results976M=128","results976M=256","results1344M=128"]

    for i in files:
        print("Starting file: " + i)
        check_file(i)
        print("Finished file: " + i)

                
test_ps()        
    
