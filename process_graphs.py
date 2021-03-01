import pickle
import matplotlib.pyplot as plt
import copy
from mpmath import *

mp.dps = 40
from math import ceil
import os
import sys
# from frodo import get_size
from tqdm import tqdm

# from estimateLWE.estimates import para_cost_scheme, flatten
# from schemes import LWE_SCHEMES
(v, _, _, _, _) = sys.version_info


# pkl_converter("results640M=256",4,2)
# pkl_converter("results976",4,2)
# pkl_converter("results1344",4,2)
def main(name):
    '''start of helper functions'''

    def pkl_converter(f, p1, p2):
        # converts a file pickled with p1 to one pickled with p2
        if not os.path.exists(f + "." + str(p1) + "pkl"):
            return False
        else:
            temp = pickle.load(open(f + "." + str(p1) + "pkl", "rb"))
            pickle.dump(temp, open(f + "." + str(p2) + "pkl", "wb"), protocol=p2)

    def dict_to_string(ps):
        string = ""
        string += ("n: " + str(ps['n']) + '\n')
        string += ("q: " + str(ps['q']) + '\n')
        string += ("chi :" + str(ps['chi']) + '\n')
        string += ("p: " + str(ps['p']) + '\n')
        string += ("t: " + str(ps['t']) + '\n')
        string += ("e: " + str(ps['e']) + '\n')
        string += ("m: " + str(ps['m']) + '\n')
        string += ("n_bar: " + str(ps['n_bar']) + '\n')
        string += ("m_bar: " + str(ps['m_bar']) + '\n')
        string += ("b: " + str(ps['b']) + '\n')
        string += ("ct: " + str(get_size(ps)[1]) + '\n')
        return string

    def dict_to_latex(ps):
        if "M=128" in name:
            pt = "128"
        elif "M=192" in name:
            pt = "192"
        elif "M=256" in name:
            pt = "256"
        else:
            pt = "pt"
        # name, n, pt, q, p, t, c, b, |ct|, |ct|/|pt| ,delta_ct
        string = "ADD NAME &"
        string += ("$" + str(ps['n']) + "$ &")
        string += ("$" + pt + "$ &")
        string += ("$ 2^{" + str(int(log(ps['q'], 2))) + "}$ &")
        string += ("$ 2^{" + str(int(log(ps['p'], 2))) + "}$ &")
        # if ps['n'] == 640:
        #    string += "$N_{12}(2.8)$&"
        #    string += "$N_{12}(2.8)$&"
        # elif ps['n'] == 976:
        #    string += "$N_{10}(2.3)$&"
        #    string += "$N_{10}(2.3)$&"
        # elif ps['n'] == 1344:
        #    string += "$N_{6}(1.4)$&"
        #    string += "$N_{6}(1.4)$&"
        string += ("$ 2^{" + str(int(log(ps['t'], 2))) + "}$ &")
        string += ("$" + str(ps['e']) + "$ &")
        string += ("$" + str(ps['b']) + "$ &")
        string += ("$ " + str(get_size(ps)[1]) + "$ &")
        string += ("$ " + str(ceil((8 * get_size(ps)[1]) / int(pt))) + "$ &")
        return string

    def get_size(ps):
        c1 = log(ps['p'], 2) * ps['n'] * ps['m_bar']
        c2 = log(ps['t'], 2) * ps['n_bar'] * ps['m_bar']
        c3 = log(ps['t'], 2) * ps['w']
        return [ceil((c1 + c2) / 8), ceil((c1 + c3) / 8)]

    def create_param(ps, q):
        if ps['n'] == 640:
            param = {
                "n": 640,
                "q": q,
                "sd": 2.75,
                "secret_distribution": "normal",
                "claimed": 103,
                "category": [1, ]}
        elif ps['n'] == 976:
            param = {
                "n": 976,
                "q": q,
                "sd": 2.3,
                "secret_distribution": "normal",
                "claimed": 150,
                "category": [3, ]}
        elif ps['n'] == 1344:
            param = {
                "n": 1344,
                "q": q,
                "sd": 1.4,
                "secret_distribution": "normal",
                "claimed": 197,
                "category": [5, ]}
        else:
            param = {}
            print("ERROR")
            return ("ERROR")

        return param

    def get_scheme(sets):
        s = {
            "name": "Frodo",
            "assumption": [
                "LWE",
            ],
            "primitive": [
                "KEM",
                "PKE",
            ],
            "params": [],
        }
        for i in range(8, 16):
            s["params"].append(create_param(sets[0], 2 ** i))
        return s

    def get_estimates(scheme):
        print(scheme)
        # scheme2 = [s for s in LWE_SCHEMES if s['name'] == 'Frodo']
        # print(scheme2)
        lwe_estimates = list(para_cost_scheme([[scheme]]))
        print(lwe_estimates)
        estimates_list = flatten([x[1] for x in lwe_estimates])
        primal = []
        dual = []
        c = 0
        for i in estimates_list:
            c += 1
            if c % 2 == 1:
                primal.append(i['cost']['Q\xe2\x80\x91Core\xe2\x80\x91Sieve']['n']['rop'])
            elif c % 2 == 0:
                dual.append(i['cost']['Q\xe2\x80\x91Core\xe2\x80\x91Sieve']['n']['rop'])
        return (primal, dual)

    def get_frodo_v(name):
        m_for_n = {"640": "128", "976": "192", "1344": "256"}
        [n] = [i for i in ["640", "976", "1344"] if i in name]
        [m] = [i for i in ["128", "192", "256"] if i in name]
        frodo_v = "FrodoKEM-" + n
        if m_for_n[n] != m:
            frodo_v += " with m = " + m
        return frodo_v

    def plotting():
        plotx = []
        ploty = []
        plotx1 = []
        ploty1 = []
        plotx2 = []
        ploty2 = []
        plotx3 = []
        ploty3 = []
        points = []
        points1 = []
        points2 = []
        points3 = []
        thresh = 50
        x = True
        for i in range(len(lens2)):
            # x = sets[i]['b'] == 4
            if save: x = fails1[i] < 1450
            if "640" in name:
                y = sets[i]['q'] <= 2 ** 15
            else:
                y = True
            # x = sets[i]['b'] > 4
            # x = sets[i]['q'] == 2**16 and sets[i]['p'] == 2**16
            if fails1[i] > thresh and x and y:
                ploty.append(lens2[i])
                plotx.append(fails1[i])
                points.append(sets[i])
        frodo_v = get_frodo_v(name)
        title = "Plot of varients of " + frodo_v + "."
        pltname = frodo_v + "-all-points.pdf"
        plot(plotx, ploty, points, title, pltname)

    def plot(x, y, data, title, pltname):
        if "640" in name:
            plt_og = [145, 9720]
        elif "976" in name:
            plt_og = [200, 15744]
        elif "1344" in name:
            plt_og = [253, 21632]
        else:
            plt_og = [0, 0]

        def update_annot(ind):
            x, y = line.get_data()
            annot.xy = (x[ind["ind"][0]], y[ind["ind"][0]])
            text = "{}, {}, {}".format(" ".join(list(map(str, ind["ind"]))),
                                       " ".join([dict_to_string(data[n]) for n in ind["ind"]]),
                                       " ".join([str(x[n]) for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)

        def point_data(ind):
            x, y = line.get_data()
            points = [dict_to_latex(data[n]) + "$2^{-" + str(int(x[n])) + "}$ \\ " for n in ind["ind"]]
            return points

        def hover(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = line.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

        def onclick(event):
            latex_str = ""
            if event.inaxes == ax:
                cont, ind = line.contains(event)
                latex_str = point_data(ind)
                for data in latex_str:
                    print(data)

        fig, ax = plt.subplots()
        line, = plt.plot(x, y, 'ro')
        plt.plot(plt_og[0], plt_og[1], 'bo')
        # plt.plot(plotx,ploty,'ro')
        annot = ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        fig.canvas.mpl_connect("motion_notify_event", hover)
        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.xlabel(r'$\delta_{ct}~(- \log_{2}(x))$')
        plt.ylabel(r'Ciphertext size (Bytes)')
        plt.title(title)
        if save:
            plt.savefig(pltname)
        else:
            plt.show()

    def remove_top(plotx, ploty, points):
        a = []
        b = []
        c = []
        for i in range(len(plotx)):
            if plotx[i] < 1500:
                a.append(plotx[i])
                b.append(ploty[i])
                c.append(points[i])
        return (a, b, c)

    def get_smallest(name):
        plotx = []
        ploty = []
        points = []
        if os.path.exists(name + "-smallest.pkl"):
            [ploty, plotx, points] = pickle.load(open(name + "-smallest.pkl", "rb"))
        else:
            for i in tqdm(range(len(lens2))):
                if "640" in name:
                    y = sets[i]['q'] <= 2 ** 15
                else:
                    y = True
                if fails1[i] > 50 and y:
                    use = True
                    for j in range(len(lens2)):
                        if "640" in name:
                            z = sets[j]['q'] <= 2 ** 15
                        else:
                            z = True
                        if (not i == j) and fails1[j] > 50 and z:
                            if lens2[j] < lens2[i] and fails1[j] >= fails1[i]:
                                use = False
                    if use:
                        ploty.append(lens2[i])
                        plotx.append(fails1[i])
                        points.append(sets[i])
            pickle.dump([ploty, plotx, points], open(name + "-smallest.pkl", "wb"))
        if save:
            (plotx, ploty, points) = remove_top(plotx, ploty, points)

        frodo_v = get_frodo_v(name)
        title = "Optimal sizes vs failure rate for varients of " + frodo_v + "."
        pltname = frodo_v + "-smallest.pdf"
        plot(plotx, ploty, points, title, pltname)

    ''' end of helper functions'''

    if v == 3:
        lists = pickle.load(open(name + ".pkl", "rb"))
    elif v == 2:
        lists = pickle.load(open(name + ".2pkl", "rb"))
    lens = lists[0]
    fails = lists[1]
    result_pairs = lists[2]
    lens1 = [lens[i] for i in range(0, len(lens)) if i % 2 == 0]
    lens2 = [lens[i] for i in range(0, len(lens)) if i % 2 == 1]
    fails1 = []
    fails2 = []
    sets = lists[3]
    for i in range(0, len(fails)):
        if i % 2 == 0:
            if fails[i] == 0:
                fails1.append(1500)
            else:
                fails1.append(-log(fails[i], 2))
        else:  # if i%2 == 0:
            if fails[i] == 0:
                fails2.append(1500)
            else:
                fails2.append(-log(fails[i], 2))

    result_pairs1 = [result_pairs[i] for i in range(0, len(result_pairs)) if i % 2 == 0]
    result_pairs2 = [result_pairs[i] for i in range(0, len(result_pairs)) if i % 2 == 1]
    '''
    print(sets[0])
    print(lens[0:1])
    print(get_size(sets[0]))
    print(sets[10])
    print(lens[20:21])
    print(get_size(sets[10]))
    print(sets[50])
    print(lens[100:101])
    print(get_size(sets[50]))
    '''
    target_ct = 2 ** -75
    best_size = 100000
    best_set = {'a': 'fail'}
    fail = 0
    for i in tqdm(range(len(sets))):
        if fails[2 * i] < target_ct:
            if min(lens[2 * i:2 * i + 2]) < best_size or (min(lens[2 * i:2 * i + 2]) == best_size and fails1[i] > fail):
                best_size = min(lens[2 * i:2 * i + 2])
                best_set = sets[i]
                fail = fails1[i]
    print(best_set)
    print(best_size)
    print(ceil((log(best_set['q'], 2) * best_set['n'] * best_set['m_bar'] + log(best_set['p'], 2) * min(
        [best_set['w'], best_set['n_bar'] * best_set['m_bar']])) / 8))
    print(fail)
    plotting()
    get_smallest(name)


save = True
names = ["results640M=128", "results976M=192", "results1344M=256"]
# main(names[0])
# main(names[1])
# main(names[2])
# main(names[3])


if __name__ == '__main__':
    i = (int(sys.argv[1]))
    names = ["results640M=128", "results976M=192", "results1344M=256", "results640M=256", "results976M=128",
             "results976M=256", "results1344M=128"]
    main(names[i])
