from numpy import *
from math import *
from random import randint
import matplotlib.pyplot as plt
import pandas as pd
#alice_base = 0 (0(1) or 180(0)) x-basis
#alice_base = 1 (90(1) or 270(0)) p-basis
#bob_base = 0 measures x (0), 1 measures p (90)
class coh_state:

    def __init__(self,x,p):
        self.x = x 
        self.p = p
        self.alpha = complex(x,p) 
        self.varmat = matrix([[1/4,0],[0,1/4]])
        self.avgphotons = pow(abs(x),2) + pow(abs(p),2)

plt.rcParams.update({'font.size': 12})

def transmit(alice_vals):

    trans_states = []
    for i in range(len(alice_vals)):
        if alice_vals[i] == 0:
            trans_states.append(coh_state(sqrt(n_sig),0))
        elif alice_vals[i] == 90:
            trans_states.append(coh_state(0,sqrt(n_sig)))
        elif alice_vals[i] == 180:
            trans_states.append(coh_state(-sqrt(n_sig),0))
        else:
            trans_states.append(coh_state(0,-sqrt(n_sig)))

    return trans_states

def receive(trans_states, bob_vals):

    measured_vals = []
    for i in range(len(bob_vals)):
        if bob_vals[i] == 0:
            measured_vals.append(random.normal((trans_states[i].x*cos(err_0)+trans_states[i].p*sin(err_0)), sqrt(trans_states[i].varmat[0,0])))
        else:
            measured_vals.append(random.normal((trans_states[i].x*sin(err_90)+trans_states[i].p*cos(err_90)), sqrt(trans_states[i].varmat[1,1])))

    return measured_vals

def alice_generate(N):
    
    rand_array = [];
    for i in range(N):
        rand_array.append(90*randint(0,3))
    
    return rand_array

def bob_generate(N):
    
    rand_array = [];
    for i in range(N):
        rand_array.append(90*randint(0,1))
    
    return rand_array

def sifting(alice_vals, bob_vals, measured_vals):

    a_key = []
    a_temp = []
    b_temp = []
    b_key = []
    pos = []
    count = 0
    for i in range(len(alice_vals)):
        if abs(alice_vals[i] - bob_vals[i]) == 0 or abs(alice_vals[i] - bob_vals[i]) == 180:
                pos.append(i)
                if alice_vals[i] == 0 or alice_vals[i] == 90:
                    a_temp.append(1)
                else:
                    a_temp.append(0)
                b_temp.append(measured_vals[i])

    for j in range(len(b_temp)):
        if b_temp[j] <= x_neg:
            b_key.append(0)
            a_key.append(a_temp[j])
        if b_temp[j] >= x_pos:
            b_key.append(1)
            a_key.append(a_temp[j])
    
    count = len(a_temp) - len(a_key)
    
    return a_key, b_key, pos, count

def noisy_channel(states):

    for i in range(len(states)):
        states[i].x = sqrt(ch_trans*n_det)*states[i].x
        states[i].p = sqrt(ch_trans*n_det)*states[i].p
        states[i].varmat = 0.25*matrix([[(ch_trans*n_det)*(1+4*prep_err) + (1-ch_trans)*ch_nois + (1-n_det)*det_nois,0],[0,(ch_trans*n_det)*(1+4*prep_err) + (1-ch_trans)*ch_nois + (1-n_det)*det_nois]])
    
    return states

def prob_dist(alice_vals,bob_vals,measurements):

    vals0 = []
    vals90 = []
    vals180 = []
    vals270 = []
    for i in range(len(alice_vals)):
        if alice_vals[i]-bob_vals[i] == 0:
            vals0.append(measurements[i])
        elif  alice_vals[i]-bob_vals[i] == 90:
            vals90.append(measurements[i])
        elif abs(alice_vals[i]-bob_vals[i]) == 180:
            vals180.append(measurements[i])
        else:
            vals270.append(measurements[i])
    print(mean(vals0))
    print(std(vals0))
    print(mean(vals90))
    print(std(vals90))
    print(mean(vals180))
    print(std(vals180))
    print(mean(vals270))
    print(std(vals270))
    #fig, axs = plt.subplots(2, 1)
    n0, x0 = histogram(vals0, bins = 150)
    n90, x90 = histogram(vals90, bins = 150)
    n180, x180 = histogram(vals180, bins = 150)
    n270, x270 = histogram(vals270, bins = 150)
    bin_centers0 = 0.5*(x0[1:]+x0[:-1])
    bin_centers90 = 0.5*(x90[1:]+x90[:-1])
    bin_centers180 = 0.5*(x180[1:]+x180[:-1])
    bin_centers270 = 0.5*(x270[1:]+x270[:-1])
    #axs[0].plot(bin_centers,n)
    plt.scatter(bin_centers0,n0,marker='.',label='$0^0$')
    plt.scatter(bin_centers90,n90,marker='x',label='$90^0$')
    plt.scatter(bin_centers180,n180,marker='+',label='$180^0$')
    plt.scatter(bin_centers270,n270,marker='^',label='$270^0$')
    leg = plt.legend()
    plt.xlabel('Quadrature $X_\phi$', fontsize = 15)
    plt.ylabel('Number of counts', fontsize = 15)
    plt.show()
    #plt.savefig('probplot.pdf')
    
def prob_dist_new(alice_vals,bob_vals,measurements):

    vals_corr = []
    vals_wrong = []
    for i in range(len(alice_vals)):
        if abs(alice_vals[i] - bob_vals[i])==0 or abs(alice_vals[i] - bob_vals[i])==180:
            vals_corr.append(measurements[i])
        else:
            vals_wrong.append(measurements[i])

    ncorr, xcorr = histogram(vals_corr, bins = 150)
    nwrong, xwrong = histogram(vals_wrong, bins = 150)
    bin_centers_corr = 0.5*(xcorr[1:]+xcorr[:-1])
    bin_centers_wrong = 0.5*(xwrong[1:]+xwrong[:-1])
    plt.scatter(bin_centers_corr,ncorr,marker='.',label='Correct')
    plt.scatter(bin_centers_wrong,nwrong,marker='.',label='Wrong')
    plt.xlabel('Quadrature $X_\phi$', fontsize = 15)
    plt.ylabel('Number of counts', fontsize = 15)
    plt.show()

def meas_var(bob_vals,measurements):
    meas_0 = [];
    meas_90 = [];
    for i in range(len(bob_vals)):
        if bob_vals[i] == 0:
            meas_0.append(measurements[i])
        else:
            meas_90.append(measurements[i])
    print(var(meas_0))
    print(var(meas_90))

def covariances(alice_vals, bob_vals, measurements):
    a_q = [];
    a_p = [];
    b_q = [];
    b_p = [];
    for i in range(len(alice_vals)):
        if alice_vals[i] == 0:
            a_q.append(sqrt(n_sig))
        elif alice_vals[i] == 180:
            a_q.append(-sqrt(n_sig))
        elif alice_vals[i] == 90:
            a_p.append(sqrt(n_sig))
        else: 
            a_p.append(-sqrt(n_sig))
    for i in range(len(bob_vals)):
        if bob_vals[i] == 0:
            b_q.append(measurements[i])
        else:
            b_p.append(measurements[i])
    print(mean(a_q))
    print(var(a_q))
    print(cov(a_q[0:15000],b_q[0:15000],False))
    print(cov(a_p[0:15000],b_p[0:15000],False))
    print(cov(a_q[0:15000],b_p[0:15000],False))
    print(pow(cov(a_q[0:15000],b_q[0:15000],False)[0,1]/(var(a_q)),2))
    # print(len(a_q))
    # print(len(a_p))
    # print(len(b_q))
    # print(len(b_p))


def DMCVQKD(N):

    alice_vals = alice_generate(N)
    bob_vals = bob_generate(N)
    states = transmit(alice_vals)
    noisy_states = noisy_channel(states)
    measurements = receive(noisy_states, bob_vals)
    alice_key , bob_key, positions, deleted_count = sifting(alice_vals, bob_vals, measurements)
    int_err = 0
    for i in range(len(alice_key)):
        if alice_key[i] != bob_key[i]:
            int_err = int_err+1
    df = pd.DataFrame({"alice phase" : alice_vals, "bob phase" : bob_vals, "measurements" : measurements})
    df.to_csv("output.csv", index=True)
    df2 = pd.DataFrame({ "alice key" : alice_key, "bob key" : bob_key})
    df2.to_csv("key.csv", index=True)
    print(len(alice_key))
    print(len(bob_key))
    print("deleted measurements=",deleted_count)
    print("error=",int_err)
    prob_dist(alice_vals,bob_vals,measurements)
    prob_dist_new(alice_vals,bob_vals,measurements)
    meas_var(bob_vals,measurements)
    print('break')
    #covariances(alice_vals, bob_vals, measurements)
    # print('done')
    # print(1-(err_count/N))
    # print(int_err/N)
    #print(alice_vals)
    #print(bob_vals)
    #print(measurements)
    #print(positions)
    #print(alice_key)
    #print(bob_key)
    #return 1-(err_count/N), int_err/N
            

#a_vals=[0,90,180,270]
#a_b = ['h','d']
#a_bits=[0,1,0,1]
#b_vals = [0,90,0,90]
#N=4
#a_bases = generate(N)
#a_bits = generate(N)
#b_bases = generate(N)
#states = transmit(a_vals)
#m_vals = receive(states,b_vals)
#print(m_vals)

print("start")
n_sig = 1
x_pos = 0.1
x_neg = -0.1
ch_trans = 1
ch_nois = 0
n_det = 1
det_nois = 0
err_0 = 0
err_90 = 0
prep_err = 0
N = 100000
DMCVQKD(N)

#error part
# N_sig = [0.1, 0.2, 0.5, 1, 2, 10, 100, 1000]
# fig, ax = plt.subplots(2,1)
# for n in N_sig:
#     n_sig = n
#     x_thresh = arange(0,10,0.1)
#     x_posarr = x_thresh
#     x_negarr = -x_thresh
#     post_eff = []
#     post_err = []
#     for j in range(len(x_thresh)):
#         x_pos = x_posarr[j]
#         x_neg = x_negarr[j]
#         eff, err = DMCVQKD(N)
#         post_eff.append(eff)
#         post_err.append(err)
    
    
#     ax[0].scatter(x_thresh,post_eff, label=fr"$N_{{sig}}$ = {n}")
    
#     ax[1].scatter(x_thresh,post_err,label=fr"$N_{{sig}}$ = {n}")

# ax[0].set_xlabel('Threshold $X_{o}$', fontsize=15)
# ax[0].set_ylabel('Post selection efficiency $p_{d}$', fontsize=15)
# ax[0].set_yscale('log')

# plt.yscale('log')
# plt.xlabel('Threshold $X_{o}$', fontsize=15)
# plt.ylabel('Bit Error Rate', fontsize=15)
# ax[1].legend(loc = 'upper right', ncol=2)
# #fig.legend(ncol = 2)
# plt.show()

# for i in range(5):
#     p_d = []
#     for j in range(len(x_thresh)):
#         x_pos = x_posarr[j]
#         x_neg = x_negarr[j]
#         err = DMCVQKD(N)
#         p_d.append(err)
#     post_err.append(p_d)
# print(post_err)
# post_mean = []
# post_var = []
# for i in range(len(x_thresh)):
#     arr_tmp = []
#     for j in range(5):
#         arr_tmp.append(post_err[j][i])
#     post_mean.append(mean(arr_tmp))
#     post_var.append(var(arr_tmp))

# print(x_thresh)
# print(post_mean)
# print(post_var)
# plt.errorbar(x_thresh,post_mean,yerr=post_var,fmt='')
# plt.subplot(211)
# plt.scatter(x_thresh,post_eff)
# plt.yscale('log')
# plt.xlabel('Threshold $X_{o}$')
# plt.ylabel('Post selection efficiency $p_{d}$')
# plt.subplot(212)
# plt.scatter(x_thresh,post_err)
# plt.yscale('log')
# plt.xlabel('Threshold $X_{o}$')
# plt.ylabel('Bit Error Rate')
# plt.show()



