#!/usr/bin/env python3.9
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import multinomial
from typing import List
import pandas as pd

###Calculate and plot the running probability of state A and state B, assuming that you initially start instate A
t_array = np.linspace(0,6,601)
file=open("probs.txt", "a")
for t in t_array:
	Pa = 1/6 + np.multiply((1-1/6), np.exp(np.multiply(-6, t)))
	Pb = 5/6 + np.multiply((0-5/6), np.exp(np.multiply(-6, t)))
	Pa_array = np.array([Pa])
	Pb_array = np.array([Pb])
	data = np.array([Pa_array, Pb_array])
	data = data.T
	np.savetxt(file, data, fmt=['%f','%f'])
file.close()

data_array = np.loadtxt("probs.txt")
state_prob = pd.DataFrame(data_array)
state_prob.plot()
plt.xlabel("Time (cs)")
plt.ylabel("Probability")
plt.show()


###Plot the time course of the system 
p_initial = np.array([1, 0])
p_transition = np.array(
    [[1/6, 5/6],
     [1/6, 5/6]]
)

def equilibrium_distribution(p_transition):
    n_states = p_transition.shape[0]
    A = np.append(
        arr=p_transition.T - np.eye(n_states),
        values=np.ones(n_states).reshape(1, -1),
        axis=0
    )
    b = np.transpose(np.array([0] * n_states + [1]))
    p_eq = np.linalg.solve(
        a=np.transpose(A).dot(A),
        b=np.transpose(A).dot(b)
    )
    return p_eq

def markov(p_initial: np.array, p_transition: np.array, sequence_length: int) -> List[int]:

    if p_initial is None:
        p_initial = equilibrium_distribution(p_transition)
    initial_state = list(multinomial.rvs(1, p_initial)).index(1)

    states = [initial_state]
    for _ in range(sequence_length - 1):
        p_tr = p_transition[states[-1]]
        new_state = list(multinomial.rvs(1, p_tr)).index(1)
        states.append(new_state)
    return states

states = markov(p_initial, p_transition, sequence_length=600)
fig, ax = plt.subplots(figsize=(12, 4))
plt.plot(states)
plt.xlabel("Time (cs)")
plt.ylabel("State")
plt.yticks([0, 1])
plt.show()

###Plot the correlation function of the signal function
def autocorr(states):
    result = np.correlate(states, states, mode='full')
    return result[result.size // 2:]

acorr = autocorr(states)
plt.plot(acorr / float(acorr.max()))
plt.xlabel("Time (cs)")
plt.ylabel("Autocorrelation")
plt.show()

