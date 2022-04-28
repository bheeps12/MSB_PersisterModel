import matplotlib.pyplot as plt
import numpy as np
import random

mean_growth_rate = 5
time_steps = 400
basal_death_prob = 0.05 #suggested: 0.05
max_cells = 100000
initial_cells = 300
death_prob= 0.2

#drugs
drug_regime = 1 #-1 for no drugs, 0 for continuous, 1 for intermittent step and 2 for zig zag

if drug_regime==0:
    drug_factor = 0.9
elif drug_regime==1:
    time_period_on = 180
    time_period_off = 30
    max_drug = 0.9
    min_drug = 0.05
elif drug_regime==2:
    time_period = 50
    max_drug = 1
    min_drug = 0.5

array_cells = np.zeros(max_cells)
array_cells[initial_cells:] = -1
growth_rates = np.zeros(max_cells)
growth_rates[0:initial_cells] = np.random.poisson(mean_growth_rate, initial_cells)
np.savetxt('growth.txt', growth_rates)
array_cells[0:initial_cells]=growth_rates[0:initial_cells]

drug_chart = np.zeros(time_steps)
population = np.zeros(time_steps)


for i in range(time_steps):
    print(i)
    if drug_regime==1:
        if (i%(time_period_on+time_period_off) < time_period_on):
            drug_factor= max_drug
        else:
            drug_factor= min_drug
    elif drug_regime==2:
        if (int(i/time_period)%2==0):
            drug_factor=(max_drug-min_drug)/time_period*(i%time_period)+min_drug
        else:
            drug_factor=(max_drug-min_drug)/time_period*(time_period-i%time_period)+min_drug
    live_indices = np.argwhere(array_cells!=-1).T[:][0]
    probabs = np.zeros(max_cells)
    probabs[live_indices] = np.random.uniform(size =np.shape(live_indices)[0])
    prob_death = np.arange(0,max_cells)[np.argwhere((probabs<=death_prob) & (probabs>0))].T
    prob_growth = np.argwhere(probabs>death_prob)
    #growth step
    array_cells[prob_growth] = array_cells[prob_growth]-1
    #death step
    for k in prob_death[:][0]:
        if drug_regime==-1:
            num = random.uniform(0,1)
            if num<basal_death_prob:
                array_cells[k]=-1
                growth_rates[k]=0
        elif drug_regime==0 or drug_regime==1 or drug_regime==2:
            drug_death_prob = drug_factor*1/(1+(growth_rates[k]/(2*mean_growth_rate))**3)*(1-basal_death_prob)+basal_death_prob
            num = random.uniform(0,1)
            if num<drug_death_prob:
                array_cells[k]=-1
                growth_rates[k]=0
    #division step
    if 1 in array_cells:
        div_indices = np.argwhere(array_cells==0).T[:][0]
        for j in div_indices:
            daughter_growth = np.random.poisson(growth_rates[j],2)
            growth_rates[j]=daughter_growth[0]
            array_cells[j]=daughter_growth[0]
            if (np.shape(np.argwhere(array_cells==-1))[0]==0):
                np.savetxt('Henlo.txt', array_cells)
                print('Insufficient max cells')
            new_loc = np.argwhere(array_cells==-1)[0]
            growth_rates[new_loc]=daughter_growth[1]
            array_cells[new_loc]= daughter_growth[1]
    if i%5==0:
        points = np.where(array_cells>-1)
        plt.hist(growth_rates[np.where(growth_rates>0)], bins = 20, edgecolor = 'black')
        plt.title('Time = '+str(i))
        plt.xlabel('Growth Periods')
        plt.ylabel('Occurrence')
        plt.savefig('Hist_'+str(i)+'.png')
        plt.close()
    if drug_regime !=-1:
        drug_chart[i]=drug_factor
    population[i]=live_indices.shape[0]
    print(live_indices.shape[0])
if drug_regime !=-1:
    plt.plot(np.arange(0,time_steps),drug_chart)
    plt.title('Drug levels over time')
    plt.xlabel('Time step')
    plt.ylabel('Drug factor')
    plt.savefig('Drugchart.png')
    plt.close()

plt.plot(np.arange(0,time_steps),population)
plt.title('Population over Time')
plt.xlabel('Time steps')
plt.ylabel('Number of cells in the population')
plt.savefig('population.png')
plt.close()