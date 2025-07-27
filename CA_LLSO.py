#**************************************************************************************************
# Implementation of Classifier-Assisted Level-Assisted Learning Swarm Optimizer (CA-LLSO)
#
# Feng-Feng Wei, Wei-Neng Chen, Qiang Yang, Jeremiah Deng, Xiao-Nan Luo, Hu Jin and Jun Zhang, 
# "A Classifier-Assisted Level-Based Learning Swarm Optimizer for Expensive Optimization", 
# IEEE Transactions on Evolutionary Computation, accepted, 2020.
#
# Modified by Congji Han for parameter calibration 
#
#**************************************************************************************************
#
# REQUIREMENTS:
# The following libraries should be intalled before running the code:
#     python v2.7.15 or newer versions
#     numpy v1.15.4 or new versions
#     scipy v1.1.0 or newer versions
#     sklearn v0.20.3
#     pyDOE v0.3.8
#
#**************************************************************************************************

import numpy as np
import time
import math
import random
import os
import uuid
import multiprocessing
import subprocess
from sklearn.ensemble import GradientBoostingClassifier
from pyDOE import *
from numpy import linalg as LA
# import warnings
# warnings.filterwarnings("ignore", category=DeprecationWarning)

np.set_printoptions(threshold=np.inf)
D = 4 # dimensionality
bound1 = 0.02 # m_road(0-0.02)
bound2 = 0.2 # m_paddy(0-0.2)
bound3 = 0.3 # m_building(0-0.3)
bound4 = 2 # cd_building(0-2)
bound = np.array([bound1, bound2, bound3, bound4])

if D>=50:
    NP = 200 # the population size
    count_l1 = 150
elif D<50:
    # NP = 100
    # count_l1 = 75
    # minimal sampling test
    NP = 10
    count_l1 = 7
max_gen = 25 # the maxinal evolution generation
max_times = 20 # the total independent run times
NL = 4 # the number of layers
phi = 0.4
pop_x = None
pop_v = np.zeros((NP, D))
pop_y = None
db_x = None
db_v = np.zeros((NP,D))
db_y = None
rank = np.zeros((NP, ))
ranked_x = np.zeros((NP, D))
ranked_v = np.zeros((NP, D))
ranked_y = np.zeros((NP, ))

offspring_rank = np.zeros((NP, ))
generated_level_1 = np.zeros((NP, D))
generated_level_1_v = np.zeros((NP, D))
eval_per_gen = int(100/max_gen)

def init_lock(l):
    global log_lock
    log_lock = l
    
def evaluate_single(row):
    m_r, m_p, m_b, cd_b = row
    unique_id = uuid.uuid4().hex
    input_file = f"ufm/readfile/data/2D/input_{unique_id}.dat"
    output_file = f"ufm/readfile/data/RMSE/output_{unique_id}.dat"
    with open(input_file, "w") as f:
        f.write(f"{m_r} {m_p} {m_b} {cd_b}\n")
    subprocess.run(["./ufm/MP", input_file, output_file])
    with open(output_file, "r") as f:
        line = f.readline().strip()
        numbers = line.split()
        rmse = float(numbers[-1])
    
    log_message = f"Input: {m_r} {m_p} {m_b} {cd_b} -> Output: {line}\n"
    global log_lock
    with log_lock:
        with open("log.dat", "a") as log:
            log.write(log_message)     
    os.remove(input_file)
    os.remove(output_file)
    return rmse

def evaluate_rmse_parallel(db_x):
    if not os.path.exists("log.dat"):
        with open("log.dat", "w") as log:
            log.write("Input (m_r, m_p, m_b, cd_b) -> Output (RMSE)\n")
    lock = multiprocessing.Lock()
    with multiprocessing.Pool(initializer=init_lock, initargs=(lock,)) as pool:
        fitness_values = pool.map(evaluate_single, db_x)
    return np.array(fitness_values)

def update(ori_x, ori_v, x_1, x_2):
    r1 = np.random.random((1,D))
    r2 = np.random.random((1,D))
    r3 = np.random.random((1,D))
    new_v = r1 * ori_v + r2 * (x_1 - ori_x) + phi * r3 * (x_2 - ori_x)
    new_x = ori_x + new_v
    for i in range(D):
        if new_x[0,i] < 0:
            new_x[0,i] = 0
        if new_x[0,i] > bound[i]:
            new_x[0,i] = bound[i]
    return new_x, new_v

if __name__ == "__main__":
    start_time = time.time()
    best_ever = np.zeros((max_times, max_gen))
    con_count = np.zeros((max_times, max_gen))
    for times in range(max_times):
        #initialization
        db_x = lhs(D, samples=NP, criterion = 'center') # LHS
        db_x[:, 0] = db_x[:, 0] * 0.02 # b1
        db_x[:, 1] = db_x[:, 1] * 0.2 # b2
        db_x[:, 2] = db_x[:, 2] * 0.3 # b3
        db_x[:, 3] = db_x[:, 3] * 2 # b4
        db_v = np.zeros((NP,D))
        db_y = evaluate_rmse_parallel(db_x)
        #iteration
        for itera in range(max_gen): 
            ORIGIN = np.concatenate((db_x,db_v,db_y.reshape(db_y.shape[0],1)),axis=1) 
            ORIGIN = ORIGIN[ORIGIN[:,-1].argsort()] 
            pop_x = ORIGIN[0:NP,0:D]
            pop_v = ORIGIN[0:NP,D:D*2]
            pop_y = ORIGIN[0:NP,-1]
            level1_count = 0
            #sort_y = np.argsort(pop_y)
            #print(sort_y)
            for i in range(NL-1):
                ranked_x[int(NP/NL)*i:int(NP/NL)*(i+1),:] = pop_x[int(NP/NL)*i:int(NP/NL)*(i+1),:]
                ranked_v[int(NP/NL)*i:int(NP/NL)*(i+1),:] = pop_v[int(NP/NL)*i:int(NP/NL)*(i+1),:]
                ranked_y[int(NP/NL)*i:int(NP/NL)*(i+1)] = pop_y[int(NP/NL)*i:int(NP/NL)*(i+1), ]
                rank[int(NP/NL)*i:int(NP/NL)*(i+1)] = i
            ranked_x[int(NP/NL)*(NL-1):,:] = pop_x[int(NP/NL)*(NL-1):,:]
            ranked_v[int(NP/NL)*(NL-1):,:] = pop_v[int(NP/NL)*(NL-1):,:]
            ranked_y[int(NP/NL)*(NL-1):,] = pop_y[int(NP/NL)*(NL-1):,]
            rank[int(NP/NL)*(NL-1):] = NL - 1
            clf = GradientBoostingClassifier().fit(ranked_x, rank) # GBC
            #update
            offspring_x = np.zeros((NP, D))
            offspring_v = np.zeros((NP, D))
            offspring_rank = np.zeros((NP, ))
            for i in range(NP-(NL-1)*int(NP/NL)):
                rl1 = np.random.randint(0, NL-1)
                rl2 = np.random.randint(0, NL-1)
                while rl1 == rl2:
                    rl2 = np.random.randint(0, NL-1)
                if rl1 > rl2:
                    tmp = rl1
                    rl1 = rl2
                    rl2 = tmp
                k1 = np.random.randint(0, int(NP/NL))
                k2 = np.random.randint(0, int(NP/NL))
                lev_1 = rl1 * int(NP/NL) + k1
                lev_2 = rl2 * int(NP/NL) + k2
                offspring_x[-(i+1),:], offspring_v[-(i+1),:] = update(ranked_x[-(i+1),:], ranked_v[-(i+1),:], ranked_x[lev_1,:], ranked_x[lev_2,:])
                offspring_rank[-(i+1),] = clf.predict(offspring_x[-(i+1),:].reshape(1, D))
                if level1_count < count_l1 and offspring_rank[-(i+1),] == 0:
                    generated_level_1[level1_count,:] = offspring_x[-(i+1),:]
                    generated_level_1_v[level1_count,:] = offspring_v[-(i+1),:]
                    level1_count = level1_count + 1

            #update the third to the last-1 level
            for i in range(2, NL - 1):
                for j in range(int(NP/NL)):
                    rl1 = np.random.randint(0, i)
                    rl2 = np.random.randint(0, i)
                    while rl1 == rl2:
                        rl2 = np.random.randint(0, i)
                    if rl1 > rl2:
                        tmp = rl1
                        rl1 = rl2
                        rl2 = tmp
                    k1 = np.random.randint(0, int(NP/NL))
                    k2 = np.random.randint(0, int(NP/NL))
                    cur = i * int(NP/NL) + j
                    lev_1 = rl1 * int(NP/NL) + k1
                    lev_2 = rl2 * int(NP/NL) + k2

                    offspring_x[cur,:], offspring_v[cur,:] = update(ranked_x[cur,:], ranked_v[cur,:], ranked_x[lev_1,:], ranked_x[lev_2,:])
                    offspring_rank[cur,] = clf.predict(offspring_x[cur,:].reshape(1, D))
                    if level1_count < count_l1 and offspring_rank[cur,] == 0:
                        generated_level_1[level1_count,:] = offspring_x[cur,:]
                        generated_level_1_v[level1_count,:] = offspring_v[cur,:]
                        level1_count = level1_count + 1

            #update the second level
            for i in range(int(NP/NL)):
                r1 = np.random.randint(0, int(NP/NL))
                r2 = np.random.randint(0, int(NP/NL))
                while r1 == r2:
                    r2 = np.random.randint(0, int(NP/NL))
                if r1 > r2:
                    tmp = r1
                    r1 = r2
                    r2 = tmp
                j = i + int(NP/NL)
                offspring_x[j,:], offspring_v[j,:] = update(ranked_x[j,:], ranked_v[j,:], ranked_x[r1,:], ranked_x[r2,:])
                offspring_rank[j,] = clf.predict(offspring_x[j,:].reshape(1, D))
                if level1_count < count_l1 and offspring_rank[j,] == 0:
                    generated_level_1[level1_count,:] = offspring_x[j,:]
                    generated_level_1_v[level1_count,:] = offspring_v[j,:]
                    level1_count = level1_count + 1
            #print(offspring_x.shape,offspring_v.shape,offspring_rank.shape)
            #copy the first level
            offspring_x[0:int(NP/NL),:] = ranked_x[0:int(NP/NL),:]
            offspring_v[0:int(NP/NL),:] = ranked_v[0:int(NP/NL),:]
            offspring_rank[0:int(NP/NL),] = 0
            # L1-exploitation
            while level1_count < count_l1:
                con_count[times,itera] = con_count[times,itera]+1
                #print(offspring_x.shape,offspring_v.shape,offspring_rank.shape)
                OFF = np.concatenate((offspring_x,offspring_v,offspring_rank.reshape(offspring_rank.shape[0],1)),axis=1)
                OFF = OFF[OFF[:,-1].argsort()]
                offspring_x = OFF[:,0:D]
                offspring_v = OFF[:,D:D*2]
                offspring_rank = OFF[:,-1]
                nl0 = np.sum(offspring_rank==0)
                nl1 = np.sum(offspring_rank==1)
                nl2 = np.sum(offspring_rank==2)
                nl3 = np.sum(offspring_rank==3)
                nl = [nl0, nl0+nl1, nl0+nl1+nl2, nl0+nl1+nl2+nl3]
                #print(level1_count, nl)
                #update the last level
                if nl3 > 0 and level1_count < count_l1:
                    for i in range(nl3):
                        k1 = np.random.randint(0, nl0+nl1+nl2)
                        k2 = np.random.randint(0, nl0+nl1+nl2)
                        while k1 == k2:
                            k2 = np.random.randint(0, nl0+nl1+nl2)
                        if (offspring_rank[k1] == offspring_rank[k2]) and nl[int(offspring_rank[k1])] == (NP-nl3):
                            lev_1 = k1
                            lev_2 = k2
                        else:
                            while offspring_rank[k1] == offspring_rank[k2]:
                                k2 = np.random.randint(0, nl0+nl1+nl2)
                            if offspring_rank[k1] > offspring_rank[k2]:
                                lev_1 = k2
                                lev_2 = k1
                            else:
                                lev_1 = k1
                                lev_2 = k2
                        offspring_x[-(i+1),:], offspring_v[-(i+1),:] = update(offspring_x[-(i+1),:], offspring_v[-(i+1),:], offspring_x[lev_1,:], offspring_x[lev_2,:])
                        offspring_rank[-(i+1),] = clf.predict(offspring_x[-(i+1),:].reshape(1, D))
                        if level1_count < count_l1 and offspring_rank[-(i+1),] == 0:
                            generated_level_1[level1_count,:] = offspring_x[-(i+1),:]
                            generated_level_1_v[level1_count,:] = offspring_v[-(i+1),:]
                            level1_count = level1_count + 1
                #update the third to the last-1 level
                for i in range(2, 3):
                    if nl2 > 0 and level1_count < count_l1:
                        for j in range(nl2):
                            cur = nl0+nl1
                            cur = cur + j
                            if nl1 == 0:
                                k1 = np.random.randint(0, nl0)
                                k2 = np.random.randint(0, nl0)
                                while r1 == r2:
                                    r2 = np.random.randint(0, nl0)
                                lev_1 = k1
                                lev_2 = k2
                            else:
                                k1 = np.random.randint(0, nl0)
                                k2 = np.random.randint(0, nl1)
                                lev_1 = k1
                                lev_2 = k2+nl0

                            offspring_x[cur,:], offspring_v[cur,:] = update(offspring_x[cur,:], offspring_v[cur,:], offspring_x[lev_1,:], offspring_x[lev_2,:])
                            offspring_rank[cur,] = clf.predict(offspring_x[cur,:].reshape(1, D))
                            if level1_count < count_l1 and offspring_rank[cur,] == 0:
                                generated_level_1[level1_count,:] = offspring_x[cur,:]
                                generated_level_1_v[level1_count,:] = offspring_v[cur,:]
                                level1_count = level1_count + 1
                            
                #update the second level
                if nl1 > 0 and level1_count < count_l1:
                    for i in range(nl1):
                        r1 = np.random.randint(0, nl0)
                        r2 = np.random.randint(0, nl0)
                        while r1 == r2:
                            r2 = np.random.randint(0, nl0)
                        j = i + nl0
                        offspring_x[j,:], offspring_v[j,:] = update(offspring_x[j,:], offspring_v[j,:], offspring_x[r1,:], offspring_x[r2,:])
                        offspring_rank[j,] = clf.predict(offspring_x[j,:].reshape(1, D))
                        if level1_count < count_l1 and offspring_rank[j,] == 0:
                            generated_level_1[level1_count,:] = offspring_x[j,:]
                            generated_level_1_v[level1_count,:] = offspring_v[j,:]
                            level1_count = level1_count + 1
            dist = np.zeros((count_l1, int(NP/NL)))
            for i in range(count_l1):
                for j in range(int(NP/NL)):
                    dist[i,j] = np.sqrt(np.sum(np.power(generated_level_1[i,:]-ranked_x[j,:], 2))) #

            max_dist = np.amax(dist, axis=1)
            #max_dist = np.sum(dist, axis=1) / int(NP/NL)
            mini_max_dist_ind = np.argsort(max_dist)
            #print(mini_max_dist_ind)
            x_tail = np.zeros((eval_per_gen, D))
            v_tail = np.zeros((eval_per_gen, D))
            x_tail[0:eval_per_gen-1,:] = generated_level_1[mini_max_dist_ind[0:eval_per_gen-1],:]
            v_tail[0:eval_per_gen-1,:] = generated_level_1_v[mini_max_dist_ind[0:eval_per_gen-1],:]
            rand_num1 = np.random.randint(eval_per_gen-1, count_l1)
            x_tail[eval_per_gen-1,:] = generated_level_1[mini_max_dist_ind[rand_num1],:]
            v_tail[eval_per_gen-1,:] = generated_level_1_v[mini_max_dist_ind[rand_num1],:]
            #print(selected_function(generated_level_1[mini_max_dist_ind,:]))
            #dynamic increase the database
            new_indices = []
            for ii in range(eval_per_gen):
                 found = False
                 for jj in range(db_x.shape[0]):
                      if np.all(x_tail[ii] == db_x[jj]):
                           found = True
                           break
                 if not found:
                      new_indices.append(ii)
                      
            if new_indices:
                 new_x = x_tail[new_indices, :] 
                 new_v = v_tail[new_indices, :] 
                 new_y = evaluate_rmse_parallel(new_x)
                 db_x = np.concatenate((db_x, new_x), axis=0)
                 db_v = np.concatenate((db_v, new_v), axis=0)
                 db_y = np.concatenate((db_y, new_y), axis=0)

            # count = 0
            # for ii in range(eval_per_gen):
            #     flag = 0
            #     for jj in range(db_x.shape[0]):
            #         if all(x_tail[ii] == db_x[jj]):
            #             flag = 1
            #             count = count + 1
            #             break
            #     if flag == 0:
            #         db_x = np.concatenate((db_x, x_tail[ii,:].reshape(1,D)), axis=0)
            #         db_v = np.concatenate((db_v, v_tail[ii,:].reshape(1,D)), axis=0)
            #         db_y = np.concatenate((db_y, evaluate_nse(x_tail[ii,:].reshape(1,D))), axis=0)
            #print(db_x.shape, db_v.shape, db_y.shape)
            best_ever[times,itera] = np.amin(db_y)
            
            best_index = np.argmin(db_y)
            best_dbx = db_x[best_index, :]
            data_to_save = np.hstack((db_x, db_y.reshape(-1, 1))) 
            filename = f"results_gen_{times}_{itera}.dat"
            np.savetxt(filename, data_to_save, fmt="%.6e", header="db_x (D columns) and fitness (last column)")
            # print(itera, min(db_y), count)
            print(itera, min(db_y))
            #print(con_count[times, itera])
        print(times, best_ever[times,-1])

    end_time = time.time()
    print('time cost: ', end_time - start_time, 's')
    best_ever = best_ever.transpose()
    best_average = np.mean(best_ever, axis=1)
    best_std = np.std(best_ever, axis=1)
    print(best_average[-1], best_std[-1])
    best_2_save = np.concatenate((best_ever, best_average.reshape(max_gen,1)),axis=1)
    #np.savetxt('new_L1explotation_f%d_D%d_eval%d_1L%d_ite%d.txt' %(fi+1, D, eval_per_gen, count_l1, max_times), best_2_save)
    #np.savetxt('convergence_speed_f%d_D%d_eval%d_1L%d_ite%d.csv' %(fi+1, D, eval_per_gen, count_l1, max_times), np.mean(con_count, axis=0).transpose())
