# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 10:39:19 2019

@author: Tom
"""

import numpy as np
import matplotlib.pyplot as plt

background = []
rate = 0.0

# VOLUME
FreeProtons = 0.668559
tankRadius = 9500 - 6.35
halfHeight = 10000 - 6.35
vol = np.pi*pow(tankRadius/1000.,2)*(2.*halfHeight/1000.)
print ("Volume of Water",vol,"m^3")

# MASS
# Input masses of components from masses.py 
nKiloTons = np.pi*pow(tankRadius/1000.,2)*(2.*halfHeight/1000.)/1000.
mass_pmt = 4900
mass_cable = 0.0
mass_conc = 2.96e6
mass_gd = 5.66e6 * 0.1
mass_rock = 7.617e6
mass_shield = 0.0
mass_tank = 8.94e4
mass_veto = 701

# ACTIVITIES
# Define Activities of each decay chain in each component from activities.py
# 238U CHAIN
activity238_pmt = 0.538
activity238_veto = 4.24
activity238_gd = 0.113

# 232TH CHAIN
activity232_pmt = 0.541
activity232_veto = 5.42
activity232_gd = 2.26

# 235U CHAIN
activity235 = 2.82
activity235_gd = 2.83

# 40K CHAIN
activity40_pmt = 0.520
activity40_veto =0.520

# 60Co CHAIN
activity60 = 19e-3

# 137Cs CHAIN
activity137 = 0.77e-3

# 222Rn CHAIN
activity222 = 2e-3

# FN CHAIN
activityFN = 0.00809

# IBD RATE
activityIBD = 9.60e-5

# CNO RATE
activityCNO = 0.00011

# EFFICIENCIES
# Input Efficiencies taken from ROOT
eff_concrete_208Tl = 0.0054
eff_concrete_210Bi = 2e-6
eff_concrete_210Tl = 0.0035 * 0.002
eff_concrete_212Bi = 0.00037
eff_concrete_212Pb = 0.00043
eff_concrete_214Bi = 0.0025
eff_concrete_214Pb = 0.00048
eff_concrete_228Ac = 0.0019
eff_concrete_234Pa = 1e-5
eff_concrete_40K = 0.00039
eff_concrete_FTFP = 0.011
eff_concrete_QBBC_EMZ = 0.011
eff_concrete_QBBC_FN = 0.012
eff_concrete_QGSP_BERT_EMV = 0.012
eff_concrete_QGSP_BERT_EMX = 0.011
eff_concrete_QGSP_BERT_FN = 0.012
eff_concrete_QGSP_BIC_FN = 0.012
eff_concrete_QGSP_FTFP = 0.011

eff_gd_207Tl = 0.000
eff_gd_208Tl = 0.051
eff_gd_210Bi = 0.027
eff_gd_210Tl = 0.049 * 0.002
eff_gd_211Bi = 0.000
eff_gd_211Pb = 0.000
eff_gd_212Bi = 0.034
eff_gd_212Pb = 0.026
eff_gd_214Bi = 0.043
eff_gd_214Pb = 0.026
eff_gd_223Fr = 0.000
eff_gd_228Ac = 0.041
eff_gd_231Th = 0.000
eff_gd_234Pa = 0.031

eff_pmt_208Tl = 0.018
eff_pmt_210Bi = 0.018
eff_pmt_210Tl = 0.019 *0.002
eff_pmt_212Bi = 0.030
eff_pmt_212Pb = 0.030
eff_pmt_214Bi = 0.023
eff_pmt_214Pb = 0.032
eff_pmt_228Ac = 0.030
eff_pmt_234Pa = 0.029
eff_pmt_40K = 0.028

eff_rock_208Tl = 0.0005
eff_rock_210Bi = 0.000
eff_rock_210Tl = 0.00024 * 0.002
eff_rock_212Bi = 3.33e-5
eff_rock_212Pb = 3.33e-5
eff_rock_214Bi = 0.0002
eff_rock_214Pb = 3.33e-5
eff_rock_228Ac = 0.0001
eff_rock_234Pa = 0.000
eff_rock_40K = 2e-5
eff_rock_FTFP = 0.012
eff_rock_EMZ_FN = 0.011
eff_rock_FN = 0.012
eff_rock_EMV = 0.011
eff_rock_EMX = 0.012
eff_rock_BERT_FN = 0.012
eff_rock_BIC = 0.012
eff_rock_FTFP_BERT_FN = 0.012

eff_veto_208Tl = 0.084
eff_veto_210Bi = 0.016
eff_veto_210Tl = 0.079 * 0.002
eff_veto_212Bi = 0.089
eff_veto_212Pb = 0.038
eff_veto_214Bi = 0.066
eff_veto_214Pb = 0.046
eff_veto_228Ac = 0.061
eff_veto_234Pa = 0.030
eff_veto_40K = 0.029

eff_water_11003 = 0.0
eff_water_210Bi = 0.025
eff_water_210Tl = 0.056
eff_water_214Bi = 0.043
eff_water_214Pb = 0.025
eff_water_40K = 0.029
eff_water_9003 = 0.0
eff_water_delayedNeutron = 0.098
eff_water_delayedPair = 0.14
eff_water_positron = 0.047

eff_CNO = 0.05

# Reshape the efficiencies into arrays based on decay chain
concrete_40K = [
        "eff_concrete_40K", eff_concrete_40K]
concrete_40K = np.array(concrete_40K)
concrete_40K = np.reshape(concrete_40K,(1,2))

concrete_U238 = [
        "eff_concrete_210Bi", eff_concrete_210Bi,
        "eff_concrete_210Tl", eff_concrete_210Tl,
        "eff_concrete_214Bi", eff_concrete_214Bi,
        "eff_concrete_214Pb", eff_concrete_214Pb,
        "eff_concrete_234Pa", eff_concrete_234Pa]
concrete_U238 = np.array(concrete_U238)
concrete_U238 = np.reshape(concrete_U238,(5,2))

concrete_TH232 = [
        "eff_concrete_228Ac", eff_concrete_228Ac,
        "eff_concrete_212Pb", eff_concrete_212Pb,
        "eff_concrete_212Bi", eff_concrete_212Bi,
        "eff_concrete_208Tl", eff_concrete_208Tl]
concrete_TH232 = np.array(concrete_TH232)
concrete_TH232 = np.reshape(concrete_TH232,(4,2))

concrete_fn = [
        "eff_concrete_FTFP", eff_concrete_FTFP,
        "eff_concrete_QBBC_EMZ", eff_concrete_QBBC_EMZ,
        "eff_concrete_QBBC_FN", eff_concrete_QBBC_FN,
        "eff_concrete_QGSP_BERT_EMV", eff_concrete_QGSP_BERT_EMV,
        "eff_concrete_QGSP_BERT_EMX", eff_concrete_QGSP_BERT_EMX,
        "eff_concrete_QGSP_BERT_FN", eff_concrete_QGSP_BERT_FN,
        "eff_concrete_QGSP_BIC_FN", eff_concrete_QGSP_BIC_FN,
        "eff_concrete_QGSP_FTFP", eff_concrete_QGSP_FTFP]
concrete_fn = np.array(concrete_fn)
concrete_fn = np.reshape(concrete_fn,(8,2))

GD_U238 = [
        "eff_gd_210Bi", eff_gd_210Bi,
        "eff_gd_210Tl", eff_gd_210Tl,
        "eff_gd_214Bi", eff_gd_214Bi,
        "eff_gd_214Pb", eff_gd_214Pb,
        "eff_gd_234Pa", eff_gd_234Pa]
GD_U238 = np.array(GD_U238)
GD_U238 = np.reshape(GD_U238,(5,2))

GD_TH232 = [
        "eff_gd_228Ac", eff_gd_228Ac,
        "eff_gd_212Pb", eff_gd_212Pb,
        "eff_gd_212Bi", eff_gd_212Bi,
        "eff_gd_208Tl", eff_gd_208Tl]
GD_TH232 = np.array(GD_TH232)
GD_TH232 = np.reshape(GD_TH232,(4,2))

GD_U235 = [
        "eff_gd_231Th", eff_gd_231Th,
        "eff_gd_223Fr", eff_gd_223Fr,
        "eff_gd_211Pb", eff_gd_211Pb,
        "eff_gd_211Bi", eff_gd_211Bi,
        "eff_gd_207Tl", eff_gd_207Tl]
GD_U235 = np.array(GD_U235)
GD_U235 = np.reshape(GD_U235,(5,2))

PMT_U238 = [
        "eff_pmt_210Bi", eff_pmt_210Bi,
        "eff_pmt_210Tl", eff_pmt_210Tl,
        "eff_pmt_214Bi", eff_pmt_214Bi,
        "eff_pmt_214Pb", eff_pmt_214Pb,
        "eff_pmt_234Pa", eff_pmt_234Pa]
PMT_U238 = np.array(PMT_U238)
PMT_U238 = np.reshape(PMT_U238,(5,2))

PMT_TH232 = [
        "eff_pmt_228Ac", eff_pmt_228Ac,
        "eff_pmt_212Pb", eff_pmt_212Pb,
        "eff_pmt_212Bi", eff_pmt_212Bi,
        "eff_pmt_208Tl", eff_pmt_208Tl]
PMT_TH232 = np.array(PMT_TH232)
PMT_TH232 = np.reshape(PMT_TH232,(4,2))

PMT_40K = [
        "eff_pmt_40K", eff_pmt_40K]
PMT_40K = np.array(PMT_40K)
PMT_40K = np.reshape(PMT_40K,(1,2))

rock_40K = [
        "eff_rock_40K", eff_rock_40K]
rock_40K = np.array(rock_40K)
rock_40K = np.reshape(rock_40K,(1,2))

rock_U238 = [
        "eff_rock_210Bi", eff_rock_210Bi,
        "eff_rock_210Tl", eff_rock_210Tl,
        "eff_rock_214Bi", eff_rock_214Bi,
        "eff_rock_214Pb", eff_rock_214Pb,
        "eff_rock_234Pa", eff_rock_234Pa]
rock_U238 = np.array(rock_U238)
rock_U238 = np.reshape(rock_U238,(5,2))

rock_TH232 = [
        "eff_rock_228Ac", eff_rock_228Ac,
        "eff_rock_212Pb", eff_rock_212Pb,
        "eff_rock_212Bi", eff_rock_212Bi,
        "eff_rock_208Tl", eff_rock_208Tl]
rock_TH232 = np.array(rock_TH232)
rock_TH232 = np.reshape(PMT_TH232,(4,2))

rock_fn = [
        "eff_rock_FTFP", eff_rock_FTFP,
        "eff_rock_EMZ", eff_rock_EMZ_FN,
        "eff_rock_FN", eff_rock_FN,
        "eff_rock_EMV", eff_rock_EMV,
        "eff_rock_EMX", eff_rock_EMX,
        "eff_rock_FN", eff_rock_FN,
        "eff_rock_BIC", eff_rock_BIC,
        "eff_rock_FTFP_BERT_FN", eff_rock_FTFP_BERT_FN]
rock_fn = np.array(rock_fn)
rock_fn = np.reshape(rock_fn,(8,2))

veto_40K = [
        "eff_veto_40K", eff_veto_40K]
veto_40K = np.array(veto_40K)
veto_40K = np.reshape(veto_40K,(1,2))

veto_U238 = [
        "eff_veto_210Bi", eff_veto_210Bi,
        "eff_veto_210Tl", eff_veto_210Tl,
        "eff_veto_214Bi", eff_veto_214Bi,
        "eff_veto_214Pb", eff_veto_214Pb,
        "eff_veto_234Pa", eff_veto_234Pa]
veto_U238 = np.array(veto_U238)
veto_U238 = np.reshape(veto_U238,(5,2))

veto_TH232 = [
        "eff_veto_228Ac", eff_veto_228Ac,
        "eff_veto_212Pb", eff_veto_212Pb,
        "eff_veto_212Bi", eff_veto_212Bi,
        "eff_veto_208Tl", eff_veto_208Tl]
veto_TH232 = np.array(veto_TH232)
veto_TH232 = np.reshape(PMT_TH232,(4,2))

water_40K = [
        "eff_water_40K", eff_water_40K]
water_40K = np.array(water_40K)
water_40K = np.reshape(water_40K,(1,2))

water_222Rn = [
        "eff_water_210Bi", eff_water_210Bi,
        "eff_water_210Tl", eff_water_210Tl,
        "eff_water_214Bi", eff_water_214Bi,
        "eff_water_214Pb", eff_water_214Pb]
water_222Rn = np.array(water_222Rn)
water_222Rn = np.reshape(water_222Rn,(4,2))

water_ibd = [
        "eff_water_delayedNeutron", eff_water_delayedNeutron,
        "eff_water_delayedPair", eff_water_delayedPair,
        "eff_water_positron", eff_water_positron]
water_ibd = np.array(water_ibd)
water_ibd = np.reshape(water_ibd,(3,2))

# Multiply efficiencies by activity and mass to get rates
# CONCRETE RATE
rate = rate + (float(concrete_40K[0,1])*activity40_pmt*mass_conc)
for i in range(5):
    rate = rate + (float(concrete_U238[i,1])*activity238_pmt*mass_conc)
for i in range(4):
    rate = rate + (float(concrete_TH232[i,1])*activity232_pmt*mass_conc)
for i in range(8):
    rate = rate + (float(concrete_fn[i,1])*activityFN*mass_conc)

# GD RATE
for i in range(5):
    rate = rate + (float(GD_U238[i,1])*activity238_gd)
for i in range(4):
    rate = rate + (float(GD_TH232[i,1])*activity232_gd)
for i in range(5):
    rate = rate + (float(GD_U235[i,1])*activity235_gd)

# PMT RATE
for i in range(5):
    rate = rate + (float(PMT_U238[i,1])*activity238_pmt*mass_pmt)
for i in range(4):
    rate = rate + (float(PMT_TH232[i,1])*activity232_pmt*mass_pmt)
rate = rate + (float(PMT_40K[0,1])*activity40_pmt*mass_pmt)

# ROCK RATE
rate = rate + (float(rock_40K[0,1])*activity40_pmt*mass_rock)
for i in range(5):
    rate = rate + (float(rock_U238[i,1])*activity238_pmt*mass_rock)
for i in range(4):
    rate = rate + (float(rock_TH232[i,1])*activity232_pmt*mass_rock)
for i in range(8):
    rate = rate + (float(rock_fn[i,1])*activityFN*mass_rock)
 
# VETO RATE
rate = rate + (float(veto_40K[0,1])*activity40_veto*mass_veto)
for i in range(5):
    rate = rate + (float(veto_U238[i,1])*activity238_veto*mass_veto)
for i in range(4):
    rate = rate + (float(veto_TH232[i,1])*activity232_veto*mass_veto)

# WATER RATE
#rate = rate + (float(water_40K[0,1])*activity40*mass_water)
for i in range(4):
    rate = rate + (float(water_222Rn[i,1])*activity222*vol)
for i in range(3):
    rate = rate + (float(water_ibd[i,1])*activityIBD)

#CNO RATE
cno_rate = eff_CNO*nKiloTons/1000*FreeProtons*3600*24

daily_rate = rate * 24 * 3600 * 0.05 * 0.0001 # Remove distance and time cuts (0.05 and 0.0001)
print("CNO_Rate:", cno_rate, "per day")
print("Daily Rate:", daily_rate, "per day")

