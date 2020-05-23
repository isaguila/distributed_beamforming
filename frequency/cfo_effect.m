close all;
clear;
clc;

fc = 5e9;
ppm = 5;
df = fc*1e-6*ppm;

phi_60 = pi/3;
phi_30 = pi/6;

t_60 = phi_60/(2*pi*df)
t_30 = phi_30/(2*pi*df)