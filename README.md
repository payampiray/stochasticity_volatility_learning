A model for learning based on the joint estimation of stochasticity and volatility 

## reference
please cite this paper if you use this code:
Piray P and Daw ND, 'A model for learning based on the joint estimation of stochasticity and volatility', 2021, Nature Communications.

## description of the models
This work addresses the problem of learning in noisy environments, in which the agent must draw inferences (e.g., about true reward rates) from observations (individual reward amounts) that are corrupted by two distinct sources of noise: process noise or volatility and observation noise or stochasticity. Volatility captures the speed by which the true value being estimated changes from trial to trial (modeled as Gaussian diffusion); stochasticity describes additional measurement noise in the observation of each outcome around its true value (modeled as Gaussian noise on each trial). The celebrated Kalman filter makes inference based on known value for both stochasticity and volatility, in which volatility and stochasticity have opposite effects on the learning rate (i.e. Kalman gain): whereas volatility increases the learning rate, stochasticity decreases the learning rate.

The learning models implemented here generalize the Kalman filter by also learning both stochasticity and volatility based on observations.
An important point is that inferences about volatility and stochasticity are mutually interdependent. But the details of the interdependence are themselves informative. From the learner’s perspective, a challenging problem is to distinguish volatility from stochasticity when both are unknown, because both of them increase the noisiness of observations. Disentangling their respective contributions requires trading off two opposing explanations for the pattern of observations, a process known in Bayesian probability theory as explaining away. This insight results in two lesion models: a stochasticity lesion model that tends to misidentify stochasticity as volatility and inappropriately increases learning rates; and a volatility lesion model that tends to misidentify volatility as stochasticity and inappropriately decreases learning rates.

## MATLAB implementation
This repository contains MATLAB code, which reproduces all figures of the reference paper.

## Python implementation
A python implementation of the models and examples are available at
https://github.com/payampiray/health-lesion-stovol

## updates:
Nov 17, 2021: No requirement of particle filter was implemented, which makes the repository independent of the MATLAB's Control System Toolbox.
Nov 17, 2021: small bug in particle filter was fixed (with no effect on results)
Nov 17, 2021: a notational problem was fixed in sim_lesioned.m
