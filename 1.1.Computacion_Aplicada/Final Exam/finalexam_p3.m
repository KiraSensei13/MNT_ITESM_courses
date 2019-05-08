%% ************************************************************************
% * AUTHOR(S) :
% *     Bruno González Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW03.m
% *
% * DESCRIPTION :
% *     Computación Aplicada (Ene 19 Gpo 1)
% *     Final Exam
% *
% * NOTES :
% *     In submitting the solution to this final exam, I (we) your name(s)
% *     affirm my (our) awareness of the standards of the Tecnológico de
% *     Monterrey Ethics Code.
% *
% * START DATE :
% *     02 May 2019
% ************************************************************************

close all, clear all, clc, format compact

%% ************************************************************************
% Problem 3: LEARNING
% The age of a specific species of shellfish is related to several physical
% characteristics. The sheet data of the Excel file shellfish.xlsx contains
% data of 4077 individuals, and their ages

% shellfish.xlsx data description:
% ---------------+------------+-------+----------------------------
% Name           | Data Type  | Meas. | Description
% ---------------+------------+-------+----------------------------
% Sex            | nominal    |       | M, F, and I (infant)
% Length         | continuous | mm    | Longest shell measurement
% Diameter       | continuous | mm    | perpendicular to length
% Height         | continuous | mm    | with meat in shell
% Whole weight   | continuous | grams | whole shellfish
% Shucked weight | continuous | grams | weight of meat
% Viscera weight | continuous | grams | gut weight (after bleeding)
% Shell weight   | continuous | grams | after being dried
% Age            | integer    |       | years
% ---------------+------------+-------+----------------------------

% load shellfish data
ssds = spreadsheetDatastore('./shellfish.xlsx');
ssds.Sheets = 1;
data = read(ssds);
ssds.Sheets = 2;
predict = read(ssds);

%% ************************************************************************
% a) Train a neural network using the information of these 4077.
%     1) This data must be divided into training, testing, and possibly
%        validation examples. Explain your decision when choosing these
%        data sets.
%     2) Explain any pre-processing done to the data.



%% ************************************************************************
% b) Using your trained neural network, determine the age of the 100
%    individuals in sheet predict. Write the results to as a column of an
%    Excel worksheet.



%% ************************************************************************
% c) Give an estimate of the expected error of your neural network on new
%    data. Explain your answer.



