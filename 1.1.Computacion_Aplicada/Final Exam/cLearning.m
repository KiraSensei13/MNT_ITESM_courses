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
% *     In submitting the solution to this final exam, We Bruno González
% *     Soria and Antonio Osamu Katagiri Tanaka affirm our awareness of the
% *     standards of the Tecnológico de Monterrey Ethics Code.
% *
% *     Thanks to The Coding Train: https://www.youtube.com/watch?v=XJ7HLz9
% *     VYz0&list=PLRqwX-V7Uu6Y7MdSCaIfsxc561QI0U0Tb&index=1
% *
% * START DATE :
% *     02 May 2019
% ************************************************************************

close all, clear all, clc, format compact
rng(31416)

%% * Problem 3: LEARNING **************************************************
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

%% ************************************************************************
% a) Train a neural network using the information of these 4077.
%     1) This data must be divided into training, testing, and possibly
%        validation examples. Explain your decision when choosing these
%        data sets.
%
%        The neural netowrk (AKA brain) is trained through a superviced
%        learning algorithm. Therefore the 4077 entries in shellfish.xlsx
%        are used for both purposes, training and testing. The reason
%        behind this desition is to consider the most amount of entries as
%        possible, therfore to train the brain with all the known data.
%
%     2) Explain any pre-processing done to the data.
%
%        Two steps were performed to pre-process the data: 1) convert the
%        imported cell-matrix data types into ordinary matrices, and
%        2) convert the categorical variables into integer values. These
%        convertions are done to ease the algorithm computation.

% load shellfish data
ssds = spreadsheetDatastore('./shellfish.xlsx');

% store the 1st sheet - DATA **********************************************
ssds.Sheets = 1;
data = read(ssds);
data_arr = zeros(height(data),width(data));

% convert cell matrix to ordinary matrix
for k=1:width(data)
   %data_varNames = data.Properties.VariableNames(k); %debuggin purposes
   %data_varNames = cell2mat(data_varNames); %debuggin purposes
    table_col = table2array(data(1:height(data),k));
    
    if isa(table_col,'cell')
        % convert categorical values to integers:
        % F = 70, M = 77, I = 73
        data_arr(1:height(data),k) = cell2mat(table_col);
    else % keep the values as they are
        data_arr(1:height(data),k) = table_col;
    end
end

% data_arr shall be used for the neural training

% let's train the brain/perceptron with the known training data in
% data_arr - Through a SUPERVISED LEARNING ALGORITHM
%    1) Provide the perceptron with inputs for which there is a known
%       answer
%    2) Ask the perceptron to guess an answer
%    3) Compute the error (Did it get the answer right or wrong?)
%    4) Adjust all the weights according to the error
%    5) Return to step 1) and repeat!

% initialize the weights randomly
nInputs = width(data)-1; % number of inputs
weights = zeros(1,nInputs);
for k=1:nInputs
    weights(1,k) = randi([-10 10]); % random numbers between -10 and 10
end
% from previous trainings, the weights are estimated to be the following:
weights = [-7.5340 5.6109 5.3689 6.3149 2.7319 -10.4802 7.2560 0.1472];

% let's create some variables to see how well brain is being trained
guess = zeros(1,length(data_arr)); % to store what the perceptron guesses
known = zeros(1,length(data_arr)); % to store the correct answers

n = 25; % let's train the brain n-times
for i=1:n
    for k=1:length(data_arr)
        % let's feed our brain/perceptron some intputs
        inputs = data_arr(k,1:width(data)-1);
        % set the known target (correct answer) to compute the error
        target = data_arr(k,width(data));
        % let's create a brain/neuron
        brain = perceptron(inputs,weights,target);
        % let's ask brain for a guess
        brain = brain.guess;
        % let's train the brain according to the previous guess
        brain = brain.train;
        % update the weights according to the training
        weights = brain.weights;
        % populate the tracking variables
        guess(1,k) = brain.Output;
        known(1,k) = brain.target;
    end
end

% Guesses are rounded to the nearest integer. Age is od type integer, as
% defined in shellfish.xlsx data description
guess = round(guess);

%% ************************************************************************
% b) Using your trained neural network, determine the age of the 100
%    individuals in sheet predict. Write the results to as a column of an
%    Excel worksheet.
%
%    The predicted data (ages) are stored in predictedAge

% store the 2nd sheet - PREDICT *******************************************
ssds.Sheets = 2;
predict = read(ssds);
predict_arr = zeros(height(predict),width(predict));

% convert cell matrix to ordinary matrix
for k=1:width(predict)
   %predict_varNames = predict.Properties.VariableNames(k); %debuggin purposes
   %predict_varNames = cell2mat(predict_varNames); %debuggin purposes
    table_col = table2array(predict(1:height(predict),k));

    if isa(table_col,'cell')
        % convert categorical values to integers:
        % F = 70, M = 77, I = 73
        predict_arr(1:height(predict),k) = cell2mat(table_col);
    else % keep the values as they are
        predict_arr(1:height(predict),k) = table_col;
    end
end

% data from predict_arr shall be used to predict the Age

% let's create a variable to store the predictions
predictedAge = zeros(length(predict_arr),1);

% let's guess predict_arr
for k=1:length(predict_arr)
    % let's feed our brain/perceptron some intputs to get a guess.
    inputs = predict_arr(k,1:size(predict_arr,2));
    % set the known target to 0 ... this is not used in this step
    target = 0;
    % let's create a brain
    brain = perceptron(inputs,weights,target);
    % let's ask brain for a guess
    brain = brain.guess;
    % populate the tracking variable
    predictedAge(k,1) = brain.Output;
end

% Predictions are rounded to the nearest integer. Age is od type integer,
% as defined in shellfish.xlsx data description
predictedAge = round(predictedAge);

% Save predictions into a CSV file
tbl = [predict array2table(predictedAge)];
writetable(tbl,'./cLearning.csv','WriteRowNames',true,'Delimiter',',');
disp(tbl)

%% ************************************************************************
% c) Give an estimate of the expected error of your neural network on new
%    data. Explain your answer.
%
%    The espected relative error is around 13% (as calculated with
%    mean(error)). We expect a single neural network to get things right
%    every time, however that's not how the brain works. The firsst neural
%    network gets things wrong often. If we trained an artificial neural
%    network to reduce the error, we would need to slow the learning rate
%    early (learningRate = 0.00001), and halt long before over-fitting
%    (n = 25). However the 1st network would still get many answers wrong.
%    Then, we would need a second neural network, that receives the 1st
%    network mistakes to drop the error significantly.

% let's create more variables to see how well brain is trained
error = abs(guess - known)./known.*100; % error shall be close to zero ...
x = 1:length(error);
% scatter the points with blue filled dots
scatter(x,error,3,[0 0.4470 0.7410],'filled');
title("Scatter Plot of the Training Relative-Error with "+n+" Trains");
xlabel('Guesses');
xlim([0 length(error)]);
ylabel('relativeError = (Guess - Target)/Target');
ylim([-15 115]);
txt = strcat('mean(relativeError) = ', num2str(mean(error)), '%');
text(0,max(error)-2,txt) % print the error's mean into the chart
