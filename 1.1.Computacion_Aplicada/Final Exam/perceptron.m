%% NEURON / PERCEPTRON
% let's implement a neuron, which does the following:
%     step 1 : sum of the weighted inputs
%         a) for every input, multiply that input by its weight
%         b) sum all the weighted inputs
%     step 2 : activation function
%         a) compute the output of the perceptron based on that sum passed
%            through an activation function
classdef perceptron
    properties

        Inputs;
        weights;
        target;
        Output;
        
    end
    methods
        
        function obj = perceptron(Inputs,weights,target)
            % class constructor
            if(nargin > 0)
                obj.Inputs = Inputs;
                obj.weights = weights;
                obj.target = target;
                obj.Output = obj.guess(); % save the preliminary guess
            end
        end
        
        function guessedObj = guess(obj)
            nInputs = size(obj.Inputs,2); % get the number of inputs
            % compute the weighted sum
            sum = 0;
            for k=1:nInputs
                sum = sum + obj.Inputs(k)*obj.weights(k);
            end
%             disp(sum); % debugging purposes
            % let the activation function to make an output guess
            guess = activationFunct(sum);
            
            % modify the Output property
            obj.Output = guess;
            % return the modified object
            guessedObj = obj;
        end
        
        function trainedObj = train(obj)
            nInputs = size(obj.Inputs,2); % get the number of inputs
            guess = obj.Output;
%             disp(target) % debugging purposes
%             disp(guess) % debugging purposes
            error = obj.target - guess;
            % value used to avoid tunning overshoots
            learningRate = 0.00001;
            train = obj.weights;
            
            % Tune all the weights
            for k=1:nInputs
                train(1,k) = ...
                    train(1,k) + error*obj.Inputs(k)*learningRate;
%                 disp(obj.Inputs(k)) % debugging purposes
            end
            
            % modify the weights property
            obj.weights = train;
            % return the modified object
            trainedObj = obj;

        end
    end
end

%% The Activation Function
function res = activationFunct(sum)
    % implement a linear trendline to make a guess
    res = 0.0074*sum + 7.4145;
end
