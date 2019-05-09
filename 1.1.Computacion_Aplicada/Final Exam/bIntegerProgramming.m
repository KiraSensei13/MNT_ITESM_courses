%% ************************************************************************
% * AUTHOR(S) :
% *     Bruno Gonz�lez Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW03.m
% *
% * DESCRIPTION :
% *     Computaci�n Aplicada (Ene 19 Gpo 1)
% *     Final Exam
% *
% * NOTES :
% *     In submitting the solution to this final exam, We Bruno Gonz�lez
% *     Soria and Antonio Osamu Katagiri Tanaka affirm our awareness of the
% *     standards of the Tecnol�gico de Monterrey Ethics Code.
% *
% * START DATE :
% *     02 May 2019
%% ************************************************************************
% This script should start with the command rng(31416), and should not
% contain any other call that initializes the state of the random number
% generator. 

rng(31416)

%% ************************************************************************
% Problem : INTEGER PROGRAMMING
% 
% An airline company is considering the purchase of new long-, medium-, and
% short-range jet passenger airplanes. The purchase price is $33.5M for
% each long-range plane, $25M for each medium-range plane, and $17.5M for
% each short-range plane. The board of directors has authorized a maximum
% of $750M for these purchases. Regardless of which planes are purchased,
% air travel of all distances is expected to be sufficiently large enough
% so that these planes would be utilized at essentially maximum capacity.
% It is estimated that the net annual profit (after subtracting capital
% recovery costs) would be $2.1M per long-range plane, $1.5M per
% medium-range plane, and $1.15M per short-range plane.
% 
% Enough trained pilots are available to the company to crew 30 new
% airplanes. If only short-range planes were purchased, facilities would be
% able to handle 40 new planes. However, each medium-plane is equivalent to
% ? ? ? short-range planes, and each long-range plane is equivalent to ? ?
% ? short-range planes in terms of their use of maintenance facilities.
% Using the preceding data, management wishes to know how many planes of
% each type should be purchased to maximize profit.

%% ************************************************************************
% a) Formulate the problem as an integer programming problem.
% b) Use intlinprog to find the solution (number of planes of each type and maximum profit)

