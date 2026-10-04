clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. 
%  For the final assignment take the ones from the specific case you are supposed to do.                  
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=102.78;
cFuel='CH4';        
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
% Make a choice for which type of solution method you want to use.

%% [2-3] Compressor

sPart = 'Compressor';

% Step 1: Determine the compressor outlet pressure.
% For Group 96 the prescribed pressure ratio P3/P2 is 7.
P3 = P2 * P3overP2;

% Step 2: Compression from state 2 to state 3 is isentropic.
% Therefore S3 = S2.
%
% Using:
% s3thermal - s2thermal = Rg*ln(P3/P2)
%
% determine the required thermal entropy contribution at state 3.
s3thermal = s2thermal + Rg*log(P3/P2);

% Step 3: Find the compressor outlet temperature corresponding
% to this thermal entropy using the NASA air-property data.
T3 = interp1(sair_a,TR,s3thermal);

% Step 4: Calculate the enthalpy of each species at T3.
for i = 1:NSp
    hi3(i) = HNasa(T3,SpS(i));
    si3(i) = SNasa(T3,SpS(i));
end

% Calculate the enthalpy of the air mixture at compressor outlet.
h3 = Yair*hi3';

% Calculate total specific entropy as a consistency check.
s3thermal_check = Yair*si3';
S3 = s3thermal_check - Rg*log(P3/Pref);

% Velocity is neglected through the compressor in this model.
v3 = 0
% Print compressor results
fprintf('\n%14s\n',sPart);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,2,3);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T2,T3);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P2/kPa,P3/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v2,v3);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h2/kJ,h3/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S2/kJ,S3/kJ);

% Compressor specific work
wc = h3-h2;
fprintf('%8s| %9.2f            [kJ/kg]\n','wc',wc/kJ);
%% [3-4] Combustor
sPart = 'Combustor';
% Step 1: Determine the mass flow rates of air and fuel.
% AF is defined as the air-fuel mass ratio.
mdot_fuel = mfurate;
mdot_air = AF*mdot_fuel;
mdot4 = mdot_air + mdot_fuel;

% Step 2: Determine the composition before combustion.
% State 3 contains the incoming air together with the added fuel.
Y3 = (mdot_air*Yair + mdot_fuel*Yfuel)/mdot4;

% Step 3: Calculate the combustion products.
% For methane complete combustion is:
% CH4 + 2 O2 -> CO2 + 2 H2O
%
% Take 1 kg of fuel as basis.
mfuel = 1;
mair = AF*mfuel;

% Convert the incoming masses to number of kmoles.
nfuel = mfuel/Mi(1);
nO2 = mair*Yair(2)/Mi(2);
nN2 = mair*Yair(5)/Mi(5);

% Complete combustion consumes 2 moles O2 for every mole CH4.
nO2_used = 2*nfuel;
nO2_left = nO2 - nO2_used;

% Product composition in the same species order:
% [CH4 O2 CO2 H2O N2]
n4 = [0 nO2_left nfuel 2*nfuel nN2];

% Convert product mole amounts back to mass fractions.
m4 = n4.*Mi;
Y4 = m4/sum(m4);

% Step 4: Calculate the stoichiometric air-fuel ratio
% and equivalence ratio.
AFst = (2*Mi(2)/Mi(1))/Yair(2);
phi = AFst/AF;

% Step 5: Determine gas constants before and after combustion.
% Rg changes because the mixture composition changes.
M3 = 1/sum(Y3./Mi);
M4 = 1/sum(Y4./Mi);

Rg3 = Runiv/M3;
Rg4 = Runiv/M4;
% Print mixture composition
fprintf('\n-------------------------------------\n');
fprintf('Combustor composition [3-4]\n');
fprintf('-------------------------------------\n');
fprintf('AF = %.2f\n',AF);
fprintf('Equivalence ratio = %.4f\n',phi);

fprintf('\nMass fractions before combustion:\n');
fprintf('Fuel = %.6f\n',Y3(1));
fprintf('O2   = %.6f\n',Y3(2));
fprintf('CO2  = %.6f\n',Y3(3));
fprintf('H2O  = %.6f\n',Y3(4));
fprintf('N2   = %.6f\n',Y3(5));
fprintf('Rg   = %.2f J/kg/K\n',Rg3);

fprintf('\nMass fractions after combustion:\n');
fprintf('Fuel = %.6f\n',Y4(1));
fprintf('O2   = %.6f\n',Y4(2));
fprintf('CO2  = %.6f\n',Y4(3));
fprintf('H2O  = %.6f\n',Y4(4));
fprintf('N2   = %.6f\n',Y4(5));
fprintf('Rg   = %.2f J/kg/K\n',Rg4);
fprintf('-------------------------------------\n');
% Step 6: Assume no pressure loss through the combustor.
P4 = P3;

% Step 7: Determine the enthalpy of the fuel entering the combustor.
% Fuel enters at the reference temperature.
hfuel = HNasa(Tref,SpS(1));

% Step 8: Apply energy conservation over the combustor.
% Energy entering with the air and fuel equals the energy
% leaving with the combustion products.
H4 = mdot_air*h3 + mdot_fuel*hfuel;
h4 = H4/mdot4;

% Step 9: Calculate the enthalpy of the product mixture
% over the complete temperature range.
hprod_a = zeros(1,NTR);

for i = 1:NSp
    hprod_a = hprod_a + Y4(i)*hia(:,i)';
end

% Step 10: Find the combustor outlet temperature.
% Interpolate the NASA product enthalpy data to find T4.
T4 = interp1(hprod_a,TR,h4);

% Step 11: Calculate the properties at state 4.
for i = 1:NSp
    hi4(i) = HNasa(T4,SpS(i));
    si4(i) = SNasa(T4,SpS(i));
end

% Check the calculated product enthalpy.
h4check = Y4*hi4';

% Calculate thermal and total entropy at state 4.
s4thermal = Y4*si4';
S4 = s4thermal - Rg4*log(P4/Pref);

% Velocity is neglected through the combustor.
v4 = 0;

% Print combustor results.
fprintf('\n%14s\n',sPart);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,3,4);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T3,T4);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P3/kPa,P4/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v3,v4);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h3/kJ,h4/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S3/kJ,S4/kJ);

%% [4-5] Turbine

sPart = 'Turbine';

% Step 1: Compressor power requirement.
% wc is the compressor specific work [J/kg] and wc is h3 - h2.
% Only air flows through the compressor.
Wcomp = mdot_air * wc;

% Step 2: Determine the mass flow through the turbine.
% After the combustor, since there is added fuel,both air and fuel products pass through the turbine.
mdot_turbine = mdot_air + mdot_fuel;

% Step 3: Apply shaft-work balance.
% The power supplied by turbine is the power required by the compressor.
% Gas loses energy at turbine so h5<h4.
%
% Wturbine = Wcomp
%
% mdot_turbine*(h4-h5) = Wcomp
%
% Therefore:
% h5 = h4 - Wcomp/mdot_turbine

h5 = h4 - Wcomp/mdot_turbine;

% Step 4: Calculate the enthalpy of the product mixture
% over the temperature range.
% The turbine contains the same product mixture as state 4.

hprod_a = zeros(1,NTR);

for i = 1:NSp
    hprod_a = hprod_a + Y4(i)*hia(:,i)';
end

% Step 5: Find turbine outlet temperature from h5.
T5 = interp1(hprod_a,TR,h5);

% Step 6: Calculate entropy properties at state 5.
for i = 1:NSp
    hi5(i) = HNasa(T5,SpS(i));
    si5(i) = SNasa(T5,SpS(i));
end

% Enthalpy check.
h5check = Y4*hi5';

% Thermal entropy of product mixture at state 5.
s5thermal = Y4*si5';

% Step 7: Turbine is assumed isentropic.
% Therefore:
%
% S5 = S4
%
% S = s_thermal - Rg*ln(P/Pref)
%
% This gives:
%
% ln(P5/P4) = (s5thermal-s4thermal)/Rg4

lnP5P4 = (s5thermal-s4thermal)/Rg4;

P5 = P4*exp(lnP5P4);

% Total entropy at state 5, used as a consistency check.
S5 = s5thermal - Rg4*log(P5/Pref);

% Velocity is neglected through the turbine.
v5 = 0;

% Print turbine results.
fprintf('\n%14s\n',sPart);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,4,5);
fprintf('-------------------------------------\n');

fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T4,T5);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P4/kPa,P5/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v4,v5);

fprintf('---  H/S    -------------------------\n');

fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h4/kJ,h5/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S4/kJ,S5/kJ);

fprintf('-------------------------------------\n');

fprintf('%8s| %9.2f            [kW]\n','Wcomp',Wcomp/kJ);

