clear 
close all


%données Cruise condition

%flight conditions
l=10.700;%Km
M=0.8;
Cp_a=1004;
gamma_a=1.4;
r_a=(Cp_a*(gamma_a-1))/gamma_a;%J/Kg-K

%Componente characteristics:
Fan_pr=1.75;
TPR_IPC=1.3;
TPR_HPC=12;
T_t4=1540;
Eff_comb=0.995;
h_25=-43124e3;%J/kg
Cpg=1147;
gamma_g=4/3;


%Losses 
 Inlet_loss = 0.96;
 Fan_is_eff = 0.91;
 IP_comp_is_eff = 0.815;
 Lo_comp = 0.99;
 HPC_is_eff= 0.89;
 Tot_p_loss_cc = 0.97;
 Poly_eff_HPt= 0.85;
 IntT_loss = 0.98;
 Poly_eff_LPt= 0.85;
 Mec_eff_HPshaft = 0.99;
 Mec_eff_LPshaft = 0.99;

%Mass flow and bleed rates :
 mred = 400;% REDUCED MASSFLOW in SECTION 2 Kg/s
 BPR = 5.5;
 Mecpower_HPshaft = 50e3;% W
 Cool_m_HPdist = 0.06; % of overall massflow in section 25, withdrawn at the exit of HP compressor
 Cool_leak_HPT = 0.04; % of overall massflow in section 25, withdrawn at the exit of HP compressor
 Cool_leak_LPT=0;


%CRUISE


%Tables creation
Section={'atm';'1';'2';'13';'16';'18';'21';'24';'25';'3';'31';'4';'41';'44';'45';'5';'6';'8'};
p0=zeros(18,1);
t0=zeros(18,1);
Massflow=zeros(18,1);
tom = table(p0,t0,Massflow);
tom.Properties.RowNames=Section;
tom.Properties.VariableUnits={'Pa','K','Kg/s'};
in=zeros(12,1);
others=table(in);% 1==yes; 0==no;
others.Properties.RowNames={'Thrust Primary Nozzle';'Primary Nozzle Chocked?';'Static Pressure Primary Nozzle';'Primary Nozzle Area';'Thrust Secondary Nozzle';'Secondary Nozzle Chocked?';'Static Pressure Secondary Nozzle';'Secondary Nozzle Area';'Thrust';'sfc';'Propulsive Efficiency';'Spec. Consumption'};



%atmosphere
tom.t0('atm')=288.15-6.5*l;
tom.p0('atm')=101325*(1-0.0225577*l)^(5.25588);

%Air entrance conduit (Point 1)
tom.p0('1')=tom.p0('atm')*(1+(gamma_a-1)*0.5*M^2)^(gamma_a/(gamma_a-1));%isentropique (??)
tom.t0('1')= tom.t0('atm')*(1+(gamma_a-1)*0.5*M^2);
V1=M*sqrt(gamma_a*r_a*tom.t0('2'));

%Flight velocity
Va=V1;

%Point 2
tom.p0('2')=tom.p0('1')*Inlet_loss;
tom.t0('2')=tom.t0('1');

%True massflow
m=(mred*(tom.p0('2')/101325))/sqrt(tom.t0('2')/288.15);
tom.Massflow('2')=m;

%back to point 1
tom.Massflow('1')=m;

%Secodary Flow

%Point 13
tom.Massflow('13')=m/(1+(1/BPR));
tom.p0('13')=Fan_pr*tom.p0('2');
t013_is=tom.t0('2')*(Fan_pr^((gamma_a-1)/gamma_a));
tom.t0('13')=((t013_is-tom.t0('2'))/ Fan_is_eff)+tom.t0('2');

%Point 16
tom.p0('16')=tom.p0('13');
tom.t0('16')=tom.t0('13');
tom.Massflow('16')=tom.Massflow('13');

%Point 18
tom.t0('18')=tom.t0('16');
tom.Massflow('18')=tom.Massflow('16');

%Critical  conditions (M=1)
pcr18=tom.p0('16')*(1-((gamma_a-1)/(gamma_a+1)))^(gamma_a/(gamma_a-1));
tcr18=2*tom.t0('16')/(gamma_a+1);

if (pcr18>=tom.p0('atm'))
	others.in('Secondary Nozzle Chocked?')=1;
    others.in('Static Pressure Secondary Nozzle')=pcr18;
    t18=tcr18;%Static Temperature Secondary Nozzle
    v18=sqrt(gamma_a*r_a*t18);%==a18
    tom.p0('18')=pcr18*(1+(v18^2/(2*Cp_a*t18)))^(gamma_a/(gamma_a-1));% p0('18')==p0('16')
    
	%Secondary Nozzle Area
	others.in('Secondary Nozzle Area')=(tom.Massflow('18')*r_a*t18)/(pcr18*v18);

	%Thrust Secondary Nozzle: chocked so addition of the term A(p-patm)
	others.in('Thrust Secondary Nozzle')=(tom.Massflow('18')*v18 - m*Va)+(others.in('Secondary Nozzle Area')*(pcr18-tom.p0('atm')));

else
    others.in('Secondary Nozzle Chocked?')=0;
    others.in('Static Pressure Secondary Nozzle')=tom.p0('atm');
    tom.p0('18')=tom.p0('1');%?
end



%Primary Flow

%Point 21
tom.Massflow('21')=m-tom.Massflow('13');
tom.p0('21')=Fan_pr*tom.p0('2');
t021_is=tom.t0('2')*(Fan_pr^((gamma_a-1)/gamma_a));
tom.t0('21')=((t013_is-tom.t0('2'))/ Fan_is_eff)+tom.t0('2');

%Point 24
tom.p0('24')=TPR_IPC*tom.p0('21');
t24is=tom.t0('21')*(TPR_IPC^((gamma_a-1)/gamma_a));
tom.t0('24')=tom.t0('21')+((t24is-tom.t0('21'))/IP_comp_is_eff);
tom.Massflow('24')=tom.Massflow('21');

%Point 25 HP compressor inlet
tom.p0('25')=Lo_comp*tom.p0('24');
tom.t0('25')=tom.t0('24'); %no work, no heat, dh0=0
tom.Massflow('25')=tom.Massflow('24');

%Point 3 HP compressor outlet
tom.p0('3')=TPR_HPC*tom.p0('25');
t03is=tom.t0('25')*(TPR_HPC^((gamma_a-1)/gamma_a));
tom.t0('3')=tom.t0('25')+((t03is-tom.t0('25'))/HPC_is_eff);
tom.Massflow('3')=tom.Massflow('25');

%rq: use gamma_g and Cpg
%Point 31 Combustion Chamber inlet
tom.p0('31')=tom.p0('3');
tom.t0('31')=tom.t0('3');
tom.Massflow('31')=tom.Massflow('3')-((Cool_m_HPdist+Cool_leak_HPT+Cool_leak_LPT)*tom.Massflow('25'));

%Point 4 Combustion Chamber outlet 
tom.p0('4')=Tot_p_loss_cc*tom.p0('31');%Pressure loss in the CC
tom.t0('4')=T_t4;

%addition of fuel
fth=((Cp_a*(tom.t0('31')-298.15))+(Cpg*(298.15-tom.t0('4'))))/(h_25+(Cpg*(tom.t0('4')-298.15))); %theoric necessary quantity of fuel
f=fth/Eff_comb; %actual quantity of fuel/ quantity of air
tom.Massflow('4')=tom.Massflow('31')*(1+f);

% %Point 41 distributor outlet/rotor inlet - first stage of the turbine
 % stator --> no work, no loss --> total pressure conserved
 tom.p0('41')=tom.p0('4');
 %introduction of cool massflow
 mcool=Cool_m_HPdist*tom.Massflow('25');
 tom.Massflow('41')=tom.Massflow('4')+mcool;
 %Heat diffusion by the cool massflow
 tom.t0('41')=(Cpg*tom.Massflow('4')*tom.t0('4')+Cp_a*mcool*tom.t0('25'))/(Cpg*tom.Massflow('4')+Cp_a*mcool);


%PAS CORRECT!!
% %Point 44 - outlet of HP turbine after restitution of the cooling air(A17) 
% tom.Massflow('44')=tom.Massflow('41')+Cool_leak_HPT*tom.Massflow('25');
% dwu44=(Mecpower_HPshaft/tom.Massflow('44'))*Mec_eff_HPshaft;
% tom.t0('44')=tom.t0('41')-(dwu44/Cpg);
% exp=(1/Poly_eff_HPt)*(gamma_g/(gamma_g-1));
% tom.p0('44')=tom.p0('41')*((tom.t0('44')/tom.t0('41'))^exp);
% 
% %Point 45 - inlet Low Pressure Turbine
% tom.p0('45')=IntT_loss*tom.p0('44');
% %no work, no heat
% tom.t0('45')=tom.t0('44');
% tom.Massflow('45')=tom.Massflow('44')+(Cool_leak_LPT*tom.Massflow('25'));

%CORRECT




tom
others









